import httpx
from fastapi.testclient import TestClient

from app.main import create_app
from app.services.pubchem import PubChemCompound, lookup_pubchem_compound


def test_pubchem_lookup_returns_isomeric_smiles(settings, monkeypatch):
    def fake_lookup(name: str) -> PubChemCompound:
        assert name == "Aspirin"
        return PubChemCompound(
            cid=2244,
            molecular_name="Aspirin",
            isomeric_smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        )

    monkeypatch.setattr("app.api.routes.parents.lookup_pubchem_compound", fake_lookup)
    with TestClient(create_app(settings)) as client:
        response = client.get("/api/parents/pubchem", params={"name": "Aspirin"})

    assert response.status_code == 200
    assert response.json() == {
        "cid": 2244,
        "molecular_name": "Aspirin",
        "isomeric_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    }


def test_pubchem_service_accepts_current_smiles_property(monkeypatch):
    response = httpx.Response(
        200,
        request=httpx.Request("GET", "https://pubchem.ncbi.nlm.nih.gov"),
        json={
            "PropertyTable": {
                "Properties": [
                    {"CID": 2244, "Title": "Aspirin", "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O"}
                ]
            }
        },
    )
    monkeypatch.setattr("app.services.pubchem.httpx.get", lambda *args, **kwargs: response)

    compound = lookup_pubchem_compound(" Aspirin ")

    assert compound.cid == 2244
    assert compound.molecular_name == "Aspirin"
    assert compound.isomeric_smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"


def test_validate_parent_returns_canonical_structure_data(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post("/api/parents/validate", json={"smiles": " C(C)O "})

    assert response.status_code == 200
    body = response.json()
    assert body["input_smiles"] == "C(C)O"
    assert body["canonical_smiles"] == "CCO"
    assert body["atom_count"] == 3
    assert body["molecular_formula"] == "C2H6O"
    assert "V2000" in body["molblock"]


def test_validate_parent_preserves_stereochemistry(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/parents/validate", json={"smiles": "N[C@@H](C)C(=O)O"}
        )

    assert response.status_code == 200
    assert "@" in response.json()["canonical_smiles"]


def test_validate_libinvent_requires_and_preserves_attachment_points(settings):
    with TestClient(create_app(settings)) as client:
        valid = client.post(
            "/api/parents/validate",
            json={"smiles": "c1ccc([*:1])cc1", "model_key": "libinvent"},
        )
        missing = client.post(
            "/api/parents/validate",
            json={"smiles": "c1ccccc1", "model_key": "libinvent"},
        )

    assert valid.status_code == 200
    assert "*" in valid.json()["canonical_smiles"]
    assert missing.status_code == 422
    assert missing.json()["error"]["context"]["reason"] == "attachment_point_required"


def test_validate_libinvent_accepts_ketcher_attachment_smiles(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/parents/validate",
            json={
                "smiles": "N#CC[C@@H]%91C1CCCC1.[*:1]%91",
                "model_key": "libinvent",
            },
        )

    assert response.status_code == 200
    assert response.json()["canonical_smiles"] == "N#CC[C@H](C1CCCC1)[*:1]"


def test_validate_linkinvent_accepts_two_pipe_separated_warheads(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/parents/validate",
            json={
                "smiles": "Oc1cncc(*)c1|*c1ccoc1",
                "model_key": "linkinvent",
            },
        )

    assert response.status_code == 200
    body = response.json()
    assert len(body["canonical_smiles"].split("|")) == 2
    assert body["canonical_smiles"].count("*") == 2
    assert body["atom_count"] == 14
    assert "V2000" in body["molblock"]


def test_validate_linkinvent_rejects_wrong_separator_and_attachment_count(settings):
    cases = {
        "Oc1cncc(*)c1.*c1ccoc1": "linkinvent_pipe_required",
        "Oc1cnccc1|*c1ccoc1": "attachment_point_count",
        "Oc1cncc(*)c1": "linkinvent_warhead_count",
    }
    with TestClient(create_app(settings)) as client:
        for smiles, reason in cases.items():
            response = client.post(
                "/api/parents/validate",
                json={"smiles": smiles, "model_key": "linkinvent"},
            )
            assert response.status_code == 422
            assert response.json()["error"]["context"]["reason"] == reason


def test_validate_parent_rejects_invalid_reaction_and_fragments(settings):
    cases = {
        "not a smiles": "parse_failed",
        "CCO.CN": "multiple_fragments",
        "CCO>>CC=O": "reaction_not_supported",
    }
    with TestClient(create_app(settings)) as client:
        for smiles, reason in cases.items():
            response = client.post("/api/parents/validate", json={"smiles": smiles})
            assert response.status_code == 422
            assert response.json()["error"] == {
                "code": "invalid_parent_smiles",
                "message": response.json()["error"]["message"],
                "context": {"reason": reason},
            }


def test_validate_parent_enforces_atom_limit(settings):
    limited = settings.__class__(
        repository_root=settings.repository_root,
        database_path=settings.database_path,
        jobs_dir=settings.jobs_dir,
        parent_max_atoms=2,
    )
    with TestClient(create_app(limited)) as client:
        response = client.post("/api/parents/validate", json={"smiles": "CCO"})

    assert response.status_code == 422
    assert response.json()["error"]["context"]["reason"] == "atom_limit_exceeded"


def test_job_creation_uses_the_same_validation_and_stores_canonical_smiles(settings):
    with TestClient(create_app(settings)) as client:
        invalid = client.post("/api/jobs", json={"parent_smiles": "CCO.CN"})
        valid = client.post("/api/jobs", json={"parent_smiles": "C(C)O"})

    assert invalid.status_code == 422
    assert valid.status_code == 202
    assert valid.json()["parent_smiles"] == "C(C)O"
    assert valid.json()["parent_smiles_canonical"] == "CCO"


def test_linkinvent_job_creation_writes_model_input_with_pipe_separator(settings):
    value = "Oc1cncc(*)c1|*c1ccoc1"
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/jobs",
            json={
                "parent_smiles": value,
                "run_conditions": {"model_key": "linkinvent"},
            },
        )

    assert response.status_code == 202
    job_id = response.json()["id"]
    saved = (settings.jobs_dir / job_id / "input" / "parent.smi").read_text()
    assert saved.count("|") == 1
    assert saved.count("*") == 2
