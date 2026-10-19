from fastapi.testclient import TestClient

from app.main import create_app


def test_upload_model_registers_and_serves_checkpoint(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/models/upload?name=external%20model&model_type=reinvent&base_model_key=reinvent",
            content=b"checkpoint-bytes",
            headers={"Content-Type": "application/octet-stream", "X-Model-Filename": "external.model"},
        )
        assert response.status_code == 201
        model = response.json()
        assert model["source"] == "transfer_learning"
        assert model["name"] == "external_model"

        listing = client.get("/api/models").json()["items"]
        assert any(item["id"] == model["id"] for item in listing)
        artifact = client.get(f"/api/models/{model['id']}/artifacts/model")
        assert artifact.status_code == 200
        assert artifact.content == b"checkpoint-bytes"
        assert "_TL_Reinvent_external_model.model" in artifact.headers["content-disposition"]


def test_upload_model_rejects_incompatible_metadata(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/models/upload?name=bad&model_type=reinvent&base_model_key=mol2mol_scaffold_generic",
            content=b"checkpoint",
        )
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "incompatible_model_metadata"


def test_upload_model_infers_metadata_and_comment_from_filename(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/models/upload",
            content=b"checkpoint-bytes",
            headers={
                "Content-Type": "application/octet-stream",
                "X-Model-Filename": "20261019_RL_Mol2MolHigh_external_agent.chkpt",
            },
        )

    assert response.status_code == 201
    model = response.json()
    assert model["name"] == "external_agent"
    assert model["source"] == "reinforcement_learning"
    assert model["model_type"] == "mol2mol"
    assert model["base_model_key"] == "mol2mol_high"


def test_upload_model_requires_metadata_for_unstructured_filename(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/models/upload",
            content=b"checkpoint-bytes",
            headers={"X-Model-Filename": "external.model"},
        )

    assert response.status_code == 422
    assert response.json()["error"]["code"] == "model_metadata_required"
