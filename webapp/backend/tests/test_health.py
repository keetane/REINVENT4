from fastapi.testclient import TestClient

from app.main import create_app


def test_health_initializes_database(settings):
    with TestClient(create_app(settings)) as client:
        response = client.get("/api/health")

    assert response.status_code == 200
    assert response.json()["status"] == "ok"
    assert response.json()["database"] == "ok"
    assert settings.database_path.is_file()


def test_openapi_is_generated(settings):
    with TestClient(create_app(settings)) as client:
        schema = client.get("/openapi.json").json()

    assert schema["info"]["title"] == "REINVENT4 Web API"
    assert "/api/health" in schema["paths"]
    assert "/api/jobs" in schema["paths"]
    assert "/api/jobs/{job_id}/cancel" in schema["paths"]
