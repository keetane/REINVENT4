from fastapi import APIRouter, Request

from app import __version__
from app.db.connection import connect
from app.schemas.health import HealthResponse


router = APIRouter(tags=["system"])


@router.get("/health", response_model=HealthResponse)
def health(request: Request) -> HealthResponse:
    settings = request.app.state.settings
    with connect(settings) as connection:
        connection.execute("SELECT 1").fetchone()
    return HealthResponse(status="ok", database="ok", version=__version__)

