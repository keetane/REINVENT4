from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app import __version__
from app.api.router import api_router
from app.api.error_handlers import register_error_handlers
from app.config import Settings, get_settings
from app.db.migrations import apply_migrations


def create_app(settings: Settings | None = None) -> FastAPI:
    resolved_settings = settings or get_settings()

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        app.state.settings.jobs_dir.mkdir(parents=True, exist_ok=True)
        app.state.settings.model_store_dir.mkdir(parents=True, exist_ok=True)
        app.state.settings.qsar_model_store_dir.mkdir(parents=True, exist_ok=True)
        apply_migrations(app.state.settings)
        yield

    application = FastAPI(
        title="REINVENT4 Web API",
        version=__version__,
        lifespan=lifespan,
    )
    application.state.settings = resolved_settings
    application.add_middleware(
        CORSMiddleware,
        allow_origins=list(resolved_settings.cors_origins),
        allow_credentials=False,
        allow_methods=["GET", "POST", "DELETE"],
        allow_headers=["Content-Type", "Accept"],
    )
    application.include_router(api_router)
    register_error_handlers(application)
    return application


app = create_app()
