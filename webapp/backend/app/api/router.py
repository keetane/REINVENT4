from fastapi import APIRouter

from app.api.routes.health import router as health_router
from app.api.routes.jobs import router as jobs_router
from app.api.routes.job_results import router as job_results_router
from app.api.routes.parents import router as parents_router


api_router = APIRouter(prefix="/api")
api_router.include_router(health_router)
api_router.include_router(parents_router)
api_router.include_router(jobs_router)
api_router.include_router(job_results_router)
