from fastapi import APIRouter

from app.api.routes.health import router as health_router
from app.api.routes.jobs import router as jobs_router
from app.api.routes.job_results import router as job_results_router
from app.api.routes.parents import router as parents_router
from app.api.routes.models import router as models_router
from app.api.routes.training import router as training_router
from app.api.routes.qsar import router as qsar_router
from app.api.routes.rl import router as rl_router


api_router = APIRouter(prefix="/api")
api_router.include_router(health_router)
api_router.include_router(parents_router)
api_router.include_router(models_router)
api_router.include_router(training_router)
api_router.include_router(qsar_router)
api_router.include_router(rl_router)
api_router.include_router(jobs_router)
api_router.include_router(job_results_router)
