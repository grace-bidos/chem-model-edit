from app.api import app
import app.routers.onboarding as onboarding_router

app.include_router(onboarding_router.router)

__all__ = ["app"]
