from __future__ import annotations

import os

from fastapi import FastAPI, HTTPException
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from redis.exceptions import RedisError

from app.errors import (
    http_exception_handler,
    overflow_error_handler,
    redis_error_handler,
    validation_exception_handler,
    value_error_handler,
)
from app.middleware import add_request_context
import app.routers.health as health_router
import app.routers.onboarding as onboarding_router
import app.routers.runtime as runtime_router
import app.routers.structures as structures_router
import app.routers.supercells as supercells_router
import app.routers.transforms as transforms_router
import app.routers.zpe as zpe_router


def _cors_allow_origins() -> list[str]:
    raw = os.getenv("CORS_ALLOW_ORIGINS")
    if raw:
        return [origin.strip() for origin in raw.split(",") if origin.strip()]
    return [
        "http://localhost:3000",
        "http://localhost:3001",
        "https://chem-model-edit.tadashi240312.workers.dev",
    ]


def _cors_allow_origin_regex() -> str | None:
    raw = os.getenv("CORS_ALLOW_ORIGIN_REGEX")
    if raw is not None:
        raw = raw.strip()
        return raw or None
    return r"^https://[a-z0-9-]+-chem-model-edit\.tadashi240312\.workers\.dev$"


def create_app() -> FastAPI:
    app = FastAPI(
        title="Chem Model API",
        version="0.2.0",
        docs_url="/api/docs",
        redoc_url="/api/redoc",
        openapi_url="/api/openapi.json",
    )

    app.add_middleware(
        CORSMiddleware,
        allow_origins=_cors_allow_origins(),
        allow_origin_regex=_cors_allow_origin_regex(),
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    app.middleware("http")(add_request_context)

    app.add_exception_handler(HTTPException, http_exception_handler)
    app.add_exception_handler(RequestValidationError, validation_exception_handler)
    app.add_exception_handler(ValueError, value_error_handler)
    app.add_exception_handler(OverflowError, overflow_error_handler)
    app.add_exception_handler(RedisError, redis_error_handler)

    app.include_router(health_router.router)
    app.include_router(structures_router.router)
    app.include_router(transforms_router.router)
    app.include_router(supercells_router.router)
    app.include_router(zpe_router.router)
    app.include_router(runtime_router.router)
    app.include_router(onboarding_router.router)

    return app


app = create_app()
