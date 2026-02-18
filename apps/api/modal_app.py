from __future__ import annotations

import os

import modal

from app.api import app as fastapi_app

image = (
    modal.Image.debian_slim(python_version="3.13")
    .pip_install_from_pyproject("apps/api/pyproject.toml")
    .add_local_dir("apps/api", remote_path="/root")
)

app = modal.App("chem-model-edit-api")
runtime_secret_name = os.environ.get(
    "MODAL_RUNTIME_SECRET_NAME",
    "chem-model-edit-runtime-gateway-dev",
)


@app.function(
    image=image,
    timeout=300,
    min_containers=1,
    secrets=[modal.Secret.from_name(runtime_secret_name)],
)
@modal.asgi_app(requires_proxy_auth=True)
def api():
    return fastapi_app
