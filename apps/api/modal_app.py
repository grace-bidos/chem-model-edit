from __future__ import annotations

import modal

from app.api import app as fastapi_app

image = (
    modal.Image.debian_slim(python_version="3.13")
    .pip_install_from_pyproject("apps/api/pyproject.toml")
    .add_local_dir("apps/api", remote_path="/root")
)

app = modal.App("chem-model-edit-api")


@app.function(
    image=image,
    timeout=300,
    min_containers=1,
)
@modal.asgi_app(requires_proxy_auth=True)
def api():
    return fastapi_app
