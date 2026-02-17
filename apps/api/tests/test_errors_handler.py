from __future__ import annotations

import json

from fastapi import HTTPException

from app.errors import http_exception_handler


def test_http_exception_handler_structured_detail_uses_safe_message() -> None:
    response = http_exception_handler(
        None,  # type: ignore[arg-type]
        HTTPException(status_code=400, detail={"reason": "invalid", "count": 2}),
    )

    body = json.loads(response.body.decode("utf-8"))
    assert response.status_code == 400
    assert body["error"]["message"] == "request failed"
    assert body["error"]["details"] == {"reason": "invalid", "count": 2}


def test_http_exception_handler_string_detail_keeps_message() -> None:
    response = http_exception_handler(
        None,  # type: ignore[arg-type]
        HTTPException(status_code=404, detail="Structure not found"),
    )

    body = json.loads(response.body.decode("utf-8"))
    assert response.status_code == 404
    assert body["error"]["message"] == "Structure not found"
    assert "details" not in body["error"]
