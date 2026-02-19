from __future__ import annotations

from io import BytesIO
import json
from types import SimpleNamespace
from typing import Any
from urllib import error as urlerror

import pytest

from services import runtime_nodes
from services.runtime_nodes import RuntimeNodeStore


def _set_required_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CONVEX_URL", "https://convex.example/")
    monkeypatch.setenv("CONVEX_DEPLOY_KEY", "deploy-key")


class _FakeResponse:
    def __init__(self, body: str) -> None:
        self._body = body

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None

    def read(self) -> bytes:
        return self._body.encode("utf-8")


def test_constructor_requires_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("CONVEX_URL", raising=False)
    monkeypatch.delenv("CONVEX_DEPLOY_KEY", raising=False)

    with pytest.raises(RuntimeError, match="CONVEX_URL and CONVEX_DEPLOY_KEY are required"):
        RuntimeNodeStore()


def test_constructor_normalizes_url_and_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    monkeypatch.setenv("RUNTIME_NODE_STORE_TIMEOUT_SECONDS", "0")

    store = RuntimeNodeStore()

    assert store._convex_url == "https://convex.example"
    assert store._timeout_seconds == 1


def test_constructor_timeout_invalid_falls_back_to_default(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    monkeypatch.setenv("RUNTIME_NODE_STORE_TIMEOUT_SECONDS", "not-a-number")

    store = RuntimeNodeStore()

    assert store._timeout_seconds == 10


def test_request_success_builds_body_and_extracts_value(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    monkeypatch.setenv("RUNTIME_NODE_STORE_TIMEOUT_SECONDS", "7")
    store = RuntimeNodeStore()

    captured: dict[str, Any] = {}

    def _fake_urlopen(req, timeout: int):
        captured["url"] = req.full_url
        captured["timeout"] = timeout
        captured["auth"] = req.get_header("Authorization")
        captured["ctype"] = req.get_header("Content-type")
        captured["body"] = json.loads(req.data.decode("utf-8"))
        return _FakeResponse('{"value":{"ok":true}}')

    monkeypatch.setattr(runtime_nodes.urlrequest, "urlopen", _fake_urlopen)

    value = store._request("mutation", "runtimeTargets:addRuntimeTarget", {"a": 1})

    assert value == {"ok": True}
    assert captured == {
        "url": "https://convex.example/api/mutation",
        "timeout": 7,
        "auth": "Convex deploy-key",
        "ctype": "application/json",
        "body": {"path": "runtimeTargets:addRuntimeTarget", "args": {"a": 1}},
    }


def test_request_handles_http_error(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    def _raise_http_error(req, timeout: int):
        raise urlerror.HTTPError(req.full_url, 403, "forbidden", hdrs=None, fp=BytesIO(b"denied"))

    monkeypatch.setattr(runtime_nodes.urlrequest, "urlopen", _raise_http_error)

    with pytest.raises(RuntimeError, match=r"convex query failed: HTTP 403 denied"):
        store._request("query", "runtimeTargets:listTargets", {})


def test_request_handles_url_error(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    def _raise_url_error(req, timeout: int):
        raise urlerror.URLError("down")

    monkeypatch.setattr(runtime_nodes.urlrequest, "urlopen", _raise_url_error)

    with pytest.raises(RuntimeError, match=r"convex mutation failed"):
        store._request("mutation", "runtimeTargets:addRuntimeTarget", {})


def test_request_handles_empty_and_non_dict_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    monkeypatch.setattr(runtime_nodes.urlrequest, "urlopen", lambda req, timeout: _FakeResponse("   "))
    assert store._request("query", "x", {}) == {}

    monkeypatch.setattr(runtime_nodes.urlrequest, "urlopen", lambda req, timeout: _FakeResponse("[1,2,3]"))
    assert store._request("query", "x", {}) == {}


def test_normalize_queue_name() -> None:
    assert runtime_nodes._normalize_queue_name("  zpe  ") == "zpe"
    with pytest.raises(ValueError, match="queue_name must be a non-empty string"):
        runtime_nodes._normalize_queue_name("   ")


def test_create_join_token_normalizes_queue_and_sets_label(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()
    monkeypatch.setattr(runtime_nodes.secrets, "token_urlsafe", lambda n: "token-123")
    monkeypatch.setattr(runtime_nodes, "_now_iso", lambda: "2026-01-01T00:00:00+00:00")
    monkeypatch.setattr(runtime_nodes, "_expires_iso", lambda ttl: "2026-01-01T00:10:00+00:00")

    captured: dict[str, Any] = {}

    def _fake_request(method: str, path: str, args: dict[str, Any]) -> dict[str, Any]:
        captured["method"] = method
        captured["path"] = path
        captured["args"] = args
        return {}

    monkeypatch.setattr(store, "_request", _fake_request)

    token = store.create_join_token(
        tenant_id="tenant-1",
        owner_user_id="user-1",
        queue_name="  qe-main  ",
        ttl_seconds=600,
        node_name_hint="node-A",
    )

    assert token.queue_name == "qe-main"
    assert token.label == "qe-main"
    assert token.token == "token-123"
    assert captured["method"] == "mutation"
    assert captured["path"] == "runtimeTargets:createJoinToken"
    assert captured["args"]["queue_name"] == "qe-main"
    assert captured["args"]["label"] == "qe-main"


def test_create_join_token_validates_ttl(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    with pytest.raises(ValueError, match="ttl_seconds must be >= 1"):
        store.create_join_token(
            tenant_id="tenant-1",
            owner_user_id="user-1",
            queue_name="zpe",
            ttl_seconds=0,
        )


def test_consume_join_token_error_mapping(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    monkeypatch.setattr(store, "_request", lambda method, path, args: [])
    with pytest.raises(KeyError, match="join token not found"):
        store.consume_join_token("abc")

    monkeypatch.setattr(store, "_request", lambda method, path, args: {"ok": False, "reason": "expired"})
    with pytest.raises(KeyError, match="expired"):
        store.consume_join_token("abc")

    monkeypatch.setattr(store, "_request", lambda method, path, args: {"ok": True, "join_token": "bad"})
    with pytest.raises(RuntimeError, match="invalid join token payload"):
        store.consume_join_token("abc")


def test_consume_join_token_success_parses_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()
    monkeypatch.setattr(runtime_nodes, "_now_iso", lambda: "2026-01-02T00:00:00+00:00")

    monkeypatch.setattr(
        store,
        "_request",
        lambda method, path, args: {
            "ok": True,
            "join_token": {
                "tenant_id": "tenant-1",
                "owner_user_id": "owner-1",
                "queue_name": "zpe",
                "created_at": "2026-01-01T00:00:00+00:00",
                "expires_at": "2026-01-01T01:00:00+00:00",
            },
        },
    )

    token = store.consume_join_token("fallback-token")

    assert token.token == "fallback-token"
    assert token.tenant_id == "tenant-1"
    assert token.owner_user_id == "owner-1"
    assert token.queue_name == "zpe"


def test_add_target_and_set_active_issue_expected_requests(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()
    monkeypatch.setattr(runtime_nodes, "uuid4", lambda: SimpleNamespace(hex="abcd"))
    monkeypatch.setattr(runtime_nodes, "_now_iso", lambda: "2026-01-03T00:00:00+00:00")

    calls: list[tuple[str, str, dict[str, Any]]] = []

    def _fake_request(method: str, path: str, args: dict[str, Any]) -> dict[str, Any]:
        calls.append((method, path, args))
        return {}

    monkeypatch.setattr(store, "_request", _fake_request)

    target = store.add_target(
        tenant_id="tenant-1",
        user_id="user-1",
        queue_name="  zpe-q  ",
        server_id="server-1",
        name="Primary",
        metadata={"cpu": 16},
    )
    store.set_active_target("tenant-1", "user-1", target.target_id)

    assert target.target_id == "qt-abcd"
    assert target.queue_name == "zpe-q"
    assert calls[0][0] == "mutation"
    assert calls[0][1] == "runtimeTargets:addRuntimeTarget"
    assert calls[0][2]["queue_name"] == "zpe-q"
    assert calls[1] == (
        "mutation",
        "runtimeTargets:setActiveTarget",
        {
            "tenant_id": "tenant-1",
            "user_id": "user-1",
            "target_id": "qt-abcd",
            "updated_at": "2026-01-03T00:00:00+00:00",
        },
    )


def test_list_targets_filters_invalid_rows(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    monkeypatch.setattr(
        store,
        "_request",
        lambda method, path, args: [
            "bad",
            {"target_id": "qt-1"},
            {
                "target_id": "qt-2",
                "tenant_id": "tenant-1",
                "user_id": "user-1",
                "queue_name": "zpe",
                "server_id": "server-1",
                "registered_at": "2026-01-04T00:00:00+00:00",
                "name": "Node",
            },
        ],
    )

    targets = store.list_targets("tenant-1", "user-1")

    assert len(targets) == 1
    assert targets[0].target_id == "qt-2"


def test_get_target_and_active_target_parsing(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()

    monkeypatch.setattr(store, "_request", lambda method, path, args: "bad")
    assert store.get_target_by_id("qt-1") is None
    assert store.get_active_target("tenant-1", "user-1") is None

    responses = {
        "runtimeTargets:getTargetById": {
            "tenant_id": "tenant-1",
            "user_id": "user-1",
            "queue_name": "zpe",
            "server_id": "server-1",
            "registered_at": "2026-01-04T00:00:00+00:00",
            "name": "Node",
        },
        "runtimeTargets:getActiveTarget": {
            "target_id": "qt-1",
            "tenant_id": "tenant-1",
            "user_id": "user-1",
            "queue_name": "zpe",
            "server_id": "server-1",
            "registered_at": "2026-01-04T00:00:00+00:00",
        },
    }

    monkeypatch.setattr(store, "_request", lambda method, path, args: responses[path])

    fetched = store.get_target_by_id("qt-1")
    active = store.get_active_target("tenant-1", "user-1")

    assert fetched is not None
    assert fetched.target_id == "qt-1"
    assert fetched.name == "Node"
    assert active is not None
    assert active.target_id == "qt-1"


def test_ensure_target_owner_checks_tenant_and_user(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    store = RuntimeNodeStore()
    target = runtime_nodes.RuntimeTarget(
        target_id="qt-1",
        tenant_id="tenant-1",
        user_id="user-1",
        queue_name="zpe",
        server_id="server-1",
        registered_at="2026-01-04T00:00:00+00:00",
    )

    monkeypatch.setattr(store, "get_target_by_id", lambda target_id: target)
    assert store.ensure_target_owner("tenant-1", "user-1", "qt-1") == target

    with pytest.raises(KeyError, match="target not found"):
        store.ensure_target_owner("tenant-2", "user-1", "qt-1")

    with pytest.raises(KeyError, match="target not found"):
        store.ensure_target_owner("tenant-1", "user-2", "qt-1")


def test_get_runtime_node_store_singleton(monkeypatch: pytest.MonkeyPatch) -> None:
    _set_required_env(monkeypatch)
    monkeypatch.setattr(runtime_nodes, "_runtime_node_store", None)

    first = runtime_nodes.get_runtime_node_store()
    second = runtime_nodes.get_runtime_node_store()

    assert first is second
