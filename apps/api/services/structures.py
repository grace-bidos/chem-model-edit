from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import import_module
from io import StringIO
import json
import os
from threading import RLock
from typing import Any, Callable, Dict, Optional, cast
from urllib import error as urlerror
from urllib import request as urlrequest
from uuid import uuid4

from ase import Atoms as ASEAtoms
from app.schemas.common import QeParameters, Structure
from services.cif import atoms_to_cif
from services.parse import extract_qe_params, parse_qe_atoms, structure_from_ase

ase_read = cast(Callable[..., Any], getattr(import_module("ase.io"), "read"))


@dataclass
class StoredStructure:
    atoms: ASEAtoms
    source: str
    cif: str
    params: QeParameters | None
    raw_input: str | None
    created_at: datetime
    tenant_id: str | None = None
    workspace_id: str | None = None


class InMemoryStructureStore:
    def __init__(self) -> None:
        self._items: Dict[str, StoredStructure] = {}
        self._lock = RLock()

    def create(
        self,
        atoms: ASEAtoms,
        source: str,
        cif: str,
        params: QeParameters | None,
        raw_input: str | None,
        *,
        tenant_id: str | None = None,
        workspace_id: str | None = None,
    ) -> str:
        now = datetime.now(timezone.utc)
        structure_id = uuid4().hex
        entry = StoredStructure(
            atoms=atoms,
            source=source,
            cif=cif,
            params=params,
            raw_input=raw_input,
            created_at=now,
            tenant_id=tenant_id,
            workspace_id=workspace_id,
        )
        with self._lock:
            self._items[structure_id] = entry
        return structure_id

    def get(
        self,
        structure_id: str,
        *,
        tenant_id: str | None = None,
        workspace_id: str | None = None,
    ) -> Optional[StoredStructure]:
        with self._lock:
            entry = self._items.get(structure_id)
        if entry is None:
            return None
        if tenant_id and entry.tenant_id and tenant_id != entry.tenant_id:
            return None
        if workspace_id and entry.workspace_id and workspace_id != entry.workspace_id:
            return None
        return entry


class ConvexStructureStore:
    def __init__(
        self,
        *,
        convex_url: str,
        convex_deploy_key: str,
        chunk_size: int = 120_000,
        timeout_seconds: int = 10,
    ) -> None:
        self._convex_url = convex_url.rstrip("/")
        self._convex_deploy_key = convex_deploy_key
        self._chunk_size = max(8_192, chunk_size)
        self._timeout_seconds = max(1, timeout_seconds)

    def _request(self, method: str, path: str, args: dict[str, Any]) -> dict[str, Any]:
        body = json.dumps(
            {"path": path, "args": args},
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode("utf-8")
        url = f"{self._convex_url}/api/{method}"
        req = urlrequest.Request(
            url,
            method="POST",
            data=body,
            headers={
                "Authorization": f"Convex {self._convex_deploy_key}",
                "Content-Type": "application/json",
            },
        )
        try:
            with urlrequest.urlopen(req, timeout=self._timeout_seconds) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
        except urlerror.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"convex {method} failed: HTTP {exc.code} {detail}") from exc
        except urlerror.URLError as exc:
            raise RuntimeError(f"convex {method} failed") from exc

        payload_obj: object = json.loads(raw) if raw.strip() else {}
        if not isinstance(payload_obj, dict):
            raise RuntimeError(f"convex {method} returned invalid payload")
        payload = cast(dict[str, Any], payload_obj)
        value = payload.get("value")
        if value is None:
            # Some Convex deployments return raw value for simple responses.
            return payload
        if not isinstance(value, dict):
            return {"value": value}
        return cast(dict[str, Any], value)

    def _query(self, path: str, args: dict[str, Any]) -> dict[str, Any]:
        return self._request("query", path, args)

    def _mutation(self, path: str, args: dict[str, Any]) -> dict[str, Any]:
        return self._request("mutation", path, args)

    def _chunks(self, text: str) -> list[str]:
        if not text:
            return []
        return [
            text[offset : offset + self._chunk_size]
            for offset in range(0, len(text), self._chunk_size)
        ]

    def create(
        self,
        atoms: ASEAtoms,
        source: str,
        cif: str,
        params: QeParameters | None,
        raw_input: str | None,
        *,
        tenant_id: str | None = None,
        workspace_id: str | None = None,
    ) -> str:
        now = datetime.now(timezone.utc).isoformat()
        structure_id = uuid4().hex
        resolved_tenant = tenant_id or "dev"
        resolved_workspace = workspace_id or resolved_tenant
        raw_chunks = self._chunks(raw_input or "")
        structure = structure_from_ase(atoms).model_dump()
        self._mutation(
            "structures:putStructureRecord",
            {
                "tenant_id": resolved_tenant,
                "workspace_id": resolved_workspace,
                "structure_id": structure_id,
                "source": source,
                "cif": cif,
                "structure": structure,
                "params": params.model_dump() if params else None,
                "created_at": now,
                "updated_at": now,
                "raw_input_chunk_count": len(raw_chunks),
            },
        )
        for index, chunk in enumerate(raw_chunks):
            self._mutation(
                "structures:putStructureRawChunk",
                {
                    "tenant_id": resolved_tenant,
                    "workspace_id": resolved_workspace,
                    "structure_id": structure_id,
                    "chunk_index": index,
                    "chunk_text": chunk,
                },
            )
        return structure_id

    def get(
        self,
        structure_id: str,
        *,
        tenant_id: str | None = None,
        workspace_id: str | None = None,
    ) -> Optional[StoredStructure]:
        payload = self._query(
            "structures:getStructureRecordById",
            {"structure_id": structure_id},
        )
        if payload.get("_id") is None:
            return None
        if tenant_id and payload.get("tenant_id") != tenant_id:
            return None
        if workspace_id and payload.get("workspace_id") != workspace_id:
            return None

        entry_tenant = cast(str, payload["tenant_id"])
        entry_workspace = cast(str, payload["workspace_id"])
        chunks_payload = self._query(
            "structures:getStructureRawChunks",
            {
                "tenant_id": entry_tenant,
                "workspace_id": entry_workspace,
                "structure_id": structure_id,
            },
        )

        chunks: list[dict[str, Any]] = []
        if isinstance(chunks_payload, list):
            chunks = cast(list[dict[str, Any]], chunks_payload)
        else:
            value = chunks_payload.get("value")
            if isinstance(value, list):
                chunks = cast(list[dict[str, Any]], value)
        raw_input = "".join(
            chunk.get("chunk_text", "")
            for chunk in sorted(chunks, key=lambda item: int(item.get("chunk_index", 0)))
        ) or None

        params_raw = payload.get("params")
        params: QeParameters | None = None
        if isinstance(params_raw, dict):
            params = QeParameters.model_validate(params_raw)

        created_at_raw = cast(str | None, payload.get("created_at"))
        created_at = (
            datetime.fromisoformat(created_at_raw)
            if created_at_raw
            else datetime.now(timezone.utc)
        )
        if created_at.tzinfo is None:
            created_at = created_at.replace(tzinfo=timezone.utc)

        cif = cast(str, payload["cif"])
        return StoredStructure(
            atoms=_atoms_from_cif(cif),
            source=cast(str, payload.get("source", "qe")),
            cif=cif,
            params=params,
            raw_input=raw_input,
            created_at=created_at,
            tenant_id=entry_tenant,
            workspace_id=entry_workspace,
        )


def _atoms_from_cif(cif_text: str) -> ASEAtoms:
    atoms = ase_read(StringIO(cif_text), format="cif")
    if not isinstance(atoms, ASEAtoms):
        raise ValueError("failed to parse CIF into ASE atoms")
    return atoms


def _build_store() -> InMemoryStructureStore | ConvexStructureStore:
    backend = os.getenv("STRUCTURE_STORE_BACKEND", "convex").strip().lower()
    if backend == "memory":
        return InMemoryStructureStore()
    if backend != "convex":
        raise ValueError("STRUCTURE_STORE_BACKEND must be 'convex' or 'memory'")

    convex_url = os.getenv("CONVEX_URL", "").strip()
    deploy_key = os.getenv("CONVEX_DEPLOY_KEY", "").strip()
    if not convex_url or not deploy_key:
        raise RuntimeError("CONVEX_URL and CONVEX_DEPLOY_KEY are required")
    chunk_size = int(os.getenv("STRUCTURE_RAW_CHUNK_SIZE", "120000"))
    timeout = int(os.getenv("STRUCTURE_STORE_TIMEOUT_SECONDS", "10"))
    return ConvexStructureStore(
        convex_url=convex_url,
        convex_deploy_key=deploy_key,
        chunk_size=chunk_size,
        timeout_seconds=timeout,
    )


_store: InMemoryStructureStore | ConvexStructureStore | None = None


def _get_store() -> InMemoryStructureStore | ConvexStructureStore:
    global _store
    if _store is None:
        _store = _build_store()
    return _store


def reload_structure_store() -> None:
    global _store
    _store = _build_store()


def create_structure_from_qe(
    content: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> tuple[str, Structure, str, QeParameters | None]:
    atoms, source = parse_qe_atoms(content)
    cif = atoms_to_cif(atoms)
    structure = structure_from_ase(atoms)
    params = extract_qe_params(content)
    structure_id = _get_store().create(
        atoms=atoms,
        source=source,
        cif=cif,
        params=params,
        raw_input=content,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
    return structure_id, structure, source, params


def get_structure_cif(
    structure_id: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> str:
    entry = _get_store().get(
        structure_id,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
    if not entry:
        raise KeyError(structure_id)
    return entry.cif


def get_structure_entry(
    structure_id: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> StoredStructure:
    entry = _get_store().get(
        structure_id,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
    if not entry:
        raise KeyError(structure_id)
    return entry


def get_structure(
    structure_id: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> Structure:
    entry = get_structure_entry(
        structure_id,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
    return structure_from_ase(entry.atoms)


def get_structure_params(
    structure_id: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> QeParameters | None:
    entry = get_structure_entry(
        structure_id,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
    return entry.params


def get_structure_raw_input(
    structure_id: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> str | None:
    entry = get_structure_entry(
        structure_id,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
    return entry.raw_input


def register_structure_atoms(
    atoms: ASEAtoms,
    source: str,
    *,
    tenant_id: str | None = None,
    workspace_id: str | None = None,
) -> str:
    cif = atoms_to_cif(atoms)
    return _get_store().create(
        atoms=atoms,
        source=source,
        cif=cif,
        params=None,
        raw_input=None,
        tenant_id=tenant_id,
        workspace_id=workspace_id,
    )
