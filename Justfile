set shell := ["bash", "-cu"]
uv_cache_dir := ".uv-cache"

@default:
  just --list

web:
  pnpm -C apps/web dev --port 3001

api:
  cd apps/api && UV_CACHE_DIR={{uv_cache_dir}} uv sync
  cd apps/api && UV_CACHE_DIR={{uv_cache_dir}} uv run uvicorn main:app --reload --port 8000

dev:
  trap 'kill 0' INT TERM
  cd apps/api && UV_CACHE_DIR={{uv_cache_dir}} uv sync
  pnpm -C apps/web dev --port 3001 &
  cd apps/api && UV_CACHE_DIR={{uv_cache_dir}} uv run uvicorn main:app --reload --port 8000 &
  wait
