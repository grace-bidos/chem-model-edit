from __future__ import annotations


def run_http_worker() -> None:
    raise RuntimeError(
        "HTTP worker runtime is retired. Use /api/runtime/* with user-managed AiiDA/Slurm."
    )


if __name__ == "__main__":
    run_http_worker()
