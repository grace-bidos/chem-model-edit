FROM python:3.13-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir uv

COPY pyproject.toml uv.lock ./

RUN uv venv /opt/venv \
    && . /opt/venv/bin/activate \
    && uv sync --frozen --no-dev --no-install-project --active

ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"

COPY . /app

EXPOSE 8080

CMD ["uvicorn", "main:app", "--host=0.0.0.0", "--port=8080"]
