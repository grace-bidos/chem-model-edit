FROM python:3.13-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends nodejs npm \
    && rm -rf /var/lib/apt/lists/*

RUN npm install --no-audit --no-fund --prefix /opt/molstar molstar@5.5.0
ENV MOLSTAR_CIF2BCIF_PATH=/opt/molstar/node_modules/molstar/lib/commonjs/cli/cif2bcif/index.js \
    NODE_BIN=node

RUN pip install --no-cache-dir uv

COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-dev

COPY . .

ENV PORT=8080
EXPOSE 8080

CMD ["uv", "run", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]
