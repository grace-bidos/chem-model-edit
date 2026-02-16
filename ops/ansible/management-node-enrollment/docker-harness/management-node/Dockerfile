FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        git \
        openssh-client \
        sudo \
        sshpass \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir ansible

WORKDIR /workspace

CMD ["sleep", "infinity"]
