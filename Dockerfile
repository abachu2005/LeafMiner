FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ssh \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml webapp/requirements.txt ./webapp/
RUN pip install --upgrade pip \
    && pip install -r webapp/requirements.txt numpy pandas matplotlib scipy requests

COPY . .

# webapp/data is where job artifacts go — mount as a volume for persistence
VOLUME ["/app/webapp/data"]
EXPOSE 8000

CMD ["python", "-m", "uvicorn", "webapp.backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
