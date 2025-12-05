# ---- base python ----
FROM python:3.13.7-alpine3.22 as python-base

ENV PYTHONUNBUFFERED=1
ENV HOME=/code


# ---- build stage with uv ----
FROM python-base as uv-base

# Install uv in builder stage only
RUN pip install --no-cache-dir uv

WORKDIR /
COPY pyproject.toml uv.lock ./

# Create a virtualenv inside /.venv and install only main deps
RUN uv venv /.venv && \
    uv pip install --no-cache --python=/.venv/bin/python -r <(uv pip compile --only-main --no-emit-project -q uv.lock)


# ---- final stage ----
FROM python-base as final

# Copy over the ready-to-use venv
COPY --from=uv-base /.venv /.venv

# Point runtime to the venv
ENV PYTHONPATH="${PYTHONPATH}:/.venv/lib/python3.13/site-packages/"
ENV PATH=/.venv/bin:$PATH

WORKDIR ${HOME}
COPY src/ ./
