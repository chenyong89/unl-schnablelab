FROM python:3.10.10

COPY src /usr/src
COPY pyproject.toml /usr/
ENV PYTHONPATH "${PYTHONPATH}:/usr/src"
WORKDIR /usr
RUN pip install .