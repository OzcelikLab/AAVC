# ---------------------------
# Dockerfile for AAVC Project
# ---------------------------

# Use an official Debian-based image
FROM ubuntu:24.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# ---------------------------
# Install system dependencies
# ---------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        gnupg2 wget curl nano lsb-release \
        python3 python3-pip \
        r-base r-base-dev \
        postgresql-16 postgresql-contrib-16 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------
# Set environment variables
# ---------------------------
ENV POSTGRES_USER=aavc_user
ENV POSTGRES_PASSWORD=aavc_pass
ENV POSTGRES_DB=aavc
ENV POSTGRES_PORT=5432
ENV PATH="/usr/lib/postgresql/16/bin/:${PATH}"

# ---------------------------
# Install R and Python packages
# ---------------------------
RUN Rscript -e "install.packages('data.table', repos='https://cran.rstudio.com/')" && \
    pip3 install --no-cache-dir pandas sqlalchemy psycopg2-binary requests

# ---------------------------
# Copy project files
# ---------------------------
WORKDIR /opt/aavc
COPY . /opt/aavc

# ---------------------------
# Initialize PostgreSQL
# ---------------------------
USER postgres
RUN /etc/init.d/postgresql start && \
    psql -c "CREATE USER ${POSTGRES_USER} WITH PASSWORD '${POSTGRES_PASSWORD}';" || true && \
    psql -c "CREATE DATABASE ${POSTGRES_DB} OWNER ${POSTGRES_USER};" || true

# ---------------------------
# Switch back to root for container commands
# ---------------------------
USER root

# ---------------------------
# Expose PostgreSQL port
# ---------------------------
EXPOSE 5432

# ---------------------------
# Default command: initialize DB & show help
# ---------------------------
CMD service postgresql start && \
    python3 write_tables.py && \
    psql -U ${POSTGRES_USER} -d ${POSTGRES_DB} -f indices.sql && \
    python3 aavc6.py --help
