#!/bin/bash
set -euo pipefail

# Variables
DB_NAME="aavc"
DB_USER="aavc_user"
DB_PASS="aavc_pass"

# Install R if missing
if ! command -v Rscript &>/dev/null; then
    sudo apt update && sudo apt install -y r-base-core
    sudo Rscript -e "install.packages('data.table', repos='https://cran.rstudio.com/')"
fi

# Install PostgreSQL if missing
if ! dpkg -l | grep -q "postgresql-16"; then
    sudo apt update
    sudo apt install -y gnupg2 wget curl
    echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" \
        | sudo tee /etc/apt/sources.list.d/pgdg.list
    curl -fsSL https://www.postgresql.org/media/keys/ACCC4CF8.asc \
        | sudo gpg --dearmor -o /etc/apt/trusted.gpg.d/postgresql.gpg
    sudo apt update
    sudo apt install -y postgresql-16 postgresql-contrib-16
    sudo systemctl enable --now postgresql
fi

# Create DB and user if missing
if ! sudo -u postgres psql -tAc "SELECT 1 FROM pg_database WHERE datname='$DB_NAME'" | grep -q 1; then
    sudo -u postgres psql -c "CREATE DATABASE $DB_NAME;"
fi

if ! sudo -u postgres psql -tAc "SELECT 1 FROM pg_roles WHERE rolname='$DB_USER'" | grep -q 1; then
    sudo -u postgres psql -c "CREATE USER $DB_USER WITH PASSWORD '$DB_PASS';"
    sudo -u postgres psql -c "ALTER DATABASE $DB_NAME OWNER TO $DB_USER;"
fi

# Get connection string
PORT=$(sudo -u postgres psql -tAc "SHOW port;")
DB_URL="postgresql://$DB_USER:$DB_PASS@localhost:$PORT/$DB_NAME"
echo "$DB_URL" > db_config.txt

# Install Python + deps
sudo apt install -y python3 python3-pip python3-venv
python3 -m venv venv
source venv/bin/activate
#install pandas
echo "Installing pandas for Python..."
pip3 install pandas
echo ""

#install sql utilites
echo "Installing sqlalchemy for Python..."
pip3 install sqlalchemy
echo ""

#install more sql utilites
echo "Installing psycopg2-binary for Python..."
pip3 install psycopg2-binary
echo ""

#install package for API queries
echo "Installing requests for Python..."
pip3 install requests
echo ""

# Run setup
python3 write_tables.py
psql "$DB_URL" -f indices.sql

echo "Setup complete!"
