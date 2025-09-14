#!/bin/bash

#install R
if command -v Rscript &> /dev/null; then
    echo "Rscript is already installed."
    echo ""
else
    echo "Rscript is not installed. Installing r-base-core..."
    sudo apt install -y r-base-core
    echo ""
fi

#install data.table
if Rscript -e "require(data.table)" &>/dev/null; then
    echo "data.table for R is already installed"
else
    echo "Installing data.table package for R..."
    Rscript -e "install.packages('data.table', repos='https://cran.rstudio.com/')"
fi
echo ""

#install postgres
if dpkg -l | grep -q "postgresql-16"; then
    echo "PostgreSQL 16 is already installed."
    echo ""
else
    #install required packages
    sudo apt update
    sudo apt install gnupg2 wget nano curl

    #add repo
    sudo sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list'

    #add repo key
    curl -fsSL https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo gpg --dearmor -o /etc/apt/trusted.gpg.d/postgresql.gpg

    #update package list
    sudo apt update

    #install postgres
    sudo apt install postgresql-16 postgresql-contrib-16

    #verify installation, start if down & check the version
    sudo systemctl start postgresql
    sudo systemctl enable postgresql
    sudo service postgresql status
    psql --version
fi

#create database, add user, grant privileges
sudo -u postgres psql -c "CREATE DATABASE aavc;"
sudo -u postgres psql -c "CREATE USER mark_antony WITH PASSWORD 'caesar';"
sudo -u postgres psql -c "ALTER DATABASE aavc OWNER TO mark_antony;"

#140183

#create db config file
PORT=$(sudo grep -E "^port =" /etc/postgresql/16/main/postgresql.conf | awk '{print $3}')
DB_URL="postgresql://mark_antony:caesar@localhost:$PORT/aavc"
echo $DB_URL > db_config.txt
echo "Configuration file created: db_config.txt"
echo "Database and user setup complete."
echo ""

####sudo apt install python3-pandas
### pip3 install pandas --break-system-packages

#install python
echo "Installing python3..."
sudo apt install -y "python3"
echo ""

#install pip
echo "Installing python3-pip..."
sudo apt install -y "python3-pip"
echo ""

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

#write tables into sql
python3 write_tables.py
echo ""

#create indices
echo "Creating indices for PostgreSQL tables..."

#read config file
CONNECTION_STRING=$(cat db_config.txt)

#parse config file
USER=$(echo $CONNECTION_STRING | sed -r 's|^[^:]*://([^:]*):.*$|\1|')
PASS=$(echo $CONNECTION_STRING | sed -r 's|^[^:]*://[^:]*:([^@]*)@.*$|\1|')
HOST=$(echo $CONNECTION_STRING | sed -r 's|^[^:]*://[^@]*@([^:]*):.*$|\1|')
PORT=$(echo $CONNECTION_STRING | sed -r 's|^[^:]*://[^@]*@[^:]*:([^/]*)/.*$|\1|')
DB=$(echo $CONNECTION_STRING | sed -r 's|^[^:]*://[^/]*/*/(.*)$|\1|')

#export password
export PGPASSWORD=$PASS

#read indexing commands
SQL_COMMANDS=$(cat indices.sql)

#execute indexing
psql -h $HOST -p $PORT -U $USER -d $DB -c "$SQL_COMMANDS"

#unset password
unset PGPASSWORD
echo ""

echo "Installation is done."
echo ""

python3 aavc.py --help
