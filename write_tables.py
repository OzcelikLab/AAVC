import os
import pandas as pd
from sqlalchemy import create_engine, Index, MetaData, Table, text
from sqlalchemy.exc import SQLAlchemyError

with open('db_config.txt', 'r') as file:
    db_url = file.read().strip()

engine = create_engine(db_url)
connection = engine.connect()


def read_csv_to_sql(file_path, engine):
    table_name = os.path.splitext(os.path.basename(file_path))[0]

    try:
        if table_name == "homopolymers1":
            df = pd.read_csv(file_path, sep="\t", low_memory=False)
            df[df.columns[0]] = df[df.columns[0]].astype(str)
            print(f"Writing '{table_name}' to PostgreSQL...")
            df.to_sql("homopolymers", engine, if_exists='replace', index=False)
            print(f"Table '{table_name}' successfully written to PostgreSQL.")
        elif "homopolymers" in table_name:
            df = pd.read_csv(file_path, sep="\t", low_memory=False)
            print(f"Writing '{table_name}' to PostgreSQL...")
            df.to_sql("homopolymers", engine, if_exists='append', index=False)
            print(f"Table '{table_name}' successfully written to PostgreSQL.")
        elif table_name == "repeats1":
            df = pd.read_csv(file_path, sep="\t", low_memory=False)
            df[df.columns[0]] = df[df.columns[0]].astype(str)
            print(f"Writing '{table_name}' to PostgreSQL...")
            df.to_sql("repeats", engine, if_exists='replace', index=False)
            print(f"Table '{table_name}' successfully written to PostgreSQL.")
        elif "repeats" in table_name:
            df = pd.read_csv(file_path, sep="\t", low_memory=False)
            print(f"Writing '{table_name}' to PostgreSQL...")
            df.to_sql("repeats", engine, if_exists='append', index=False)
            print(f"Table '{table_name}' successfully written to PostgreSQL.")
        else:
            df = pd.read_csv(file_path, sep="\t", low_memory=False)
            print(f"Writing '{table_name}' to PostgreSQL...")
            df.to_sql(table_name, engine, if_exists='replace', index=False)
            print(f"Table '{table_name}' successfully written to PostgreSQL.")
    except SQLAlchemyError as e:
        print(f"Error occurred while writing '{table_name}' to PostgreSQL: {e}")

for file_name in os.listdir("db/"):
    if file_name.endswith('.txt'):
        file_path = os.path.join("db/", file_name)
        read_csv_to_sql(file_path, engine)
        print("\n")

engine.dispose()
