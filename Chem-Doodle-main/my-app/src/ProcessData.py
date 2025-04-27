#!/usr/bin/env python3

from pyspark.sql import SparkSession
from pyspark.sql.functions import split, col, udf
from pyspark.sql.types import BinaryType, StringType
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_fingerprint_blob(smiles_str):
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            bit_str = fp.ToBitString()
            arr = np.array([int(bit) for bit in bit_str], dtype=np.uint8)
            blob = np.packbits(arr)
            return blob.tobytes()
    except Exception as e:
        print(f"Error processing {smiles_str}: {e}")
    return None

fingerprint_udf = udf(smiles_to_fingerprint_blob, BinaryType())

def main():
    spark = SparkSession.builder \
        .appName("CID-SMILES-Ingestion") \
        .config("spark.cassandra.connection.host", "127.0.0.1") \
        .config("spark.jars.packages", "com.datastax.spark:spark-cassandra-connector_2.12:3.2.0") \
        .getOrCreate()

    df = spark.read.text("data/CID-SMILES")
    df = df.withColumn("parts", split(col("value"), r"\s+")) \
           .withColumn("cid", col("parts")[0].cast("int")) \
           .withColumn("smiles", col("parts")[1]) \
           .drop("value", "parts")

    df = df.withColumn("fingerprint", fingerprint_udf(col("smiles")))
    df = df.filter(col("fingerprint").isNotNull())

    df.write \
      .format("org.apache.spark.sql.cassandra") \
      .options(keyspace="drug_discovery", table="compounds") \
      .mode("append") \
      .save()

    spark.stop()

if __name__ == "__main__":
    main()
