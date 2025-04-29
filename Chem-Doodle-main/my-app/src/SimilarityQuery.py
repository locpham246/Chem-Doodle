#!/usr/bin/env python3
import sys
import json
import time
import numpy as np
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, udf
from pyspark.sql.types import DoubleType
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_fingerprint_blob(smiles_str):
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            arr = np.array([int(bit) for bit in fp.ToBitString()], dtype=np.uint8)
            return np.packbits(arr).tobytes()
    except Exception as e:
        print(f"Error processing {smiles_str}: {e}", file=sys.stderr)
    return None

# Precompute the query fingerprint once
query_smiles = None
query_blob = None
query_bits = None

def compute_similarity(blob):
    candidate_bits = np.unpackbits(np.frombuffer(blob, dtype=np.uint8))
    intersection = int(np.sum(np.logical_and(query_bits, candidate_bits)))
    union = int(np.sum(np.logical_or(query_bits, candidate_bits)))
    return float(intersection) / union if union else 0.0

def main():
    spark = SparkSession.builder \
        .appName("CassandraSimilarityQuery") \
        .config("spark.cassandra.connection.host", "96.125.118.113") \
        .config("spark.cassandra.connection.port", "9042") \
        .config("spark.jars.packages", "com.datastax.spark:spark-cassandra-connector_2.12:3.2.0") \
        .getOrCreate()
    
    spark.sparkContext.setLogLevel("ERROR")

    df = spark.read \
        .format("org.apache.spark.sql.cassandra") \
        .options(keyspace="drug_discovery", table="compounds") \
        .load() \
        .filter("cid > 0 AND cid < 10000000") \
        .limit(10000)

    # register UDF after query_bits is set
    compute_similarity_udf = udf(compute_similarity, DoubleType())

    df.cache()
    start_time = time.time()

    df_with_similarity = df.withColumn("similarity", compute_similarity_udf(col("fingerprint")))
    result_df = df_with_similarity \
        .filter(col("similarity") > 0.005) \
        .orderBy(col("similarity").desc()) \
        .select("cid") \
        .limit(10000)

    # collect the top CIDs
    cids = [row.cid for row in result_df.collect()]

    # output JSON array of CIDs
    print(json.dumps({ "cids": cids }))

    spark.stop()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python SimilarityQuery.py \"<query_smiles>\"", file=sys.stderr)
        sys.exit(1)

    query_smiles = sys.argv[1]
    query_blob = smiles_to_fingerprint_blob(query_smiles)
    if query_blob is None:
        print(f"Failed to compute fingerprint for query: {query_smiles}", file=sys.stderr)
        sys.exit(1)
    query_bits = np.unpackbits(np.frombuffer(query_blob, dtype=np.uint8))
    main()
