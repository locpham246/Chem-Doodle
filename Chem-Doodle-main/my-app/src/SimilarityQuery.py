#!/usr/bin/env python3
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, udf
from pyspark.sql.types import DoubleType
import time
import numpy as np
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
        print(f"Error processing {smiles_str}: {e}")
    return None

query_smiles = "c1ccccc1"
query_blob = smiles_to_fingerprint_blob(query_smiles)
if query_blob is None:
    raise ValueError("Failed to compute fingerprint for query compound")
query_bits = np.unpackbits(np.frombuffer(query_blob, dtype=np.uint8))
print(f"Query SMILES: {query_smiles}")
print(f"Query fingerprint length (in bits): {len(query_bits)}")

def compute_similarity(blob):
    candidate_bits = np.unpackbits(np.frombuffer(blob, dtype=np.uint8))
    intersection = int(np.sum(np.logical_and(query_bits, candidate_bits)))
    union = int(np.sum(np.logical_or(query_bits, candidate_bits)))
    return float(intersection) / union if union else 0.0

compute_similarity_udf = udf(compute_similarity, DoubleType())

def main():
    spark = SparkSession.builder \
        .appName("CassandraSimilarityQuery") \
        .config("spark.cassandra.connection.host", "96.125.118.113") \
        .config("spark.cassandra.connection.port", "9042") \
        .config("spark.jars.packages", "com.datastax.spark:spark-cassandra-connector_2.12:3.2.0") \
        .getOrCreate()

    df = spark.read \
        .format("org.apache.spark.sql.cassandra") \
        .options(keyspace="drug_discovery", table="compounds") \
        .load() \
        .filter("cid > 0 AND cid < 10000000") \
        .limit(10000000)

    record_count = df.count()
    print(f"Number of records read: {record_count}")

    df.cache()
    start_time = time.time()

    df_with_similarity = df.withColumn("similarity", compute_similarity_udf(col("fingerprint")))
    result_df = df_with_similarity.filter(col("similarity") > 0.2).orderBy(col("similarity").desc())

    result_count = result_df.count()
    print(f"Number of compounds with similarity > 0.2: {result_count}")
    result_df.show(10, truncate=False)

    elapsed = time.time() - start_time
    print(f"Similarity query took {elapsed:.2f} seconds")

    spark.stop()

if __name__ == "__main__":
    main()
