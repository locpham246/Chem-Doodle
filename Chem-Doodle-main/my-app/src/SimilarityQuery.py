# src/spark_query.py

from pyspark.sql import SparkSession
from pyspark.sql.functions import col
from pyspark.sql.types import DoubleType
import time
import numpy as np
import pandas as pd
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

query_smiles = "c1ccccc1"
query_blob = smiles_to_fingerprint_blob(query_smiles)
if query_blob is None:
    raise ValueError("Failed to compute fingerprint for query compound")
else:
    query_bits = np.unpackbits(np.frombuffer(query_blob, dtype=np.uint8))
    print(f"Query SMILES: {query_smiles}")
    print(f"Query fingerprint length (in bits): {len(query_bits)}")

from pyspark.sql.functions import pandas_udf

@pandas_udf(DoubleType())
def compute_similarity(candidate_blob_series: pd.Series) -> pd.Series:
    query_bits = np.unpackbits(np.frombuffer(query_blob, dtype=np.uint8))
    similarities = []
    for i, blob in enumerate(candidate_blob_series):
        candidate_bits = np.unpackbits(np.frombuffer(blob, dtype=np.uint8))
        intersection = np.sum(np.logical_and(query_bits, candidate_bits))
        union = np.sum(np.logical_or(query_bits, candidate_bits))
        sim = float(intersection) / union if union != 0 else 0.0
        if i < 3:
            print(f"Candidate {i}: intersection={intersection}, union={union}, similarity={sim}")
        similarities.append(sim)
    return pd.Series(similarities)

def main():
    spark = SparkSession.builder \
        .appName("CassandraSimilarityQuery") \
        .config("spark.cassandra.connection.host", "127.0.0.1") \
        .config("spark.jars.packages", "com.datastax.spark:spark-cassandra-connector_2.12:3.2.0") \
        .getOrCreate()
    
    print("Reading data from Cassandra...")
    
    df = spark.read \
    .format("org.apache.spark.sql.cassandra") \
    .options(keyspace="drug_discovery", table="compounds") \
    .load() \
    .filter("cid > 0 AND cid < 10000000")

    df = df.limit(10000000)
    
    record_count = df.count()
    print(f"Number of records read from Cassandra (limited): {record_count}")
    
    df.cache()
    
    start_time = time.time()
    
    print("Computing similarities...")
    df_with_similarity = df.withColumn("similarity", compute_similarity(col("fingerprint")))
    
    result_df = df_with_similarity.filter(col("similarity") > 0.2).orderBy(col("similarity").desc())
    
    result_count = result_df.count()
    print(f"Number of compounds with similarity > 0.5: {result_count}")
    
    result_df.show(10, truncate=False)
    
    end_time = time.time()
    print("Similarity query took {:.2f} seconds".format(end_time - start_time))
    
    spark.stop()

if __name__ == "__main__":
    main()

