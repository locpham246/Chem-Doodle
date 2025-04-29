import sys
import json
import time
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem
import numpy as np 
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, udf
from pyspark.sql.types import DoubleType

def calculate_fp(smiles_str):
    """Calculates RDKit Morgan Fingerprint from SMILES. Returns None on error."""
    if not smiles_str:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            return fp
        else:
            return None
    except Exception as e:
        print(f"RDKit Error processing SMILES '{smiles_str}': {e}", file=sys.stderr)
    return None

query_smiles = None
query_fp = None

def compute_structure_similarity(candidate_smiles):
    if query_fp is None or candidate_smiles is None:
        return 0.0

    candidate_fp = calculate_fp(candidate_smiles)
    if candidate_fp is None:
        return 0.0

    try:
        similarity = DataStructs.TanimotoSimilarity(query_fp, candidate_fp)
        return float(similarity)
    except Exception as e:
        print(f"Similarity Calculation Error: {e}", file=sys.stderr)
        return 0.0

def main():
    """Main logic to connect to Spark, read data, calculate similarity, and output results."""
    spark = SparkSession.builder \
        .appName("CassandraStructureSimilarity") \
        .config("spark.cassandra.connection.host", "96.125.118.113") \
        .config("spark.cassandra.connection.port", "9042") \
        .config("spark.jars.packages", "com.datastax.spark:spark-cassandra-connector_2.12:3.2.0") \
        .getOrCreate()

    spark.sparkContext.setLogLevel("ERROR")
    print("Spark log level set to ERROR.", file=sys.stderr) 

    print("Reading data from Cassandra...", file=sys.stderr)
    try:
        df = spark.read \
            .format("org.apache.spark.sql.cassandra") \
            .options(keyspace="drug_discovery", table="compounds") \
            .load() \
            .select("cid", "smiles") \
            .filter(col("smiles").isNotNull() & (col("smiles") != "")) \
            .filter("cid > 0 AND cid < 10000000") \
            .limit(2000)

        count = df.count()
        print(f"Loaded {count} rows with non-null SMILES.", file=sys.stderr)
        if count == 0:
             print("No data loaded, exiting.", file=sys.stderr)
             print(json.dumps({"results": []})) 
             spark.stop()
             sys.exit(0)


    except Exception as e:
         print(f"Error reading from Cassandra: {e}", file=sys.stderr)
         print(json.dumps({"error": f"Cassandra read failed: {e}"})) 
         spark.stop()
         sys.exit(1)

    if query_fp is None:
         print(f"Error: Query fingerprint is invalid. Cannot proceed.", file=sys.stderr)
         print(json.dumps({"error": "Invalid query SMILES."})) 
         spark.stop()
         sys.exit(1)

    compute_similarity_udf = udf(compute_structure_similarity, DoubleType())
    print("Spark UDF registered.", file=sys.stderr)

    print("Calculating similarities...", file=sys.stderr)
    start_time = time.time()

    try:
        df_with_similarity = df.withColumn(
            "similarity",
            compute_similarity_udf(col("smiles"))
        )

        similarity_threshold = 0.2 
        result_limit = 12

        print(f"Filtering results with similarity >= {similarity_threshold}", file=sys.stderr)
        result_df = df_with_similarity \
            .filter(col("similarity") >= similarity_threshold) \
            .orderBy(col("similarity").desc()) \
            .select("cid", "similarity") \
            .limit(result_limit)

        print("Collecting results...", file=sys.stderr)
        results = [
            {"cid": row.cid, "similarity": round(row.similarity, 4)}
            for row in result_df.collect()
        ]
        result_count = len(results)
        print(f"Found {result_count} similar compounds.", file=sys.stderr)
        print(f"Calculation took {time.time() - start_time:.2f} seconds.", file=sys.stderr)

        print(json.dumps({"results": results}))

    except Exception as e:
        print(f"Error during Spark similarity calculation: {e}", file=sys.stderr)
        print(json.dumps({"error": f"Spark processing failed: {e}"}))

    spark.stop()
    print("Spark session stopped.", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python StructureSimilarity.py \"<query_smiles>\"", file=sys.stderr)
        sys.exit(1)

    query_smiles = sys.argv[1]
    print(f"Received query SMILES: {query_smiles}", file=sys.stderr)

    query_fp = calculate_fp(query_smiles)

    if query_fp is None:
        print(f"Failed to compute fingerprint for query SMILES: {query_smiles}. Please provide a valid SMILES string.", file=sys.stderr)
        print(json.dumps({"error": "Invalid query SMILES provided."}))
        sys.exit(1) 

    print("Query fingerprint calculated successfully.", file=sys.stderr)
    main()