#!/usr/bin/env python3
import sys
import json
import time
# RDKit produces a lot of stderr warnings we often want to ignore in production
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem
import numpy as np # Although numpy isn't directly used here now, Spark might need it implicitly
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, udf
from pyspark.sql.types import DoubleType

# --- Silence RDKit messages ---
# You might want this in production to keep stderr clean
# rdBase.DisableLog('rdApp.*')
# print("RDKit logging disabled.", file=sys.stderr) # Optional: confirm disabling

# --- RDKit Fingerprint Function ---
# Calculate Morgan Fingerprint directly from SMILES
def calculate_fp(smiles_str):
    """Calculates RDKit Morgan Fingerprint from SMILES. Returns None on error."""
    if not smiles_str: # Handle empty or null SMILES
        return None
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol:
            # Using Morgan fingerprint as in the previous script (radius=2, 2048 bits)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            return fp
        else:
            # Optionally log invalid SMILES that don't produce a molecule object
            # print(f"Warning: Could not create molecule from SMILES: {smiles_str}", file=sys.stderr)
            return None
    except Exception as e:
        # Log other RDKit errors
        print(f"RDKit Error processing SMILES '{smiles_str}': {e}", file=sys.stderr)
    return None

# --- Precompute Query Fingerprint ---
# These will be set when the script starts
query_smiles = None
query_fp = None

# --- UDF for Similarity Calculation ---
def compute_structure_similarity(candidate_smiles):
    """Calculates Tanimoto similarity between precomputed query_fp and candidate_smiles."""
    # Check if global query_fp is valid and candidate SMILES is provided
    if query_fp is None or candidate_smiles is None:
        return 0.0

    # Calculate fingerprint for the candidate molecule from the database
    candidate_fp = calculate_fp(candidate_smiles)
    if candidate_fp is None:
        # If candidate SMILES is invalid or fails fingerprinting, similarity is 0
        return 0.0

    try:
        # Calculate Tanimoto Similarity (returns float between 0.0 and 1.0)
        similarity = DataStructs.TanimotoSimilarity(query_fp, candidate_fp)
        return float(similarity)
    except Exception as e:
        # Log unexpected errors during similarity calculation
        print(f"Similarity Calculation Error: {e}", file=sys.stderr)
        return 0.0

# --- Main Function ---
def main():
    """Main logic to connect to Spark, read data, calculate similarity, and output results."""
    spark = SparkSession.builder \
        .appName("CassandraStructureSimilarity") \
        .config("spark.cassandra.connection.host", "96.125.118.113") \
        .config("spark.cassandra.connection.port", "9042") \
        .config("spark.jars.packages", "com.datastax.spark:spark-cassandra-connector_2.12:3.2.0") \
        .getOrCreate()

    # Set Spark log level to avoid excessive messages polluting stdout/stderr
    spark.sparkContext.setLogLevel("ERROR")
    print("Spark log level set to ERROR.", file=sys.stderr) # Optional confirmation

    # --- Read SMILES from Cassandra ---
    # Make sure your table 'compounds' has 'cid' and 'smiles' columns
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

        # Let's see how many rows we have after loading and filtering
        count = df.count()
        print(f"Loaded {count} rows with non-null SMILES.", file=sys.stderr)
        if count == 0:
             print("No data loaded, exiting.", file=sys.stderr)
             print(json.dumps({"results": []})) # Output empty results
             spark.stop()
             sys.exit(0)


    except Exception as e:
         print(f"Error reading from Cassandra: {e}", file=sys.stderr)
         print(json.dumps({"error": f"Cassandra read failed: {e}"})) # Output error JSON
         spark.stop()
         sys.exit(1)

    # --- Register UDF ---
    # Make sure query_fp was successfully calculated before registering the UDF
    if query_fp is None:
         print(f"Error: Query fingerprint is invalid. Cannot proceed.", file=sys.stderr)
         print(json.dumps({"error": "Invalid query SMILES."})) # Output error JSON
         spark.stop()
         sys.exit(1)

    # Define the UDF that Spark will use to call our Python function
    compute_similarity_udf = udf(compute_structure_similarity, DoubleType())
    print("Spark UDF registered.", file=sys.stderr)

    # --- Calculate Similarity ---
    # Cache DataFrame for potentially better performance if reused
    # df.cache() # Caching might be less useful if only used once like here
    print("Calculating similarities...", file=sys.stderr)
    start_time = time.time()

    try:
        df_with_similarity = df.withColumn(
            "similarity",
            compute_similarity_udf(col("smiles")) # Apply UDF to the 'smiles' column
        )

        # --- Filter and Sort Results ---
        # Define your desired similarity threshold and result limit
        similarity_threshold = 0.2 # Example: Find compounds >= 70% similar
        result_limit = 12

        print(f"Filtering results with similarity >= {similarity_threshold}", file=sys.stderr)
        result_df = df_with_similarity \
            .filter(col("similarity") >= similarity_threshold) \
            .orderBy(col("similarity").desc()) \
            .select("cid", "similarity") \
            .limit(result_limit)

        # --- Collect Results ---
        print("Collecting results...", file=sys.stderr)
        # Collect results into a list of dictionaries: [{"cid": ..., "similarity": ...}, ...]
        results = [
            {"cid": row.cid, "similarity": round(row.similarity, 4)}
            for row in result_df.collect()
        ]
        result_count = len(results)
        print(f"Found {result_count} similar compounds.", file=sys.stderr)
        print(f"Calculation took {time.time() - start_time:.2f} seconds.", file=sys.stderr)

        # --- Output JSON ---
        # Output format includes both CID and similarity score
        print(json.dumps({"results": results}))

    except Exception as e:
        print(f"Error during Spark similarity calculation: {e}", file=sys.stderr)
        # Output error JSON if Spark processing fails
        print(json.dumps({"error": f"Spark processing failed: {e}"}))

    # Stop the Spark session
    spark.stop()
    print("Spark session stopped.", file=sys.stderr)


# --- Script Entry Point ---
if __name__ == "__main__":
    # Check if SMILES argument is provided
    if len(sys.argv) < 2:
        print("Usage: python StructureSimilarity.py \"<query_smiles>\"", file=sys.stderr)
        sys.exit(1)

    query_smiles = sys.argv[1]
    print(f"Received query SMILES: {query_smiles}", file=sys.stderr)

    # Calculate fingerprint for the query SMILES right at the start
    query_fp = calculate_fp(query_smiles)

    if query_fp is None:
        print(f"Failed to compute fingerprint for query SMILES: {query_smiles}. Please provide a valid SMILES string.", file=sys.stderr)
        # Output error JSON if query is invalid
        print(json.dumps({"error": "Invalid query SMILES provided."}))
        sys.exit(1) # Exit if query is invalid

    print("Query fingerprint calculated successfully.", file=sys.stderr)
    # Call the main function to execute the Spark job
    main()