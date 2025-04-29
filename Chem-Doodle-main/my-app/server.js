/* eslint-env node */
/* global process */
import express from 'express';
import { spawn } from 'child_process';

const app = express();

// --- Middleware ---
// Optional: Basic request logger
app.use((req, res, next) => {
  console.log(`--> Incoming Request: ${req.method} ${req.url}`);
  next();
});

const PORT = process.env.PORT || 5001; // Using port 5001

// Built-in JSON body parser - IMPORTANT: place before routes that need req.body
app.use(express.json());


// --- Function to run Fingerprint Similarity Script ---
function runSimilarity(smiles) {
  return new Promise((resolve, reject) => {
    const py = spawn('python3', ['./src/SimilarityQuery.py', smiles]); // Original script
    let stdout = '';
    let stderr = '';

    py.stdout.on('data', data => { stdout += data.toString(); });
    py.stderr.on('data', data => { stderr += data.toString(); });

    py.on('close', code => {
       // Optional Debug logging
      console.log(`\n--- [Fingerprint Sim] Debug Info ---`);
      console.log(`[Fingerprint Sim] Python script exited with code: ${code}`);
      // console.log('[Fingerprint Sim] --- Raw stdout received: ---\n', stdout); // Uncomment for deep debug
      // console.log('[Fingerprint Sim] --- Raw stderr received: ---\n', stderr); // Uncomment for deep debug
      // console.log(`[Fingerprint Sim] --------------------------\n`);


      if (code !== 0) {
         // Try to parse stderr for JSON error first
         try {
             const errorJson = JSON.parse(stderr.trim());
             if (errorJson.error) return reject(new Error(errorJson.error));
         } catch(e) { /* Ignore parsing error */ }
        return reject(new Error(stderr.trim() || `Fingerprint script exited with code ${code}`));
      }

      try {
        const lines = stdout.trim().split('\n');
        const lastLine = lines[lines.length - 1];
        if (!lastLine) throw new Error("Fingerprint script stdout was empty.");

        const parsedJson = JSON.parse(lastLine); // Expects {"cids": [...]}
        if (parsedJson.cids === undefined) {
             if (parsedJson.error) throw new Error(`Python script returned error: ${parsedJson.error}`);
             throw new Error("Parsed JSON from script is missing 'cids' key.");
        }
        const cids = parsedJson.cids;
        console.log("[Fingerprint Sim] Successfully parsed cids from last line of stdout.");
        resolve(cids); // Resolve with the CIDs array
      } catch (err) {
        console.error("[Fingerprint Sim] Failed to parse LAST LINE of stdout as JSON or missing 'cids'.");
        console.error("[Fingerprint Sim] Raw stdout was:\n", stdout);
        reject(new Error(`Invalid JSON or format from Fingerprint script's last line: ${err.message}`));
      }
    });
  });
}


// --- NEW: Function to run Structure (SMILES) Similarity Script ---
function runSmilesSimilarity(smiles) {
  return new Promise((resolve, reject) => {
    // Spawn the NEW Python script
    const py = spawn('python3', ['./src/SimilaritySMILES.py', smiles]);
    let stdout = '';
    let stderr = '';

    py.stdout.on('data', data => { stdout += data.toString(); });
    py.stderr.on('data', data => { stderr += data.toString(); });

    py.on('close', code => {
      // Optional Debug logging
      console.log(`\n--- [SMILES Sim] Debug Info ---`);
      console.log(`[SMILES Sim] Python script exited with code: ${code}`);
      // console.log('[SMILES Sim] --- Raw stdout received: ---\n', stdout); // Uncomment for deep debug
      // console.log('[SMILES Sim] --- Raw stderr received: ---\n', stderr); // Uncomment for deep debug
      // console.log(`[SMILES Sim] --------------------------\n`);


      if (code !== 0) {
        // Handle script errors (non-zero exit code)
         // Try to parse stderr for JSON error first
         try {
             const errorJson = JSON.parse(stderr.trim());
             if (errorJson.error) return reject(new Error(errorJson.error));
         } catch(e) { /* Ignore parsing error */ }
        return reject(new Error(stderr.trim() || `SMILES Sim script exited with code ${code}`));
      }

      // Handle successful exit (code 0)
      try {
        const lines = stdout.trim().split('\n');
        const lastLine = lines[lines.length - 1];
        if (!lastLine) throw new Error("SMILES Sim script stdout was empty.");

        // Parse the last line expecting {"results": [...]}
        const parsedJson = JSON.parse(lastLine);

        // Check if the expected 'results' key exists
        if (parsedJson.results === undefined) {
             if (parsedJson.error) throw new Error(`Python script returned error: ${parsedJson.error}`);
             throw new Error("Parsed JSON from script is missing 'results' key.");
        }
        const results = parsedJson.results; // Extract the results array
        console.log("[SMILES Sim] Successfully parsed results from last line of stdout.");
        resolve(results); // Resolve with the results array [{cid:..., similarity:...}]

      } catch (err) {
        console.error("[SMILES Sim] Failed to parse LAST LINE of stdout as JSON or missing 'results'.");
        console.error("[SMILES Sim] Raw stdout was:\n", stdout);
        reject(new Error(`Invalid JSON or format from SMILES Sim script's last line: ${err.message}`));
      }
    });
  });
}


// --- API Routes ---

// Existing route for Fingerprint Similarity
app.post('/api/similarity', async (req, res) => {
  console.log("!!! Handling /api/similarity request (Fingerprint) !!!");
  const { smiles } = req.body;
  if (!smiles) {
    return res.status(400).json({ error: 'Missing SMILES in request body' });
  }
  console.log(`Received request for SMILES (Fingerprint route): ${smiles}`);

  try {
    const cids = await runSimilarity(smiles); // Calls the original function
    res.json({ cids }); // Returns { cids: [...] }
  } catch (err) {
    console.error('Fingerprint Similarity script error:', err);
    res
      .status(500)
      .json({ error: 'Fingerprint Similarity script failed', detail: err.message });
  }
});

// NEW route for Structure (SMILES) Similarity
app.post('/api/similarity_smiles', async (req, res) => {
  console.log("!!! Handling /api/similarity_smiles request !!!");
  const { smiles } = req.body;
  if (!smiles) {
    console.error("Missing SMILES for /api/similarity_smiles");
    return res.status(400).json({ error: 'Missing SMILES in request body' });
  }
  console.log(`Received request for SMILES (SMILES route): ${smiles}`);

  try {
    const results = await runSmilesSimilarity(smiles); // Calls the NEW function
    res.json({ results }); // Returns { results: [{cid:..., similarity:...}, ...] }

  } catch (err) {
    console.error('SMILES Similarity script error:', err);
    // Send back the error message captured from the promise rejection
    res
      .status(500)
      .json({ error: 'SMILES Similarity script failed', detail: err.message });
  }
});


// --- Start Server ---
app.listen(PORT, () => {
  console.log(`Server listening on port ${PORT}`);
});