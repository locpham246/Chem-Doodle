/* eslint-env node */
/* global process */
import express from 'express'
import { spawn } from 'child_process'

const app = express()
app.use((req, res, next) => {
  console.log(`--> Incoming Request: ${req.method} ${req.url}`); // Log method and URL
  // Log headers if needed for debugging (can be verbose)
  // console.log('--> Headers:', req.headers);
  next(); // IMPORTANT: Pass control to the next middleware/handler
});

const PORT = process.env.PORT || 5001

// built-in JSON body parser
app.use(express.json());

/**
 * Runs the SimilarityQuery.py script with the given SMILES
 * and returns a Promise that resolves to an array of CIDs.
 */
function runSimilarity(smiles) {
  return new Promise((resolve, reject) => {
    const py = spawn('python3', ['./src/SimilarityQuery.py', smiles]);
    let stdout = '';
    let stderr = '';

    py.stdout.on('data', data => {
      stdout += data.toString();
    });
    py.stderr.on('data', data => {
      stderr += data.toString();
    });

    py.on('close', code => {
      // --- START DEBUGGING ADDITION ---
      console.log(`\n--- Debug Info ---`);
      console.log(`Python script exited with code: ${code}`);
      console.log('--- Raw stdout received: ---');
      console.log('>>> START <<<');
      console.log(stdout); // Print the raw stdout content
      console.log('>>> END <<<');
      console.log('--- Raw stderr received: ---');
      console.log('>>> START <<<');
      console.log(stderr); // Print the raw stderr content too
      console.log('>>> END <<<');
      console.log(`----------------------\n`);
      // --- END DEBUGGING ADDITION ---

      if (code !== 0) {
        // Keep original error handling for non-zero exit code
        return reject(new Error(stderr.trim() || `Exited with code ${code}`));
      }

      // Now, try to parse the stdout we just logged
      try {
        // --- MODIFICATION START ---
        // 1. Trim whitespace and split stdout into an array of lines
        const outputLines = stdout.trim().split('\n');

        // 2. Get the very last line from the array
        const lastLine = outputLines[outputLines.length - 1];

        console.log("Attempting to parse last line of stdout:", lastLine); // Optional: for debugging

        if (!lastLine) {
          // Handle cases where the script might somehow produce no output
          throw new Error("Python script produced no output on stdout.");
        }

        // 3. Parse ONLY the last line as JSON
        const { cids } = JSON.parse(lastLine);
        // --- MODIFICATION END ---

        console.log("Successfully parsed JSON from last line.");
        resolve(cids); // Resolve with the extracted array of CIDs

      } catch (err) {
        // Error handling if parsing the last line fails
        console.error("Failed to parse the LAST LINE of stdout as JSON.");
        console.error("Raw stdout was:\n", stdout); // Log raw output for context
        reject(new Error(`Invalid JSON from script's last line: ${err.message}`));
      }
    });
  });
}

app.post('/api/similarity', async (req, res) => {
  const { smiles } = req.body;
  if (!smiles) {
    return res.status(400).json({ error: 'Missing SMILES in request body' });
  }

  console.log(`Received request for SMILES: ${smiles}`); // Add this log

  try {
    const cids = await runSimilarity(smiles);
    res.json({ cids });
  } catch (err) {
    console.error('Similarity script error:', err);
    res
      .status(500)
      .json({ error: 'Similarity script failed', detail: err.message });
  }

});

app.listen(PORT, () => {
  console.log(`Server listening on port ${PORT}`);
});
