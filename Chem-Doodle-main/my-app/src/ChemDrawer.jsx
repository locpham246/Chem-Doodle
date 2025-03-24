import React, { useEffect, useRef, useState } from 'react';
import data from './assets/data.json';  // data.json: array of objects with cid and smiles
import OCL from 'openchemlib/full.js';

// --- Helper Functions ---

// Canonicalize SMILES using OpenChemLib while preserving wildcard tokens.
// We temporarily replace "*" and "+" with placeholders so that OCL can process the rest.
// Then, we dearomatize the molecule so that benzene rings are represented with explicit double bonds (e.g. C1=CC=CC=C1)
const canonicalizeSmilesUsingOCL = (smiles) => {
  const starPlaceholder = '__STAR__';
  const plusPlaceholder = '__PLUS__';
  // Replace wildcards with placeholders.
  const placeholderSmiles = smiles.replace(/\*/g, starPlaceholder).replace(/\+/g, plusPlaceholder);
  try {
    const mol = OCL.Molecule.fromSmiles(placeholderSmiles);
    // Force the molecule to display in Kekulé form (explicit double bonds) by dearomatizing it.
    mol.dearomatize();
    let canonical = mol.getCanonizedSmiles();
    // Restore the wildcards.
    canonical = canonical.replace(new RegExp(starPlaceholder, 'g'), '*')
                         .replace(new RegExp(plusPlaceholder, 'g'), '+');
    return canonical;
  } catch (err) {
    console.error("Error canonicalizing SMILES:", err);
    return smiles;
  }
};

// Convert a (canonicalized) SMILES with wildcards to a RegExp pattern.
// "*" becomes ".*" (matches any sequence) and "+" becomes "." (matches exactly one character).
const convertSmilesToRegex = (smiles) => {
  const escapeChar = (char) => {
    if (char === '*' || char === '+') return char;
    return char.replace(/[-\/\\^$?.()|[\]{}]/g, '\\$&');
  };
  const escaped = Array.from(smiles).map(escapeChar).join('');
  const pattern = escaped.replace(/\*/g, '.*').replace(/\+/g, '.');
  return new RegExp(`^${pattern}$`, 'i');
};

const ChemDrawer = () => {
  const containerId = "jsme_container";
  const jsmeRef = useRef(null);
  const [smiles, setSmiles] = useState("");
  const [matches, setMatches] = useState(null); // null indicates no search performed yet

  useEffect(() => {
    window.jsmeOnLoad = () => {
      console.log("jsmeOnLoad called");
      if (window.JSApplet && window.JSApplet.JSME && typeof window.JSApplet.JSME === 'function') {
        try {
          const jsmeInstance = new window.JSApplet.JSME(
            containerId,
            "800px",
            "450px",
            { options: "query" }
          );
          console.log("JSME initialized.");
          jsmeRef.current = jsmeInstance;
        } catch (err) {
          console.error("Error initializing JSME:", err);
        }
      } else {
        console.error("Could not find JSME constructor. Check the JSME file.");
      }
    };

    if (!document.getElementById("jsme-script")) {
      const script = document.createElement("script");
      script.id = "jsme-script";
      script.src = "/jsme/jsme.nocache.js";
      script.async = true;
      script.onerror = () => {
        console.error("Failed to load JSME script from /jsme/jsme.nocache.js");
      };
      document.body.appendChild(script);
    }
  }, []);

  const handleParseSmiles = () => {
    if (jsmeRef.current) {
      const currentSmiles = jsmeRef.current.smiles();
      setSmiles(currentSmiles);
      console.log("SMILES:", currentSmiles);
      // Clear previous search results.
      setMatches(null);
    }
  };

  const handleFindSimilarity = () => {
    if (!smiles) {
      console.error("No SMILES available!");
      return;
    }
    // Canonicalize the drawn SMILES.
    const canonicalDrawnSmiles = canonicalizeSmilesUsingOCL(smiles);
    console.log("Canonical Drawn SMILES:", canonicalDrawnSmiles);

    // Determine if wildcards are present.
    const hasWildcards = /[\*\+]/.test(canonicalDrawnSmiles);
    let smilesMatchFunction;
    if (hasWildcards) {
      // Use regex matching.
      const regex = convertSmilesToRegex(canonicalDrawnSmiles);
      smilesMatchFunction = (candidateSmiles) => regex.test(candidateSmiles);
      console.log("Using regex pattern:", regex);
    } else {
      // Use exact match.
      smilesMatchFunction = (candidateSmiles) => candidateSmiles === canonicalDrawnSmiles;
    }

    // Filter candidate molecules based solely on SMILES matching.
    const found = data.filter(candidate => {
      let candidateCanonical;
      try {
        candidateCanonical = canonicalizeSmilesUsingOCL(candidate.smiles);
      } catch (err) {
        console.error("Error canonicalizing candidate SMILES:", err);
        candidateCanonical = candidate.smiles;
      }
      return smilesMatchFunction(candidateCanonical);
    });

    setMatches(found);
    console.log("Found candidates:", found);
  };

  return (
    <div style={{ textAlign: "center", marginTop: "2rem" }}>
      <h2>JSME Chemical Drawer</h2>
      <div id={containerId} style={{ margin: "0 auto" }}></div>
      <button onClick={handleParseSmiles} style={{ marginTop: "1rem" }}>
        Parse to SMILES
      </button>
      <button onClick={handleFindSimilarity} style={{ marginTop: "1rem", marginLeft: "1rem" }}>
        Find Similarity
      </button>
      {smiles && (
        <p>
          <strong>SMILES Output:</strong> {smiles}
        </p>
      )}
      {matches !== null && (
        <div style={{ marginTop: "1rem" }}>
          <h3>Matching Molecules</h3>
          <div style={{ display: "flex", flexWrap: "wrap", justifyContent: "center" }}>
            {matches.length > 0 ? (
              matches.map(candidate => (
                <div
                  key={candidate.cid}
                  style={{
                    border: "1px solid #ccc",
                    padding: "0.5rem 1rem",
                    margin: "0.5rem",
                    borderRadius: "4px",
                    cursor: "pointer"
                  }}
                  onClick={() => console.log("Clicked on CID:", candidate.cid)}
                >
                  CID: {candidate.cid}
                </div>
              ))
            ) : (
              <div
                style={{
                  border: "1px solid #ccc",
                  padding: "0.5rem 1rem",
                  margin: "0.5rem",
                  borderRadius: "4px"
                }}
              >
                No CID found that matches this structure.
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
};

export default ChemDrawer;