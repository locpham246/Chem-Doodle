import React, { useEffect, useRef, useState } from "react";


export default function ChemDrawer() {
  const containerId = "jsme_container";
  const jsmeRef = useRef(null);
  const [smiles, setSmiles] = useState("");
  // 'matches' state will now hold an array of {cid: num} OR {cid: num, similarity: num}
  const [matches, setMatches] = useState(null);
  const [isLoading, setIsLoading] = useState(false); // Added loading state
  const [error, setError] = useState(null); // Added error state
  const [, setTheme] = useState("light"); // Keep theme logic if needed

  useEffect(() => {
    // Theme logic remains the same
    const saved = localStorage.getItem("theme") || "light";
    setTheme(saved);
    document.body.setAttribute("data-theme", saved);

    // JSME initialization logic remains the same
    window.jsmeOnLoad = () => {
      if (!jsmeRef.current && window.JSApplet?.JSME) {
        try {
          jsmeRef.current = new window.JSApplet.JSME(
            containerId, "750px", "450px", { options: "query" }
          );
        } catch (err) {
          console.error("JSME init error:", err);
          setError("Failed to load chemical drawer."); // Set error state
        }
      }
    };

    const existing = document.getElementById("jsme-script");
    if (!existing) {
      const script = document.createElement("script");
      script.id = "jsme-script";
      script.src = "/jsme/jsme.nocache.js";
      script.async = true;
      script.onload = () => {
        if (window.jsmeOnLoad) window.jsmeOnLoad();
      };
      document.body.appendChild(script);
    } else {
      if (window.JSApplet?.JSME) {
        window.jsmeOnLoad();
      }
    }
  }, []);

  // Parse SMILES from drawer
  const handleParseSmiles = () => {
    if (jsmeRef.current) {
      const currentSmiles = jsmeRef.current.smiles();
      setSmiles(currentSmiles);
      setMatches(null); // Clear matches when structure changes
      setError(null); // Clear errors
    }
  };

  // --- Fingerprint Similarity (UPDATED for consistency) ---
  const handleFindFingerprint = async () => {
    if (!smiles) {
      console.error("No SMILES!");
      setError("Please generate a SMILES string first."); // Use setError
      return;
    }

    // --- Added Loading/Error State Handling ---
    setIsLoading(true); // Set loading
    setMatches(null);   // Clear previous matches
    setError(null);     // Clear previous errors
    // --- End Added ---

    try {
      const res = await fetch('/api/similarity', { // Endpoint for fingerprint
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles })
      });

      // --- Improved Error Handling ---
      if (!res.ok) {
        let errorData;
        try {
            // Try to parse potential JSON error body
            errorData = await res.json();
        } catch (e) {
            // If no JSON body, use the status text
            errorData = { detail: res.statusText };
        }
        console.error("Fingerprint API error:", errorData);
        // Set the error state for the UI
        setError(`Fingerprint search failed: ${errorData?.detail || errorData?.error || res.statusText}`);
        // No need to return here, finally block will handle loading state
      } else {
        // --- Success Handling (remains the same) ---
        // Expects { cids: [...] }
        const { cids } = await res.json();
        // Map to the format needed by the display logic [{cid: ...}]
        const top = cids.slice(0, 12).map(cid => ({ cid }));
        setMatches(top);
        console.log("Fingerprint candidates:", top);
      }
      // --- End Improved Error Handling ---

    } catch (e) {
      console.error("Fetch failed (Fingerprint):", e);
      // Set error state for fetch errors
      setError("Failed to fetch fingerprint similarity. Check network or server logs.");
    } finally {
      // --- Added Finally Block ---
      // Ensure loading state is turned off whether request succeeded or failed
      setIsLoading(false);
      // --- End Added ---
    }
  };

  // --- NEW: Structure Similarity (SMILES based Backend Call) ---
  const handleFindSmilesSimilarity = async () => {
    if (!smiles) {
      console.error("No SMILES!");
      setError("Please generate a SMILES string first.");
      return;
    }

    setIsLoading(true); // Set loading
    setMatches(null);   // Clear previous matches
    setError(null);     // Clear previous errors

    try {
      // Fetch from the NEW backend endpoint
      const res = await fetch('/api/similarity_smiles', { // *** NEW ENDPOINT ***
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles }) // Send the current SMILES
      });

      if (!res.ok) {
        // Try to parse error JSON, otherwise use status text
        let errorData;
        try {
            errorData = await res.json();
        } catch (e) {
            errorData = { detail: res.statusText };
        }
        console.error("SMILES Similarity API error:", errorData);
        // Use errorData.error if available from Python script, otherwise detail/statusText
        setError(`SMILES similarity search failed: ${errorData?.error || errorData?.detail || res.statusText}`);
        return; // Important: Return early on error
      }

      // Expects { results: [{cid: ..., similarity: ...}, ...] } from the new Python script
      const { results } = await res.json();
      setMatches(results); // Store the array [{cid:..., similarity:...}] directly
      console.log("SMILES Similarity candidates:", results);

    } catch (e) {
      console.error("Fetch failed (SMILES Similarity):", e);
      setError("Failed to fetch SMILES similarity. Check network or server logs.");
    } finally {
      setIsLoading(false); // Clear loading regardless of success/error
    }
  };


  return (
    <>
      <div className="chem-drawer-container" id="chem-drawer-container">
        <h1>JSME Chemical Drawer</h1>
        <div id={containerId} className="chem-canvas" />
        <div className="button-container">
          <button onClick={handleParseSmiles} className="parse-btn">
            Parse to SMILES
          </button>
          {/* Removed client-side matching button */}
          {/* <button onClick={handleFindSimilarity} className="parse-btn">
            Find Similarity For Matching
          </button> */}
          <button onClick={handleFindFingerprint} className="parse-btn" disabled={isLoading}>
            Find by Fingerprint
          </button>
          {/* --- NEW BUTTON --- */}
          <button onClick={handleFindSmilesSimilarity} className="parse-btn" disabled={isLoading}>
            Find by Structure (SMILES)
          </button>
        </div>

        {/* Display Loading Indicator */}
        {isLoading && <div className="loading-indicator">Searching...</div>}

        {/* Display Error Message */}
        {error && <div className="error-message">Error: {error}</div>}


        {smiles && !isLoading && ( // Show SMILES only when not loading
          <div className="match-results-container1">
            <p><strong>SMILES Output:</strong> {smiles}</p>
          </div>
        )}

        {/* Updated Results Display */}
        {matches !== null && !isLoading && ( // Show matches only when not loading
          <div className="matches-section" style={{ position: 'relative', zIndex: 1001 }}>
            <div className="match-results-container2">
              <h3>Matching Molecules {matches.length > 0 ? `(${matches.length} found)` : ''}</h3>
              <div style={{
                display: "flex",
                flexWrap: "wrap",
                justifyContent: "center"
              }}>
                {matches.length > 0 ? (
                  // Map over matches, check if similarity score exists
                  matches.map((match) => (
                    <a
                      key={match.cid} // Use CID as key
                      href={`/compound/${match.cid}`} // Link remains the same
                      target="_blank"
                      rel="noopener noreferrer"
                      style={{
                        textDecoration: "none",
                        color: "inherit",
                        border: "1px solid #ccc",
                        padding: "0.5rem 1rem",
                        margin: "0.5rem",
                        borderRadius: "4px",
                        display: "inline-block",
                        cursor: "pointer",
                        textAlign: "center" // Center text
                      }}
                    >
                      {/* Display CID and Similarity if available */}
                      <div>CID: {match.cid}</div>
                      {match.similarity !== undefined && (
                         <div style={{ fontSize: '0.8em', marginTop: '0.2em' }}>
                           Similarity: {match.similarity.toFixed(3)}
                         </div>
                      )}
                    </a>
                  ))
                ) : (
                  // Message when no matches are found
                  <div style={{
                    border: "1px solid #ccc",
                    padding: "0.5rem 1rem",
                    margin: "0.5rem",
                    borderRadius: "4px"
                  }}>
                    No similar structures found for the given criteria.
                  </div>
                )}
              </div>
            </div>
          </div>
        )}
      </div>
    </>
  );
}
