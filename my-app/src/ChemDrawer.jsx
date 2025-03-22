import React, { useEffect, useRef, useState } from 'react';

const ChemDrawer = () => {
  const containerId = "jsme_container";
  const jsmeRef = useRef(null); // Store the JSME instance here
  const [smiles, setSmiles] = useState("");

  useEffect(() => {
    // Global callback JSME calls once it loads
    window.jsmeOnLoad = () => {
      console.log("jsmeOnLoad called");
      console.log("window.JSApplet:", window.JSApplet);

      if (
        window.JSApplet && 
        window.JSApplet.JSME && 
        typeof window.JSApplet.JSME === 'function'
      ) {
        try {
          // Pass "query" in options to enable wildcard usage
          const jsmeInstance = new window.JSApplet.JSME(
            containerId, 
            "480px", 
            "400px", 
            { options: "query" }
          );
          console.log("JSME initialized.");

          // Store the instance in a ref so we can call .smiles() later
          jsmeRef.current = jsmeInstance;
        } catch (err) {
          console.error("Error initializing JSME:", err);
        }
      } else {
        console.error("Could not find JSME constructor. Check the JSME file.");
      }
    };

    // Dynamically load JSME script if not already loaded
    if (!document.getElementById("jsme-script")) {
      const script = document.createElement("script");
      script.id = "jsme-script";
      script.src = "/jsme/jsme.nocache.js"; // Must exist in public/jsme/
      script.async = true;
      script.onerror = () => {
        console.error("Failed to load JSME script from /jsme/jsme.nocache.js");
      };
      document.body.appendChild(script);
    }
  }, []);

  // Handler for "Parse to SMILES" button
  const handleParseSmiles = () => {
    if (jsmeRef.current) {
      // getSmiles() returns the SMILES string from the current drawing
      const currentSmiles = jsmeRef.current.smiles();
      setSmiles(currentSmiles);
      console.log("SMILES:", currentSmiles);
    }
  };

  return (
    <div style={{ textAlign: "center", marginTop: "2rem" }}>
      <h2>JSME Chemical Drawer</h2>

      {/* The JSME applet replaces this div */}
      <div id={containerId} style={{ margin: "0 auto" }}></div>

      {/* Button to parse the structure into SMILES */}
      <button onClick={handleParseSmiles} style={{ marginTop: "1rem" }}>
        Parse to SMILES
      </button>

      {/* Display the SMILES below */}
      {smiles && (
        <p>
          <strong>SMILES Output:</strong> {smiles}
        </p>
      )}
    </div>
  );
};

export default ChemDrawer;
