import React, { useEffect, useRef, useState } from "react";
import data from "./assets/data.json";  // Array of objects with { cid, smiles }
import OCL from "openchemlib/full.js";

// --- Helper Functions ---

const canonicalizeSmilesUsingOCL = (smiles) => {
  const starPlaceholder = '__STAR__';
  const plusPlaceholder = '__PLUS__';
  const placeholderSmiles = smiles.replace(/\*/g, starPlaceholder).replace(/\+/g, plusPlaceholder);
  try {
    const mol = OCL.Molecule.fromSmiles(placeholderSmiles);
    mol.dearomatize();
    let canonical = mol.getCanonizedSmiles();
    canonical = canonical
      .replace(new RegExp(starPlaceholder, "g"), "*")
      .replace(new RegExp(plusPlaceholder, "g"), "+");
    return canonical;
  } catch (err) {
    console.error("Error canonicalizing SMILES:", err);
    return smiles;
  }
};

const convertSmilesToRegex = (smiles) => {
  const escapeChar = (char) => {
    if (char === "*" || char === "+") return char;
    return char.replace(/[-\/\\^$?.()|[\]{}]/g, '\\$&');
  };
  const escaped = Array.from(smiles).map(escapeChar).join("");
  const pattern = escaped.replace(/\*/g, ".*").replace(/\+/g, ".");
  return new RegExp(`^${pattern}$`, "i");
};

const ChemDrawer = () => {
  const containerId = "jsme_container";
  const jsmeRef = useRef(null);
  const [smiles, setSmiles] = useState("");
  const [matches, setMatches] = useState(null);
  const [theme, setTheme] = useState("light");

  useEffect(() => {
    const savedTheme = localStorage.getItem("theme") || "light";
    setTheme(savedTheme);
    document.body.setAttribute("data-theme", savedTheme);

    window.jsmeOnLoad = () => {
      console.log("jsmeOnLoad called");

      if (window.JSApplet && window.JSApplet.JSME && typeof window.JSApplet.JSME === "function") {
        try {
          const jsmeInstance = new window.JSApplet.JSME(containerId, "750px", "450px", { options: "query" });
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
      setMatches(null);
    }
  };

  const handleFindSimilarity = () => {
    if (!smiles) {
      console.error("No SMILES available!");
      return;
    }
    const canonicalDrawnSmiles = canonicalizeSmilesUsingOCL(smiles);
    console.log("Canonical Drawn SMILES:", canonicalDrawnSmiles);

    const hasWildcards = /[\*\+]/.test(canonicalDrawnSmiles);
    let smilesMatchFunction;
    if (hasWildcards) {
      const regex = convertSmilesToRegex(canonicalDrawnSmiles);
      smilesMatchFunction = (candidateSmiles) => regex.test(candidateSmiles);
      console.log("Using regex pattern:", regex);
    } else {
      smilesMatchFunction = (candidateSmiles) => candidateSmiles === canonicalDrawnSmiles;
    }

    const found = data.filter((candidate) => {
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

  const toggleTheme = () => {
    const newTheme = theme === "light" ? "dark" : "light";
    setTheme(newTheme);
    localStorage.setItem("theme", newTheme);
    document.body.setAttribute("data-theme", newTheme);
  };

  return (
    <>
      <header className="header">
        <img 
          src={theme === "light" ? "./src/assets/logo1.png" : "./src/assets/logo2.png"} 
          alt="Logo" 
          className="logo" 
        />
        <nav className="nav-links">
          <a href="/">Home</a>
          <a href="#chem-drawer-container">Drawer</a>
          <a href="/about">About Us</a>
          <a href="/contact">Contact</a>
        </nav>
        <button className="toggle-btn" onClick={toggleTheme}>
          {theme === "light" ? "Dark" : "Light"}
        </button>
      </header>

      <div className="chem-drawer-container" id="chem-drawer-container">
        <h1>JSME Chemical Drawer</h1>
        <div id={containerId} className="chem-canvas"></div>
        <div className="button-container">
          <button onClick={handleParseSmiles} className="parse-btn">
            Parse to SMILES
          </button>
          <button onClick={handleFindSimilarity} className="parse-btn" style={{ marginLeft: "1rem" }}>
            Find Similarity
          </button>
        </div>

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
                matches.map((candidate) => (
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
                    borderRadius: "4px",
                  }}
                >
                  No CID found that matches this structure.
                </div>
              )}
              
            </div>
          </div>
        )}
      </div>
    </>
  );
};

export default ChemDrawer;
