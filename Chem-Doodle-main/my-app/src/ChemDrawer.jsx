import React, { useEffect, useRef, useState } from "react";
import { Link } from "react-router-dom";
import data from "./assets/data.json";  // Array of objects with { cid, smiles }
import OCL from "openchemlib/full.js";

// --- Helper Functions ---

const canonicalizeSmilesUsingOCL = (smiles) => {
  const starPlaceholder = "__STAR__";
  const plusPlaceholder = "__PLUS__";
  const placeholderSmiles = smiles
    .replace(/\*/g, starPlaceholder)
    .replace(/\+/g, plusPlaceholder);
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
    return char.replace(/[-\/\\^$?.()|[\]{}]/g, "\\$&");
  };
  const escaped = Array.from(smiles).map(escapeChar).join("");
  const pattern = escaped.replace(/\*/g, ".*").replace(/\+/g, ".");
  return new RegExp(`^${pattern}$`, "i");
};

export default function ChemDrawer() {
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
      if (
        window.JSApplet &&
        window.JSApplet.JSME &&
        typeof window.JSApplet.JSME === "function"
      ) {
        try {
          const instance = new window.JSApplet.JSME(
            containerId,
            "750px",
            "450px",
            { options: "query" }
          );
          jsmeRef.current = instance;
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
      script.onerror = () =>
        console.error("Failed to load JSME script from /jsme/jsme.nocache.js");
      document.body.appendChild(script);
    }
  }, []);

  const handleParseSmiles = () => {
    if (jsmeRef.current) {
      const s = jsmeRef.current.smiles();
      setSmiles(s);
      setMatches(null);
    }
  };

  const handleFindSimilarity = () => {
    if (!smiles) {
      console.error("No SMILES available!");
      return;
    }
    const canonical = canonicalizeSmilesUsingOCL(smiles);
    const wildcard = /[\*\+]/.test(canonical);
    const matcher = wildcard
      ? (() => {
          const re = convertSmilesToRegex(canonical);
          return (cand) => re.test(cand);
        })()
      : (cand) => cand === canonical;

    const found = data.filter((d) => {
      let c;
      try {
        c = canonicalizeSmilesUsingOCL(d.smiles);
      } catch {
        c = d.smiles;
      }
      return matcher(c);
    });
    setMatches(found);
  };

  const toggleTheme = () => {
    const next = theme === "light" ? "dark" : "light";
    setTheme(next);
    localStorage.setItem("theme", next);
    document.body.setAttribute("data-theme", next);
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
          <Link to="/">Home</Link>
          <Link to="#chem-drawer-container">Drawer</Link>
          <Link to="/about">About Us</Link>
          <Link to="/contact">Contact</Link>
        </nav>
        <button className="toggle-btn" onClick={toggleTheme}>
          {theme === "light" ? "Dark" : "Light"}
        </button>
      </header>

      <div className="chem-drawer-container" id="chem-drawer-container">
        <h1>JSME Chemical Drawer</h1>
        <div id={containerId} className="chem-canvas" />
        <div className="button-container">
          <button onClick={handleParseSmiles} className="parse-btn">
            Parse to SMILES
          </button>
          <button onClick={handleFindSimilarity} className="parse-btn">
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
            <div
              style={{
                display: "flex",
                flexWrap: "wrap",
                justifyContent: "center",
              }}
            >
              {matches.length > 0 ? (
                matches.map((c) => (
                  <Link
                    key={c.cid}
                    to={`/compound/${c.cid}`}
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
                    }}
                  >
                    CID: {c.cid}
                  </Link>
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
}
