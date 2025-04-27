// src/ChemDrawer.jsx
import React, { useEffect, useRef, useState } from "react";
import data from "./assets/data.json";  // Array of objects with { cid, smiles }
import OCL from "openchemlib/full.js";

const canonicalizeSmilesUsingOCL = (smiles) => {
  const starPlaceholder = "__STAR__", plusPlaceholder = "__PLUS__";
  const placeholderSmiles = smiles
    .replace(/\*/g, starPlaceholder)
    .replace(/\+/g, plusPlaceholder);
  try {
    const mol = OCL.Molecule.fromSmiles(placeholderSmiles);
    mol.dearomatize();
    let canonical = mol.getCanonizedSmiles();
    return canonical
      .replace(new RegExp(starPlaceholder, "g"), "*")
      .replace(new RegExp(plusPlaceholder, "g"), "+");
  } catch {
    return smiles;
  }
};

const convertSmilesToRegex = (smiles) => {
  const escapeChar = (c) =>
    c === "*" || c === "+" ? c : c.replace(/[-/\\^$?.()|[\]{}]/g, "\\$&");
  const escaped = Array.from(smiles).map(escapeChar).join("");
  const pattern = escaped.replace(/\*/g, ".*").replace(/\+/g, ".");
  return new RegExp(`^${pattern}$`, "i");
};

export default function ChemDrawer() {
  const containerId = "jsme_container";
  const jsmeRef = useRef(null);
  const [smiles, setSmiles] = useState("");
  const [matches, setMatches] = useState(null);
  const [, setTheme] = useState("light");

  useEffect(() => {
    const saved = localStorage.getItem("theme") || "light";
    setTheme(saved);
    document.body.setAttribute("data-theme", saved);

    window.jsmeOnLoad = () => {
      if (!jsmeRef.current && window.JSApplet?.JSME) {
        try {
          jsmeRef.current = new window.JSApplet.JSME(
            containerId, "750px", "450px", { options: "query" }
          );
        } catch (err) {
          console.error("JSME init error:", err);
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

  const handleParseSmiles = () => {
    if (jsmeRef.current) {
      setSmiles(jsmeRef.current.smiles());
      setMatches(null);
    }
  };

  const handleFindSimilarity = () => {
    if (!smiles) return console.error("No SMILES!");
    const canon = canonicalizeSmilesUsingOCL(smiles);
    const matcher = /[*+]/.test(canon)
      ? convertSmilesToRegex(canon)
      : canon;
    const found = data.filter(d => {
      const cand = canonicalizeSmilesUsingOCL(d.smiles);
      return matcher instanceof RegExp ? matcher.test(cand) : cand === matcher;
    });
    const limitedResults = found.slice(0, 12); //limit 12 datas

    setMatches(limitedResults);
    console.log("Found candidates:", limitedResults);
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
          <button onClick={handleFindSimilarity} className="parse-btn">
            Find Similarity For Matching
          </button>
          <button onClick={handleFindSimilarity} className="parse-btn">
            Find Similarity For Fingerprint
          </button>
        </div>

        {smiles && (
          <div className="match-results-container1">
            <p>
              <strong>SMILES Output:</strong> {smiles}
            </p>
          </div>
        )}

        {matches !== null && (
            <div
            className="matches-section"
            style={{ position: 'relative', zIndex: 1001 }}
          >
            <div className="match-results-container2">
              <h3>Matching Molecules</h3>
              <div style={{
                display: "flex",
                flexWrap: "wrap",
                justifyContent: "center"
              }}>
                {matches.length > 0 ? (
                  matches.map(({ cid }) => (
                    <a
                      key={cid}
                      href={`/compound/${cid}`}
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
                        cursor: "pointer"
                      }}
                    >
                      CID: {cid}
                    </a>
                  ))
                ) : (
                  <div style={{
                    border: "1px solid #ccc",
                    padding: "0.5rem 1rem",
                    margin: "0.5rem",
                    borderRadius: "4px"
                  }}>
                    No CID found that matches this structure.
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
