import React, { useEffect, useState } from "react";
import { useParams } from "react-router-dom";
import "./App.css"; // so we pick up .chem-drawer-container, header, etc.

export default function CompoundPage() {
  const { cid } = useParams();
  const [info, setInfo] = useState(null);
  const [error, setError] = useState(null);

  useEffect(() => {
    async function fetchData() {
      try {
        // Fetch IUPACName & SMILES
        const propUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/IUPACName,CanonicalSMILES/JSON`;
        const propRes = await fetch(propUrl);
        const propJson = await propRes.json();
        const props = propJson.PropertyTable?.Properties?.[0] || {};

        // Fetch synonyms (we’ll use the first one as “name”)
        const synUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/synonyms/JSON`;
        const synRes = await fetch(synUrl);
        const synJson = await synRes.json();
        const syns = synJson.InformationList?.Information?.[0]?.Synonym || [];

        setInfo({
          name: syns[0] || `CID ${cid}`,
          iupac: props.IUPACName,
          smiles: props.CanonicalSMILES,
          synonyms: syns,
        });
      } catch (e) {
        console.error(e);
        setError("Could not fetch PubChem data.");
      }
    }
    fetchData();
  }, [cid]);

  return (
    <>
      <header className="header">
        {/* keep your nav/logo here if you’d like site-wide nav */}
      </header>

      <div className="chem-drawer-container">
        <h1 style={{ textAlign: "center", marginTop: 0 }}>
          {info?.name || `CID ${cid}`}
        </h1>

        {error && <p style={{ color: "red", textAlign: "center" }}>{error}</p>}
        {!info && !error && <p style={{ textAlign: "center" }}>Loading…</p>}

        {info && (
          <div style={{ textAlign: "center", padding: "1rem" }}>
            <img
              src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG?image_size=300x300`}
              alt={info.name}
              style={{ margin: "1rem auto", maxWidth: "300px" }}
            />
            <p>
              <strong>IUPAC Name:</strong> {info.iupac}
            </p>
            <p>
              <strong>SMILES:</strong> {info.smiles}
            </p>
            <p>
              <strong>Synonyms:</strong>{" "}
              {info.synonyms.slice(0, 5).join(", ")}
              {info.synonyms.length > 5 ? " …" : ""}
            </p>
          </div>
        )}

        <a
          href="/"
          onClick={() => window.location.reload()}
          style={{
            position: "absolute",
            bottom: "1rem",
            left: "1rem",
            textDecoration: "none",
            padding: "0.5rem 1rem",
            border: "1px solid #ccc",
            borderRadius: "4px",
            background: "#f0f0f0",
          }}
        >
          ← Back to Home
        </a>
      </div>
    </>
  );
}
