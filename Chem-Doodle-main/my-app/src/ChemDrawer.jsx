import React, { useEffect, useRef, useState } from "react";

const ChemDrawer = () => {
  const containerId = "jsme_container";
  const jsmeRef = useRef(null);
  const [smiles, setSmiles] = useState("");
  const [theme, setTheme] = useState("light");

  useEffect(() => {
    const savedTheme = localStorage.getItem("theme") || "light";
    setTheme(savedTheme);
    document.body.setAttribute("data-theme", savedTheme);

    window.jsmeOnLoad = () => {
      console.log("jsmeOnLoad called");

      if (
        window.JSApplet &&
        window.JSApplet.JSME &&
        typeof window.JSApplet.JSME === "function"
      ) {
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
    }
  };

  const toggleTheme = () => {
    const newTheme = theme === "light" ? "dark" : "light";
    setTheme(newTheme);
    localStorage.setItem("theme", newTheme);
    document.body.setAttribute("data-theme", newTheme);
  };

  return (
    <>
      {/*header */}
      <header className="header">
        <img src="./src/assets/logo1.png" alt="Logo" className="logo" />
        <button className="toggle-btn" onClick={toggleTheme}>
          {theme === "light" ? "Dark" : "Light"}
        </button>
      </header>

      {/* body */}
      <div className="chem-drawer-container">
        <h1>JSME Chemical Drawer</h1>
        <div id={containerId} className="chem-canvas"></div>
        <button onClick={handleParseSmiles} className="parse-btn">
          Parse to SMILES
          </button>

  {smiles && (
    <p>
      <strong>SMILES Output:</strong> {smiles}
    </p>
  )}
</div>

    </>
  );
};

export default ChemDrawer;
