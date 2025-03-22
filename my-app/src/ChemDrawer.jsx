// src/ChemDrawer.jsx
import React, { useEffect } from 'react';

const ChemDrawer = () => {
  const containerId = "jsme_container";

  useEffect(() => {
    // Define the global callback that JSME will call when loaded.
    window.jsmeOnLoad = () => {
      console.log("jsmeOnLoad called");
      console.log("window.JSApplet:", window.JSApplet);

      if (window.JSApplet && typeof window.JSApplet.JSME === 'function') {
        try {
          // Pass an options object as the fourth parameter.
          new window.JSApplet.JSME(containerId, "380px", "340px", { theme: "light" });
          console.log("JSME initialized successfully via window.JSApplet.JSME.");
        } catch (err) {
          console.error("Error while initializing JSME:", err);
        }
      } else {
        console.error("JSME constructor not found in window.JSApplet. Please verify the JSME file.");
      }
    };

    if (!document.getElementById("jsme-script")) {
      const script = document.createElement('script');
      script.id = "jsme-script";
      script.src = "/jsme/jsme.nocache.js";
      script.async = true;
      script.onerror = () => {
        console.error("Failed to load JSME script from /jsme/jsme.nocache.js");
      };
      document.body.appendChild(script);
    }
  }, []);

  return <div id={containerId}></div>;
};

export default ChemDrawer;
