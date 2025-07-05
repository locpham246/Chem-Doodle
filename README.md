#**Chem-Doodle: Chemical Structure Drawing & Similarity Search Web App.**

This project is a web application for drawing chemical structures and searching for similar compounds using two approaches:</br>

Fingerprint-based Similarity.</br>
SMILES-based Structure Similarity.</br>
JSME Chemical Editor Integration, so users can draw chemical structures in the browser using the JSME molecular editor.</br>

# **Key Features:** </br>
**SMILES Parsing:** </br>
The drawn structure can be converted to a SMILES string for further processing.
**Fingerprint Similarity:** </br>
Submits the SMILES to a backend (Node.js + Python), which executes a fingerprint similarity search against a compound database.</br>
Results are returned as a ranked list of candidate compound IDs (CIDs).</br>
**SMILES Similarity:** </br>
Submits the SMILES to a backend endpoint, which runs a Python script to compute structure similarity.</br>
Results include both compound IDs and similarity scores.</br>
**Modern React Frontend:** </br>
Built with React and Vite for fast development, with clear state handling for loading, errors, and results. </br>
**Error Handling:** </br>
Both client and server provide user-friendly error messages and handle failed backend calls gracefully.</br>

#**How It Works** </br>
**Draw Structure:** </br>
Use the in-browser editor to sketch a molecule. 

**Parse to SMILES:** </br>
Convert the drawing into a SMILES string.

**Run Similarity Search:** </br>
Choose either fingerprint or SMILES-based similarity.</br>
The app sends your SMILES to the backend, which invokes Python scripts to perform the actual search.</br>
Results are displayed back in the UI. </br>

#**Tech Stack:** </br>
**Frontend:** </br>
React, Vite, JSME (JavaScript Molecular Editor). </br>
**Backend:** </br>
Node.js/Express server.</br>
Python scripts for chemical similarity searching.
