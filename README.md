Chem-Doodle: Chemical Structure Drawing & Similarity Search Web App.

This project is a web application for drawing chemical structures and searching for similar compounds using two approaches:

Fingerprint-based Similarity
SMILES-based Structure Similarity
Key Features
JSME Chemical Editor Integration:
Users can draw chemical structures in the browser using the JSME molecular editor.

SMILES Parsing:
The drawn structure can be converted to a SMILES string for further processing.

Similarity Search (Backend-Powered):

Fingerprint Similarity:
Submits the SMILES to a backend (Node.js + Python), which executes a fingerprint similarity search against a compound database.
Results are returned as a ranked list of candidate compound IDs (CIDs).
SMILES Similarity:
Submits the SMILES to a backend endpoint, which runs a Python script to compute structure similarity.
Results include both compound IDs and similarity scores.
Modern React Frontend:
Built with React and Vite for fast development, with clear state handling for loading, errors, and results.

Robust Error Handling:
Both client and server provide user-friendly error messages and handle failed backend calls gracefully.

How It Works
Draw Structure:
Use the in-browser editor to sketch a molecule.

Parse to SMILES:
Convert the drawing into a SMILES string.

Run Similarity Search:

Choose either fingerprint or SMILES-based similarity.
The app sends your SMILES to the backend, which invokes Python scripts to perform the actual search.
Results are displayed back in the UI.
Tech Stack
Frontend:
React, Vite, JSME (JavaScript Molecular Editor)
Backend:
Node.js/Express server
Python scripts for chemical similarity searching
