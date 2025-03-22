// src/App.jsx
import React from 'react';
import ChemDrawer from './ChemDrawer';
import './index.css';

function App() {
  return (
    <div style={{ textAlign: 'center', marginTop: '2rem' }}>
      <h1>Chemical Drawer Demo</h1>
      <p>Draw your chemical structure below:</p>
      <ChemDrawer />
    </div>
  );
}

export default App;