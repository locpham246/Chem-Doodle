// src/App.jsx  (or wherever your header lives)
import React, { useState, useEffect } from 'react';
import { Routes, Route, Link } from 'react-router-dom';

// Import your logo files!
import logoLight from './assets/logo1.png';
import logoDark  from './assets/logo2.png';

import ChemDrawer   from './ChemDrawer';
import CompoundPage from './CompoundPage';

export default function App() {
  const [theme, setTheme] = useState('light');
  useEffect(() => {
    const saved = localStorage.getItem('theme') || 'light';
    setTheme(saved);
    document.body.setAttribute('data-theme', saved);
  }, []);
  const toggleTheme = () => {
    const next = theme === 'light' ? 'dark' : 'light';
    setTheme(next);
    localStorage.setItem('theme', next);
    document.body.setAttribute('data-theme', next);
  };

  return (
    <>
      <header className="header" style={{ zIndex: 10 }}>
        {/* Wrap logo in a Link so clicking it goes home */}
        <Link to="/">
          <img
            src={theme === 'light' ? logoLight : logoDark}
            alt="Molecule Map"
            className="logo"
          />
        </Link>
        <nav className="nav-links">
          <Link to="/">Home</Link>
          <Link to="/">Drawer</Link>
          <Link to="/about">About Us</Link>
          <Link to="/contact">Contact</Link>
        </nav>
        <button className="toggle-btn" onClick={toggleTheme}>
          {theme === 'light' ? 'Dark' : 'Light'}
        </button>
      </header>

      <div style={{ paddingTop: '80px' }}>
        <Routes>
          <Route path="/" element={<ChemDrawer />} />
          <Route path="/compound/:cid" element={<CompoundPage />} />
        </Routes>
      </div>
    </>
  );
}
