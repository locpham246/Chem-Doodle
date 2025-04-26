// src/App.jsx
import React, { useState, useEffect } from 'react';
import { Routes, Route, Link } from 'react-router-dom';
import ChemDrawer from './ChemDrawer';
import CompoundPage from './CompoundPage';

function App() {
  // Move theme state here
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
      <header className="header">
        <img src={theme === 'light' ? './src/assets/logo1.png' : './src/assets/logo2.png'} alt="Logo" className="logo" />
        <nav className="nav-links">
          <Link to="/">Home</Link>
          <Link to="/">Drawer</Link>
          <Link to="/about">About Us</Link>
          <Link to="/contact">Contact</Link>
        </nav>
        {/* here’s your theme toggle, now in the shared header */}
        <button className="toggle-btn" onClick={toggleTheme}>
          {theme === 'light' ? 'Dark' : 'Light'}
        </button>
      </header>

      <div style={{ paddingTop: '60px' }}>
        <Routes>
          <Route path="/" element={<ChemDrawer />} />
          <Route path="/compound/:cid" element={<CompoundPage />} />
        </Routes>
      </div>
    </>
  );
}

export default App;
