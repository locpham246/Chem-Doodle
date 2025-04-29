// src/pages/Profile.jsx
import React from 'react';
import { auth } from '../firebase';
import { signOut } from 'firebase/auth';
import { useNavigate } from 'react-router-dom';

export default function Profile() {
  const navigate = useNavigate();

  const handleLogout = async () => {
    await signOut(auth);
    navigate('/');
  };

  return (
    <div className="profile-container">
      <h2>Welcome {auth.currentUser?.email}</h2>
      <button onClick={handleLogout}>Logout</button>
    </div>
  );
}
