// src/firebase.js
import { initializeApp } from "firebase/app";
import { getAuth } from "firebase/auth";
import { getFirestore } from "firebase/firestore";

const firebaseConfig = {
    apiKey: "AIzaSyBVvoRcsXvLOWj5xMwGrjZN6HQ9qi7tleA",
    authDomain: "moleculemap.firebaseapp.com",
    projectId: "moleculemap",
    storageBucket: "moleculemap.firebasestorage.app",
    messagingSenderId: "593696250136",
    appId: "1:593696250136:web:0631c9527bcebbae57e66d",
    measurementId: "G-5Q4MHT2JLE"
  };

const app = initializeApp(firebaseConfig);
export const auth = getAuth(app);
export const db = getFirestore(app);
