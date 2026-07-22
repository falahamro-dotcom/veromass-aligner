@echo off
REM ─────────────────────────────────────────────────────────────────────────────
REM  VeroMass Aligner  v1.7.0
REM  VeroMass / MoleculeID Platform — Standalone Utility
REM  Double-click this file to launch VeroMass Aligner.
REM ─────────────────────────────────────────────────────────────────────────────

cd /d "%~dp0"
python VeroMass_Aligner.py

REM  If Python is not on PATH, try the Microsoft Store Python location:
REM  "%LOCALAPPDATA%\Microsoft\WindowsApps\python.exe" VeroMass_Aligner.py

pause
