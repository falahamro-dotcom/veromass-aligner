@echo off
REM ─────────────────────────────────────────────────────────────────────────────
REM  VeroMass Aligner
REM  VeroMass / MoleculeID Platform — Standalone Utility
REM
REM  This launcher always shows a black console window for its whole
REM  run — that's inherent to .bat files (they run inside a real console
REM  host), not fixable from in here. For normal use, double-click
REM  Start_VeroMass_Aligner_Windows.vbs instead — it launches the exact
REM  same thing with zero console flash. Use THIS .bat only if that one
REM  doesn't work and you need to see real error output.
REM ─────────────────────────────────────────────────────────────────────────────

cd /d "%~dp0"
python VeroMass_Aligner.py

REM  If Python is not on PATH, try the Microsoft Store Python location:
REM  "%LOCALAPPDATA%\Microsoft\WindowsApps\python.exe" VeroMass_Aligner.py

pause
