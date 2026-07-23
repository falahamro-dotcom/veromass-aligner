' VeroMass Aligner — silent standalone launcher.
'
' Double-click this file to launch VeroMass Aligner with NO black
' console/PowerShell window flashing. The old .bat launcher always
' showed one, because a .bat file runs inside a real console host
' (cmd.exe) for its whole lifetime — there's no way to avoid that from
' inside a .bat itself. This VBScript instead calls pythonw.exe (the
' windowless Python build — same one veromass-bridge/launcher.py uses
' for exactly this reason) via WScript.Shell.Run with a hidden window
' style, and wscript.exe (what runs a double-clicked .vbs by default)
' has no console of its own either — so nothing ever flashes.
'
' If this doesn't launch anything (e.g. pythonw isn't on PATH), use
' Start_VeroMass_Aligner_Windows.bat instead — it shows real error
' output in a console window, which this silent version deliberately
' does not.

Set fso = CreateObject("Scripting.FileSystemObject")
scriptDir = fso.GetParentFolderName(WScript.ScriptFullName)

Set shell = CreateObject("WScript.Shell")
shell.CurrentDirectory = scriptDir
shell.Run "pythonw """ & scriptDir & "\VeroMass_Aligner.py""", 0, False
