# LINCELLAUT
C code for studying the dynamics of linear cellular automaton under some modulus.

# Compiling Libraries for ORBITVIS
This project contains a shared library used for another program, [ORBITVIS](https://github.com/Cocoatwix/ORBITVIS). In order for that program to function properly, some code from this repo must be compiled then dragged into ORBITVIS' directory. For Windows, I've included the file "makeshare.ps1", a Powershell script for building the necessary files and putting them in one neat folder. It assumes the C compiler installed is gcc. 
You may need to bypass Powershell's default security to run the script. I find the easiest way to do that is by running `Set-ExecutionPolicy -ExecutionPolicy bypass -Scope Process`.

Check the "documentation" directory for guidance on LINCELLAUT's usage and for explanations on the inner workings of the program and its libraries.