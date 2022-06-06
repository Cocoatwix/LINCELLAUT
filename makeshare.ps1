<#
A simple Powershell script for automating the build process
for the .o and .so files needed for ORBITVIS

May 5, 2022
#>

<# 
The following resources were used as a reference:
https://www.prodjim.com/how-to-comment-in-powershell
https://linuxhint.com/create-edit-text-files-powershell/
https://devblogs.microsoft.com/scripting/powertip-new-lines-with-powershell/
https://docs.microsoft.com/en-us/powershell/scripting/learn/deep-dives/everything-about-if?view=powershell-7.2
https://docs.microsoft.com/en-us/previous-versions/windows/it-pro/windows-powershell-1.0/ee177015(v=technet.10)
https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.management/copy-item?view=powershell-7.2
https://devblogs.microsoft.com/scripting/learn-four-ways-to-use-powershell-to-create-folders/
https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_variables?view=powershell-7.2
#>

$OBJPATH = "to-orbitvis/objects"

#Create required directories and files
if (!(Test-Path to-orbitvis))
{
	#mkdir to-orbitvis
	New-Item -Path to-orbitvis -ItemType directory
}

if (!(Test-Path to-orbitvis/objects))
{
	New-Item -Path $OBJPATH -ItemType directory
}

if (!(Test-Path to-orbitvis/README.txt))
{
	New-Item to-orbitvis/README.txt
}

#I should add a way to check if these commands execute correctly
#I should also add a way to check whether these files need to be created
# (so we only compile them if needed)
#Compile required object files
gcc -Wall -Wextra -pedantic -c -o $OBJPATH/modular.o libraries/modular.c
gcc -Wall -Wextra -pedantic -c -o $OBJPATH/linalg.o libraries/linalg.c
gcc -Wall -Wextra -pedantic -c -o $OBJPATH/cycles.o libraries/cycles.c
gcc -Wall -Wextra -pedantic -c -o $OBJPATH/factors.o libraries/factors.c

gcc -fPIC -shared -o $OBJPATH/orbitvis.so libraries/orbitvis.c $OBJPATH/modular.o $OBJPATH/linalg.o $OBJPATH/cycles.o $OBJPATH/factors.o

Set-Content to-orbitvis/README.txt "This folder contains the .o and .so files needed to run ORBITVIS. Simply drag this folder into your ORBITVIS directory, then specify where ORBITVIS can find the .o and .so files in the .config file. `r`n`r`nBy default, the files were placed in the objects subfolder, so the objects key in the .config file should probably have to-orbitvis/objects as its value."