# ProjetLong
PROJET PYTHON
==================================
Author
-------
Jaysen SAWMYNADEN


Python 2.7
---------------------------------------------

Files includ
---------------
comparative_modeling_Jaysen.py:
Script builds backbone and peptidic bond of a sequence protein based to a protein template. 
It computes also the RMSD between alpha carbons before gap and first and last carbon alpha of the fragment of interest.
The script needs:
		_A Pir file wich countains alignment
		_the ID of the protein request
		_the ID of the template protein
		_One PDB file which contains the 3D structure
It needs also the TOP500 folder in the same directory where the script is launched !
Output
---------------
	PDB file
	log file: all the data about the fragments


Modules for Python 2.7
--------------------------
sqlite3
sys
numpy as np
os
re
subprocess
math
Bio
time
random
Bio


How to use this program
----------------------------------
python comparative_modeling_Jaysen.py glob.pir 1bin 1lh1 1bin_chainA.pdb
    Type in the Bash:
    $ python [arg1] [arg2] [arg3] [arg4] [arg5]

arg1: comparative_modeling_Jaysen.py.py
arg2: /Directory/glob.pir
arg3: string
arg4: string
arg5: PDB

