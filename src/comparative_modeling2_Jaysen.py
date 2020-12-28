#!/usr/bin/python
# -*- coding: utf-8 -*-
import sqlite3
import sys
import numpy as np
import os
import re
import subprocess
import math
import Bio
from Bio.PDB import *
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from copy import deepcopy
import time
from random import *

class Atom:
	""""
	Class composed:
	-Name of the type of amino acid
	-Coordonnee
	-ID
	"""
	def __init__(self, nom, x, y, z, idCA, nb_res, chain):
		self.name = nom
		self.X = x
		self.Y = y
		self.Z = z
		self.ID = idCA
		self.NbRes = nb_res
		self.Chain = chain

class CA:
	""""
	Class composed:
	-Name of the type of amino acid
	-Coordonnee of N, CA, C and O group
	"""
	def __init__(self, AA, N):
		self.Name = AA
		self.PosiN = N
		self.PosiCA = []
		self.PosiC = []
		self.PosiO = []

#Residu(id_res,name_res,

path = os.getcwd()
path = path+"/top500H/"
SupPDB = os.listdir(path)
AminoAcid = {'ALA':'A', 'ARG':'R', 'ASP':'D', 'ASN':'N', 'CYS':'C',
'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L',
'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S','THR':'T',
'TRP':'W','TYR':'Y','VAL':'V'}
start = time.time()

#################################################################
#																#
#			Function to check the arguments						#
#																#
#################################################################
# Check the arguments given
# Vérifie les arguments donné à l'interpréteur Python
def help_me():
	print"If you need help, there is a README. The value of radiusn\
\tpython comparative_modeling_Jaysen.py glob.pir 1bin 1lh1 1bin_chainA.pdb\n\n"

#CodeTemplate = sys.argv[2]#"1bin"
#CodeQuery = sys.argv[3] #"1lh1"
#Fichier = sys.argv[1] #"glob.pir"
#PDB = sys.argv[4] #'1bin_chainA.pdb'
def verification_argv():
    """
    Check the arguments given to python:
	For help you can use the option -help

    An error message will be on screen and the programme the programme'll close if:
        Number of arguments is different of (4)
        #2 argument: The file doesn't exist
        #2 argument: The files doesn't have extension .pdb
        #3 argument: The chain doesn't have the good format:
		all: for all the protein
		A: analyse only the chain A
		A-B-C: analyse the chain A B C


    Exemples:
        $ python script.py
        ERREUR: ERREUR: The program needs a PDB file and the number of points per sphere.

        $ python script.py file.pdb
        ERREUR #2 argument: The value must be an integer > 0. Why don't you choose 92,
	it's a good number. :)

        $ python script.py doesnt_exist.pdb 92
        ERREUR #3 argument: The file "doesnt_exist.pdb" can't found.


        $ python script.py fichierPDB 92
        ERREUR #4 argument: A PDB file is required
        Add an extension *.pdb in the end of the file's name if it's a PDB.
    """
    if len(sys.argv) != 5:
	help_me()
        sys.exit("ERREUR: The program needs a PDB, Pir file and the ID of the template and request protein.")

    if not os.path.exists(sys.argv[1]):
        sys.exit("ERREUR #3 argument: The file \"%s\" doesn't exists." % sys.argv[1])

    if not re.search('.((p|P)(d|D)(b|B))$', sys.argv[4]):
        sys.exit("ERREUR #4 argument: A PDB file is required.\n\
Add an extension *.pdb in the end of the file's name if it's a PDB.")

#################################################################
#				Function to launch command bash					#
#################################################################
def bash_commande(cmd):
    subprocess.Popen(cmd, shell=True, executable='/bin/bash')

#################################################################
#	I) Parse all pdb in the folder to create data base			#
#################################################################
#################################################################
#	1) Fill the table C_ALPHA 									#
#################################################################
def feed_CA(DicoCA):
	"""
	Fill the table C_ALPHA
	_DicoCA: dictionary of Atom object
	"""
	test = ""
	for elemt in DicoCA.keys():
	#print [CA.name, CA.X, CA.Y, CA.Z,CA.ID]
		cursor.execute("""
			INSERT INTO C_ALPHA(name,X,Y,Z,id_CA, Nb, chain)
			VALUES(?,?,?,?,?,?,?)
			""",
				(DicoCA[elemt].name,
				DicoCA[elemt].X,
				DicoCA[elemt].Y,
				DicoCA[elemt].Z,
				DicoCA[elemt].ID,
				DicoCA[elemt].NbRes,
				DicoCA[elemt].Chain)
		)

#################################################################
#	2) Compute the distance between CA of the same chain		#
#################################################################

def distCA(DicoCA,SizeFrag, nom, sequence):
	"""
	Compute the distance between the first C alpha and the other
	C alpha of the chain. The function fills also the tables fragment
	and dist.
	Arguments:
		_DicoCA: Dictionnary which countains Atom object
		_SizeFrag: the fragment's size
		_nom: PDB's code
		_sequence: the protein's sequence
	"""
	for i in xrange(len(DicoCA)-1):
		if i+SizeFrag < len(sequence):
			end = SizeFrag + i
		else:
			end = -1
		cursor.execute("""
		INSERT INTO fragment(id_frag, sequence, id_CA1, id_CA2, indice)
		VALUES(?,?,?,?,?)
		""",(
			'frag_'+str(DicoCA[i].ID),
			sequence[i:end],
			DicoCA[i].ID,
			DicoCA[i].ID,
			0
			)
		)
		for j in xrange(1,SizeFrag):
			if(i+j not in DicoCA.keys()):
				break
	#We fill the table fragment with the feature of the other C alpha
			cursor.execute("""
			INSERT INTO fragment(id_frag, sequence, id_CA1, id_CA2, indice)
			VALUES(?,?,?,?,?)
			""",(
				'frag_'+str(DicoCA[i].ID),
				sequence[i:end],
				DicoCA[i].ID,
				DicoCA[i+j].ID,
				j
				)
			)
	#We fill the table dist
			X = (DicoCA[i].X - DicoCA[i+j].X)**2
			Y = (DicoCA[i].Y - DicoCA[i+j].Y)**2
			Z = (DicoCA[i].Z - DicoCA[i+j].Z)**2
			cursor.execute("""
			INSERT INTO distance(dist,id_CA1,id_CA2)
			VALUES(?,?,?)
			""",(
				round(np.sqrt(X+Y+Z),2),
				DicoCA[i].ID,
				DicoCA[i+j].ID)
			)


#################################################################
#	3) Function to parse the pdb																	#
#################################################################

def parse_pdb_data_base(SupPDB, path, AminoAcid, i):
	"""
	Parse the PDB and create a dictionary of Atom objet
	composed:
	Arguments:
		_SupPDB: List of all the PDB in the folder top500H
		_path: path to the curent directory
		_AminioAcid: dictionary converts amino acide 3 letters
		to amino acid 1 letters
		_i (INTEGER): access to a specific PDB in SupPDB
	"""
	SizeFrag = 30 #taille des fragments request
	nom = path+SupPDB[i]
	sequence = ""
	chain = []
	cpt = 0	#compte le nombre de fragment généré
	DicoCA = {} #contains the list of Carbon alpha
	with open(nom, "r") as filin:
		for line in filin:
			if line[:4] == 'ATOM' and 'CA' in line[12:16].split():
				splitted_line = [
				line[:6], #'ATOM'
				line[6:11], #Atom's number
				line[12:16], #Atom's type (CA,NH,CZ...)
				line[17:20], #Amino acid's name
				line[21], #Chain's name
				line[22:26], #Residue's number
				line[30:38], #cartesian coordinates:
				line[38:46], #cartesian coordinates:
				line[46:54], #cartesian coordinates:
				line[77] #Atom's name
				]
#################################################################
#					When we are in a new chain					#
#################################################################
				if line[21] not in chain:
					if len(sequence) > 0:
						feed_CA(DicoCA)
						distCA(DicoCA,SizeFrag, SupPDB[i], sequence)
						sequence = ""
						chain = []
						cpt = 0
						DicoCA = {}
#we save the previous fragment (if it exist) before we erase it to make the fragment
					chain.append(line[21])
#If the amino acid is known
					if(splitted_line[3] in AminoAcid.keys()):
						AA = AminoAcid[splitted_line[3]]
					else:
						AA = 'X'
					sequence = AA
							#Atom(nom, x, y, z, idCA, Residue's number, chain)
					DicoCA[0] = Atom(
								AA,
								float(splitted_line[6]),
								float(splitted_line[7]),
								float(splitted_line[8]),
								SupPDB[i]+"_"+str(cpt)+"_"+splitted_line[4],
								splitted_line[5],
								splitted_line[4]
								)
#We fill the data base with the features of the C alpha
					cpt += 1
#################################################################
#		When we are in the same chain					 		#
#################################################################
				else:
#If the amino acid is known
					if(splitted_line[3] in AminoAcid.keys()):
						AA = AminoAcid[splitted_line[3]]
					else:
						AA = 'X'
					sequence += AA
							#Atom(nom, x, y, z, idCA, Residue's number, chain)
					DicoCA[cpt] = Atom(
								AA,
								float(splitted_line[6]),
								float(splitted_line[7]),
								float(splitted_line[8]),
								SupPDB[i]+"_"+str(cpt)+"_"+splitted_line[4],
								splitted_line[5],
								splitted_line[4]
								)
#We fill the data base with the features of the C alpha
					cpt += 1
	feed_CA(DicoCA)
	distCA(DicoCA,SizeFrag, SupPDB[i], sequence)



#################################################################
#	0) Create database if it doesn't exit						#
#################################################################
#def create_database():
conn = sqlite3.connect('jaysen_base.db')
print("Create database...")
conn = sqlite3.connect('jaysen_base.db')
cursor = conn.cursor()
cursor.execute("""
		CREATE TABLE IF NOT EXISTS C_ALPHA(
		name STRING,
		X DOUBLE,
		Y DOUBLE,
		Z DOUBLE,
		id_CA STRING PRIMARY KEY UNIQUE,
		Nb STRING,
		chain STRING
		)
	""")
conn.commit()
cursor.execute("""
	CREATE TABLE IF NOT EXISTS fragment(
	id_frag STRING,
	sequence STRING,
	id_CA1 REFERENCES C_ALPHA(id_CA),
	id_CA2 REFERENCES C_ALPHA(id_CA),
	indice INTEGER NOT NULL,
	PRIMARY KEY(id_frag, id_CA1, id_CA2)
	)
	""")
conn.commit()
cursor.execute("""
	CREATE TABLE IF NOT EXISTS distance(
	dist DOUBLE,
	id_CA1 REFERENCES C_ALPHA(id_CA),
	id_CA2 REFERENCES C_ALPHA(id_CA)
	)
""")
conn.commit()
for i in xrange(len(SupPDB)):
	parse_pdb_data_base(SupPDB, path, AminoAcid, i)

print "Curent time since the lunch of the software: {0:.2f} s".format(time.time() - start)
#################################################################
#	II) Parse the template and compare to the query sequence	#
#################################################################

#################################################################
#	1) Parse the template and Check miss query					#
#################################################################
def parse_template2(PDB, AminoAcid, MissCA):
	"""
Parse a PDB file and compares it with a sequence query
	Arguments:
		_PDB: PDB file
		_AminioAcid: Dictionary which converts alphabet 3
		letters into alphabets 1 letter
		_query: the query sequence in string format
	Return:
		_DicoCA: Dictionary which contains the coordinates of
		Carbon alpha
		_missing: Dictionary which contains the Carbon alpha'number
		of carbon alpha missing
		_BoxDimension: The box's dimension of the PDB
	"""
	KeyMiss = 0 #count the number of missing fragment
	cpt = 0	#count the number of C alpha met
	debut = 0
	#and the rest of the principal chain (N, CO)
	DicoCA = {} #contains the list of Carbon alpha
	PeptideCA = {}
	BoxDimension =''
	chain = []
	splitted = []
	NbResid = 0
	FlagN = False;
	FlagC = False;
	FlagO = False;
	with open(PDB, "r") as filin:
		for line in filin:
			if line[:6] == 'CRYST1':
				BoxDimension = line
			if line[:4] == 'ATOM':
				splitted = [
				line[:6], #'ATOM'
				line[6:11], #Atom's number
				line[12:16], #Atom's type (CA,NH,CZ...)
				line[17:20], #Amino acid's name
				line[21], #Chain's name
				line[22:26], #Residue's number
				line[30:38], #cartesian coordinates:
				line[38:46], #cartesian coordinates:
				line[46:54], #cartesian coordinates:
				line[77] #Atom's name
				]
		#################################################################
		#					When we are in a new chain					#
		#################################################################
				if str(line[21]) not in chain:
					chain.append(str(line[21]))
					cpt = int(line[22:26])
					debut = int(line[22:26])
					NbResid = int(line[22:26])
				if NbResid != int(line[22:26]):
					FlagN = False
					FlagC = False
					FlagO = False
					NbResid = int(line[22:26])
		#################################################################
		#		When we are in the same chain					 		#
		#################################################################
				if splitted[2].split()[0] == 'N' and FlagN is False:
					if splitted[3] in AminoAcid.keys():
						AA = AminoAcid[splitted[3]]
					else:
						AA = 'X'
						print "The amino acid {0} is unknown !".format(int(splitted[5]))
					PeptideCA[cpt - debut] = CA(AA,
												[
												float(splitted[6]),
												float(splitted[7]),
												float(splitted[8])
												]
											)
					FlagN = True
				if splitted[2].split()[0] == 'CA':
					PeptideCA[cpt - debut].PosiCA = [
													float(splitted[6]),
													float(splitted[7]),
													float(splitted[8])
													]
					if splitted[3] in AminoAcid.keys():
						AA = AminoAcid[splitted[3]]
					else:
						AminoAcid = 'X'
						print "The amino acid {0} is unknown !".format(int(splitted[5]))
					DicoCA[cpt - debut] = [
									AA,
									float(splitted[6]),
									float(splitted[7]),
									float(splitted[8]),
									]
				if splitted[2].split()[0] == 'C' and FlagC is False:
					PeptideCA[cpt - debut].PosiC = [
													float(splitted[6]),
													float(splitted[7]),
													float(splitted[8])
													]
					FlagC = True
				if splitted[2].split()[0] == 'O' and FlagO is False:
					PeptideCA[cpt - debut].PosiO = [
													float(splitted[6]),
													float(splitted[7]),
													float(splitted[8])
													]
					FlagO = True
					if MissCA[KeyMiss][0] == (cpt-debut):
						cpt = MissCA[KeyMiss][1] + 1
						KeyMiss += 1
					else:
						cpt += 1
					if KeyMiss not in MissCA.keys():
						KeyMiss = 0
	return DicoCA, PeptideCA, BoxDimension


#################################################################
def computeDist(list1, list2):
	X = (list1[1]-list2[1])**2
	Y = (list1[2]-list2[2])**2
	Z = (list1[3]-list2[3])**2
	return round(np.sqrt(X+Y+Z),2)



#################################################################
#					potential lennard jones						#
#################################################################

def CoordMissFrag(SeqGive, BoxDimension,path):
	"""
	backup the CA's coordinate of each fragment and
	storaged them in the dictionary DicoCoord
	Arguments:
		_SeqGive: List of the fragment's ID
		_BoxDimension: The box's dimension of the template's PDB
		_path: Where are localizated the PDBs
	Return
		_DicoCoord: dictionary of biopython object
		_OtherInfo: dictionary that countains the beging, end
		the chain's name and the name of the PDB where the
		fragment was from.
	"""
	DicoCoord = {}
	OtherInfo = {}
	i = 0
	for Frag in SeqGive:
		cursor.execute("""
			SELECT C_ALPHA.Nb, C_ALPHA.chain
			FROM C_ALPHA, fragment
			WHERE C_ALPHA.id_CA = fragment.id_CA1
			AND fragment.indice = 0
			AND fragment.id_frag = ?
			""",(Frag,)
		)
		resp = cursor.fetchone()
		BegingPDB = int(resp[0])
		EndPDB = BegingPDB + len(MissCA[clef][2]) -1
		NamePDB = Frag.split("_")[1]
		ChainPDB = str(resp[1])
		BackupFromPDB(NamePDB, BegingPDB, EndPDB, ChainPDB, BoxDimension,path)
		p = PDBParser()
		#Loading of the fragment's structure
		DicoCoord[i] = p.get_structure('X', 'tmpFrag.pdb')
		OtherInfo[i] = [BegingPDB, EndPDB, ChainPDB, NamePDB]
		i += 1
	return DicoCoord, OtherInfo

def TinyGap(MissFrag, query, TemplateCA, constraint, strict, template):
	beg = MissFrag[0] #Indice where begings the gap
	end = MissFrag[1] #Indice where the gap ended
	SeqMiss = query[beg-2:end+2+1]
	SizeFrag = len(SeqMiss)
	check = list(range(beg-2, end+2+1))
	check.remove(beg+1)
	elemt = TemplateCA.keys()
	FinalSize = len(template[beg-2:end+2+1])
	for elemt in check:
		if elemt not in TemplateCA.keys():
			print "{0}".format(elemt)
			print "Can't build a model for this gap"
			return [], [], len(template[beg-2:end+2+1])
	ListAA1 = TemplateCA[beg-1]
	ListAA2 = TemplateCA[beg-2]
	ListAA1B = TemplateCA[end+1]
	ListAA2B = TemplateCA[end+2]
	###########
	dist1 = computeDist(ListAA2,TemplateCA[end+1])
	SizeDist1 = len(template[beg-2:end+1+1])
	###########
	dist2 = computeDist(ListAA2,ListAA2B)
	SizeDist2 = len(template[beg-2:end+2+1])
	###########
	cursor.execute(
		"""
		SELECT SUBSTR(fragment.sequence,0, ?), fragment.id_frag
		FROM distance, fragment
		WHERE distance.dist between ? AND ?
		AND fragment.id_CA1 = distance.id_CA1
		AND (distance.id_CA2 = fragment.id_CA2
		AND fragment.indice = ? )
		AND fragment.id_CA1 IN
			(SELECT fragment.id_CA1
			FROM distance, fragment
			WHERE distance.dist between ? AND ?
			AND fragment.id_CA1 = distance.id_CA1
			AND (distance.id_CA2 = fragment.id_CA2
			AND fragment.indice = ? )
			)
			""",
				(FinalSize + 1,
				dist1 - constraint,
				dist1 - constraint,
				SizeDist1 -1,
				dist2 - constraint,
				dist2 + constraint,
				SizeDist2 - 1,
				)
	)
	SeqGive = []
	Score = []
	tmp = 0
	testcursor = cursor
	for row in cursor:
		tmp = pairwise2.align.globalxx(\
		SeqMiss,str(row[0]), score_only=1)/float(FinalSize)
		if tmp >= strict:
			SeqGive.append(str(row[1]))
			Score.append(round(tmp,2))
		elif str(row[0])[0] == SeqMiss[0] and tmp > 0.2:
			SeqGive.append(str(row[1]))
			Score.append(round(tmp,2))
	if strict == 0:
		return SeqGive, Score, FinalSize, [beg-2, end+2, template[beg-2:end+2+1], len(template[beg-2:end+2+1])]
	elif len(Score) == 0:
			#we retry with a lower strict
		return TinyGap(MissFrag, query, TemplateCA, constraint+0.05, strict-0.1, template)
	else:
		return SeqGive, Score, len(template[beg-2:end+2+1])


def SmallGap(MissFrag, query, TemplateCA, constraint, strict, template):
	"""
	Arguments:
		_MissFrag: List countains the position before (beg) and after (end)
		the gap, the sequence (between beg-2 and end+2) and its size
		_query: the sequence in amino acid
		_TemplateCA: dictionary which countains the positions and the name of
		all the CA known
		_constraint:
	"""
	beg = MissFrag[0] #Indice where begings the gap
	end = MissFrag[1] #Indice where the gap ended
	SeqMiss = query[beg-2:end+2+1]
	SizeFrag = len(SeqMiss)
	check = list(range(beg-4, beg+1))
	check += list(range(end, end+4+1))
	ListClee = TemplateCA.keys()
	for elemt in check:
		if elemt not in TemplateCA.keys():
			print "Gap is too close of an other gap..."
			print "A model less strict based on the distance between position {0} - {1} and the sequence will be generated.".format(beg-2, end+2)
			return TinyGap(MissFrag, query, TemplateCA, constraint, strict, template)
	#ListAA1 is the second CA with a known position before the gap
	ListAA1 = TemplateCA[beg-1]
	ListAA2 = TemplateCA[beg-2]
	ListAA4 = TemplateCA[beg-4]
	ListAA1B = TemplateCA[end+1]
	ListAA2B = TemplateCA[end+2]
	ListAA3B = TemplateCA[end+3]
	###########
	dist1 = computeDist(ListAA2,TemplateCA[end])
	SizeDist1 = len(template[beg-2:end+1])
	###########
	dist2 = computeDist(ListAA2,ListAA2B)
	SizeDist2 = len(template[beg-2:end+2+1])
	###########
	dist3 = computeDist(ListAA2,ListAA4)
	SizeDist3 = len(template[beg-4:beg-2+1])
	############################################
	dist1B = computeDist(ListAA2,TemplateCA[beg])
	SizeDist1B = len(template[beg-2:beg+1])
	###########
	dist2B = computeDist(TemplateCA[end],ListAA3B)
	SizeDist2B = len(template[end:end+3+1])
	###########
	FinalSize = len(template[beg-2:end+2+1])
	cursor.execute(
		"""
		SELECT SUBSTR(fragment.sequence,0, ?), fragment.id_frag
		FROM distance, fragment
		WHERE distance.dist between ? AND ?
		AND fragment.indice = ?
		AND fragment.id_CA1 = distance.id_CA1
		AND fragment.id_CA1 IN
			(SELECT distance.id_CA1
			FROM distance, fragment
			WHERE distance.dist between ? AND ?
			AND fragment.indice = ?
			AND fragment.id_CA1 = distance.id_CA1
			AND fragment.id_CA1 IN (
				SELECT distance.id_CA2
				FROM distance, fragment
				WHERE distance.dist between ? AND ?
				AND fragment.indice = ?
				AND distance.id_CA2 = fragment.id_CA2
				AND fragment.id_CA2 IN (
					SELECT distance.id_CA1
					FROM distance, fragment
					WHERE distance.dist between ? AND ?
					AND (distance.id_CA2 = fragment.id_CA2
					AND fragment.indice = ?
					AND distance.id_CA2 IN (
						SELECT distance.id_CA1
						FROM distance, fragment
						WHERE distance.dist between ? AND ?
						AND (distance.id_CA2 = fragment.id_CA2
						AND fragment.indice = ? )
							)
						)
					)
				)
			)""",(FinalSize + 1,
					dist1B -constraint,
					dist1B + constraint,
					SizeDist1B - 1,
					dist2 - constraint,
					dist2 + constraint,
					SizeDist2 - 1,
					dist3 - constraint,
					dist3 + constraint,
					SizeDist3 - 1,
					dist1 - constraint,
					dist1 + constraint,
					SizeDist1 - 1,
					dist2B - constraint,
					dist2B + constraint,
					SizeDist2B - 1,
				)
	)
	SeqGive = []
	Score = []
	tmp = 0
	for row in cursor:
		tmp = pairwise2.align.globalxx(\
		SeqMiss,str(row[0]), score_only=1)/float(FinalSize)
		if tmp >= strict:
			SeqGive.append(str(row[1]))
			Score.append(round(tmp,2))
		elif str(row[0])[0] == SeqMiss[0] and tmp > 0.2:
			SeqGive.append(str(row[1]))
			Score.append(round(tmp,2))
	if strict == 0:
		return SeqGive, Score, FinalSize, [beg-2, end+2, template[beg-2:end+2+1], len(template[beg-2:end+2+1])]
	elif len(Score) == 0:
			#we retry with a lower strict
		return SmallGap(MissFrag, query, TemplateCA, constraint+0.05, strict-0.1, template)
	else:
		return SeqGive, Score, len(template[beg-2:end+2+1])


def OutPutPDB(indice,Coord,BoxDimension):
	cpt = 0
	BB = 1
	nb = 1
	NomFichier = "tmp.pdb"
	Fichier = open(NomFichier, "w")
	Fichier.write("REMARK    GENERATED BY JAYSEN\n")
	Fichier.write("TITLE     structure based from pdb t=   "+str(nb)+"\n")
	Fichier.write("REMARK    THIS IS A SIMULATION BOX\n")
	Fichier.write(str(BoxDimension))
	Fichier.write("MODEL "+str(nb)+"\n")
	for ligne in Coord:
		X = round(ligne[0],2)
		Y = round(ligne[1],2)
		Z = round(ligne[2],2)
		cpt += 1
		Fichier.write("%-6s%5d "%('ATOM', cpt))
		Fichier.write("%4s%1s%3s "%(str("CA"), " ",str("C")[:3]))
		Fichier.write("%1s%4d%1s   "%('A', cpt, " ") )
		Fichier.write("%8.3f%8.3f%8.3f%6.2f%6.2f"%(X,Y,Z, 1.00 ,42.00) )
		Fichier.write("          %2s%2s\n"%(str(cpt), 0))
	Fichier.write("TER\n")
	Fichier.close()
	return MDAnalysis.Universe("tmp.pdb")


def	RotatePosition(CoordRotated, indice,\
 SplitTemp, clef, TemplatePDB, PeptideCA, OtherInfo):
 	"""
	Rotate the carbon alpha chain of the fragment. The
	rotation is made to minimize the RMSD between the
	first and last fragment's C alpha and the C alpha
	before the gap.
	Arguments:
		_CoordRotated: dictionary countains the coordinates
		of the fragment's C alpha
		_indice a key value for the dictionary CoordRotated
		_SplitTemp: The list of the template's sequence
		between each gap
		_clef: indice for the list SplitTemp
		_TemplateCA: MDanalysis objet
		_PeptideCA: dictionary
		_OtherInfo: dictionary that countains the beging,
		end, the chain's name and the name of the PDB where
		the fragment was from.
	"""
	#gap: where are the gap in the template's structure
	beging = ""
	for i in xrange(len(SplitTemp[0:clef+1])):
		beging += SplitTemp[i]
	end = len(beging) + 1 + 2
	beging = len(beging) - 2
	gap = [beging, end]
	#Loading of the fragment's structure
	TmpModel = DicoCoord[indice][0]
	TmpChain = TmpModel[OtherInfo[indice][2]]
	TmpResBeg = TmpChain[(' ',OtherInfo[indice][0],' ')]
	TmpResEnd = TmpChain[(' ',OtherInfo[indice][1],' ')]
	TmpAtomBeg = TmpResBeg['CA']
	TmpAtomEnd = TmpResEnd['CA']
	#TmpAtomBeg.get_coord()
	#TmpAtomEnd.get_coord()
	#Loading of the template's structure
	#We add +1 for all indie because biopython beging to 1.
	TemplateBeg = TemplatePDB[0]['A'][(' ',gap[0],' ')]['CA']
	TemplateEnd = TemplatePDB[0]['A'][(' ',gap[1],' ')]['CA']
	#TemplateBeg.get_coord()
	#TemplateEnd.get_coord()
	sup = Superimposer()
	sup.set_atoms([TemplateBeg,TemplateEnd],[TmpAtomBeg,TmpAtomEnd])
	#print sup.rotran
	#print sup.rms
	tmp = TmpModel[OtherInfo[indice][2]]
	sup.apply(tmp)
	#for i in xrange(len(TmpChain)):
	#	print tmp[(' ',OtherInfo[indice][0]+i,' ')]['N'].get_coord()
	#print 'value of RMSD: {0}'.format(sup.rms)
	return tmp, sup.rms

#################################################################
#					potential lennard jones						#
#################################################################
def lennardJones(TemplateCA, Xx, Yy, Zz, MissFrag):
	LJ = 0
	clach = False
	for atom in TemplateCA.keys():
#Don't compute the potential where the fragment will be
		if atom >= MissFrag[0]-2 and atom <= MissFrag[1]+2:
			continue
		X1 = TemplateCA[atom][1]
		Y1 = TemplateCA[atom][2]
		Z1 = TemplateCA[atom][3]
		dist = round(np.sqrt(
			np.square(X1-Xx) +
			np.square(Y1-Yy) +
			np.square(Z1-Zz)
			),2)
		if dist <= 1:
			clach = True
			#r is equal to 1
			LJ += 4*0.0738*((4.315**12)-2*(4.315**6))
		elif dist <= 6:
			LJ += 4*0.0738*(((4.315/dist)**12)-2*((4.315/dist)**6))
	if clach == True:
		print 'steric clash !'
	return LJ



def BestFragment(DicoCoord, MissFrag, OtherInfo,\
 TemplateCA, PeptideCA, SplitTemp, clef, TemplatePDB):
	FinalFrag = []
	beg = MissFrag[0]
	end = MissFrag[1]
	CoordRotated = {}
	MatrixRotate = {}
	RMSD = {}
	print "Determination of the best framgent..."
	for indice in xrange(len(DicoCoord)):
		CoordRotated[indice], RMSD[indice] = RotatePosition(CoordRotated,
												indice,
												SplitTemp,
												clef,
												TemplatePDB,
												PeptideCA,
												OtherInfo)
	Coord = CoordRotated
	PLJ = []
	#We don't take the first and the last atoms because
	#of they will be used to replace atoms before and end
	#the gap
	minRMSD = []
	print "compute the potential Lennard Jones between Calpha"
	message = "{0}\t{1}\t{2}".format("PDB","RMSD","potential Lennard Jones")
	print message
	build_info(message)
	for frag in xrange(len(Coord)):
		LJ = 0
		if RMSD[frag] is None:
			#print "The fragment is not fitted"
			PLJ.append(9999999999)
			continue
		for j in xrange(len(Coord[frag])):
			TmpChain = DicoCoord[frag][0][OtherInfo[frag][2]]
			TmpRes = TmpChain[(' ',OtherInfo[frag][0]+j,' ')]
			TmpAtom = TmpRes['CA'].get_coord()
			Xx = round(TmpAtom[0],2)
			Yy = round(TmpAtom[1],2)
			Zz = round(TmpAtom[2],2)
			LJ += lennardJones(TemplateCA, Xx, Yy, Zz, MissFrag)
		message = "{0}\t{1}\t{2}".format(OtherInfo[frag][3],\
		round(RMSD[frag],2),round(LJ,2))
		print message
		build_info(message)
		PLJ.append(round(LJ,2))
		minRMSD.append(round(RMSD[frag],2))
	PLJ = np.array(PLJ)
	minRMSD = np.array(minRMSD)
	#Look only fragments which have a RMSD below to 1 A
	if len(np.where(minRMSD <= 1)[0]) > 0:
		IndMinRMSD = np.where(minRMSD <= 1)[0]
		indLowLJ = np.where(PLJ == min(PLJ[IndMinRMSD]))[0][0]
		#indLowLJ = np.where(minRMSD == min(minRMSD))[0][0]
	else:
		indLowLJ = np.where(PLJ <= min(PLJ))[0][0]
	message = "Fragment chosen ({0}) has a energie of {1} kJ/mol \
 and RMSD of {2} A²".format(SeqGive[indLowLJ],\
	round(PLJ[indLowLJ],2)\
	, round(RMSD[indLowLJ],2))
	print message
	build_info(message+"\n")
	#Several minum value ofr PLJ we use the rmsd to
	#choose the best fragment
	return Coord[indLowLJ], OtherInfo[indLowLJ], RMSD[indLowLJ]


def BackupFromPDB(NamePDB, BegingPDB, EndPDB, ChainPDB, BoxDimension, path):
	Fichier = open("tmpFrag.pdb", "w")
	Fichier.write("MODEL 1\n")
	Fichier.write(str(BoxDimension))
	DicoRequest = {} #contains the list of Carbon alpha
	#and the rest of the principal chain (N, CO)
	splitted = []
	NbResid = 0
	FlagN = False;
	FlagC = False;
	FlagO = False;
	with open(path+NamePDB, "r") as filin:
		for line in filin:
			if line[:4] == 'ATOM' and int(line[22:26]) >= BegingPDB \
			and  int(line[22:26]) <= EndPDB and ChainPDB == line[21]:
				splitted = [
				line[:6], #'ATOM'
				line[6:11], #Atom's number
				line[12:16], #Atom's type (CA,NH,CZ...)
				line[17:20], #Amino acid's name
				line[21], #Chain's name
				line[22:26], #Residue's number
				line[30:38], #cartesian coordinates:
				line[38:46], #cartesian coordinates:
				line[46:54], #cartesian coordinates:
				line[77] #Atom's name
				]
				if NbResid != int(line[22:26]):
					NbResid = int(line[22:26])
					FlagN = False
					FlagC = False
					FlagO = False
				if splitted[2].split()[0] == 'N' and FlagN is False:
					if splitted[3] in AminoAcid.keys():
						AA = AminoAcid[splitted[3]]
					else:
						AA = 'X'
						print "The amino acid {0} is unknown !".format(int(splitted[5]))
					DicoRequest[int(line[22:26])] = CA(AA,np.array(
												[
												float(splitted[6]),
												float(splitted[7]),
												float(splitted[8])
												])
											)
					Fichier.write(line)
					FlagN = True
				if splitted[2].split()[0] == 'CA':
					DicoRequest[int(line[22:26])].PosiCA = np.array([
													float(splitted[6]),
													float(splitted[7]),
													float(splitted[8])
													])
					Fichier.write(line)
				if splitted[2].split()[0] == 'C' and FlagC is False:
					DicoRequest[int(line[22:26])].PosiC = np.array([
													float(splitted[6]),
													float(splitted[7]),
													float(splitted[8])
													])
					Fichier.write(line)
					FlagC = True
				if splitted[2].split()[0] == 'O' and FlagO is False:
					DicoRequest[int(line[22:26])].PosiO = np.array([
													float(splitted[6]),
													float(splitted[7]),
													float(splitted[8])
													])
					Fichier.write(line)
					FlagO = True
	Fichier.write("TER\n")
	Fichier.close()
	return DicoRequest

def ReplaceCoord(FinalFrag, PeptideCA, query,\
ClefFrag, SplitTemp, MissCA):
#PeptideCA beging to 0.
	for elment in MissCA.keys():
		beging = ""
		for i in xrange(len(SplitTemp[0:elment+1])):
			beging += SplitTemp[i]
		end = len(beging) + MissCA[elment][3] + 2
		beging = len(beging) - 2
		gap = [beging-1, end]
		print gap
		IndTemp = beging
			#ClefFrag[clef]
		for indice in xrange(ClefFrag[elment][0],ClefFrag[elment][1]+1):
			FinalFrag[elment][(' ',indice,' ')]
			PeptideCA[IndTemp] = CA(query[IndTemp],\
			list(FinalFrag[elment][(' ',indice,' ')]['N'].get_coord()))
			PeptideCA[IndTemp].PosiCA = list(FinalFrag[elment][(' ',indice,' ')]['CA'].get_coord())
			PeptideCA[IndTemp].PosiC = list(FinalFrag[elment][(' ',indice,' ')]['C'].get_coord())
			PeptideCA[IndTemp].PosiO = list(FinalFrag[elment][(' ',indice,' ')]['O'].get_coord())
			IndTemp += 1
	return PeptideCA

def WriteLigne(Fichier, ligne, cpt, elemt, AA, AtomType):
	X = round(ligne[0],2)
	Y = round(ligne[1],2)
	Z = round(ligne[2],2)
	Fichier.write("%-6s%5d "%('ATOM', cpt))
	Fichier.write("%4s%1s%3s "%(AtomType, " ",AtomType[0][:3]))
	Fichier.write("%1s%4d%1s   "%('A', elemt, " ") )
	Fichier.write("%8.3f%8.3f%8.3f%6.2f%6.2f"%(X,Y,Z, 1.00 ,42.00) )
	Fichier.write("          %2s%2s\n"%(str(cpt), 0))


def OutPutFinalModel(FinalModel, BoxDimension, PDB, CodeQuery):
	cpt = 1
	NomFichier = "NewModel"+CodeQuery
	Fichier = open(NomFichier, "w")
	Fichier.write("REMARK    GENERATED BY JAYSEN\n")
	Fichier.write("TITLE     structure based from pdb t=   "+str(1)+"\n")
	Fichier.write("REMARK    THIS IS A SIMULATION BOX\n")
	Fichier.write(str(BoxDimension))
	Fichier.write("MODEL "+str(1)+"\n")
	for elemt in xrange(len(FinalModel)):
		WriteLigne(Fichier, FinalModel[elemt].PosiN, cpt,\
	 	elemt+1,FinalModel[elemt].Name, "N")
		WriteLigne(Fichier, FinalModel[elemt].PosiCA, cpt,\
	 	elemt+1, FinalModel[elemt].Name, "CA")
		WriteLigne(Fichier, FinalModel[elemt].PosiC, cpt,\
	 	elemt+1, FinalModel[elemt].Name, "C")
		WriteLigne(Fichier, FinalModel[elemt].PosiO, cpt,\
	 	elemt+1, FinalModel[elemt].Name, "O")
		cpt += 1
	Fichier.write("TER\n")
	Fichier.close()


def GapTemplate(query, template, Fichier):
	TempSeq = ""
	QuerySeq = ""
	NameSeq = ""
	with open(Fichier, "r") as filin:
		for line in filin:
			if line.split(":")[0] == "structureX":
				NameSeq = line.split(":")[1]
				continue
			if NameSeq == template and line[0] != ">":
				TempSeq += line[:-1]
			if NameSeq == query and line[0] != ">":
				 QuerySeq += line[:-1]
	tmpQuery = ""
	tmpTempSeq = ""
	for i in xrange(min( len(QuerySeq), len(TempSeq)) ):
		if (QuerySeq[i] == "-" and TempSeq[i] == "-"):
			continue
		else:
			tmpQuery += QuerySeq[i]
			tmpTempSeq += TempSeq[i]
	return tmpQuery[:-1], tmpTempSeq[:-1]

def WhereGapBegings(sequence):
	flagbeg = False
	MissCA = {}
	clef = 0
	for indice in xrange(len(sequence)):
		if (sequence[indice] == "-" and flagbeg == False):
			beging = indice-1
			flagbeg = True
		if (sequence[indice] != "-" and flagbeg == True):
			end = indice
			flagbeg = False
			if(beging >= 2 and end <= len(sequence)-2):
				MissCA[clef] = [
								beging, end,
								sequence[beging-2:end+2+1],
								end-beging-1
								]
				clef += 1
	return MissCA

#################################################################
#						log info								#
#################################################################
def build_info(mess):
	filout = open('build_info', 'a')
	filout.write(mess+"\n")
	filout.close()


#################################################################
#								MAIN							#
#################################################################
Fichier = sys.argv[1] #"glob.pir"
CodeTemplate = sys.argv[2]#"1bin"
CodeQuery = sys.argv[3] #"1lh1"
PDB = sys.argv[4] #'1bin_chainA.pdb'
bash_commande("rm build_info")
parser = PDBParser()
structure = parser.get_structure('BioTempl', PDB)
TemplatePDB = structure

query, template = GapTemplate(CodeQuery, CodeTemplate, Fichier)
i = 0
RMSD = {}
while template[i] == "-":
 	template = template[i+1:]
	query = query[i+1:]

SplitTemp = template.split('-')
SplitTemp = filter(None, SplitTemp)
MissCA = WhereGapBegings(template)
SplitQuery = template.split('-')
SplitQuery = filter(None, SplitQuery)
GapQuery = WhereGapBegings(query)
TemplateCA, PeptideCA, BoxDimension = parse_template2(PDB, AminoAcid, MissCA)
print 'There are {0} gaps encountered in the template, that can be modelized'.format( len(MissCA.keys()) )
FinalFrag = {}
ClefFrag = {}
OtherInfo = {}
DicoCoord = {}
for clef in MissCA.keys():
	beg  = MissCA[clef][0]
	end  = MissCA[clef][1]
	print 'Looking for fragments for the gap: {0}.'.format(clef+1)
	SeqGive, Score, SizeFrag =  SmallGap(MissCA[clef],query,TemplateCA,\
								0.05,0.3,template)
	if len(SeqGive) == 0:
		del MissCA[clef]
		continue
	DicoCoord, OtherInfo = CoordMissFrag(SeqGive, BoxDimension,path)
	message = 'For the gap {0}, {1} fragments are found'.format(clef+1, len(Score))
	print message
	build_info(message+"\n")
	FinalFrag[clef],ClefFrag[clef],  RMSD[clef]= BestFragment(DicoCoord,\
									MissCA[clef],OtherInfo,TemplateCA,PeptideCA,\
									SplitTemp,clef,TemplatePDB)

PeptideCA = ReplaceCoord(FinalFrag, PeptideCA, query,\
ClefFrag, SplitTemp, MissCA)

OutPutFinalModel(PeptideCA, BoxDimension, PDB, CodeQuery)
print "The software software made a model in a {0:.2f} min".format((time.time() - start)/60.0)
