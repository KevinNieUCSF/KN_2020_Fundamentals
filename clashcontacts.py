#! /usr/bin/env python3
import sys
import os 
import Bio.PDB
import numpy as np
import math
#kyle is big stinky poopy
try:
    os.chdir('/Users/kevin/KN_2020_Fundamentals/')
except:
    pass
def read_input():
    """ Read the PDB-formatted input file and return a structure """
    from Bio.PDB import PDBParser
    while True:
        try:
            input_file=str(input('PDB Filename:'))
            parser=PDBParser()
            structure=parser.get_structure("sample",input_file)
            break
        except:
            print('Invalid Filename, Please Try Again')
            continue
    return structure
def magnitude(distvec):
    return math.sqrt(sum(pow(element, 2) for element in distvec))

def calculate_contacts(structure):
    """ For every pair of residues in the structure, find if the residues
    or any of the alternate locations make contact.  If so add them to
    the list of contacts.
    [NOTE: not sure if this should be a list or a dictionary...] """
    clash=-2
    for atom in structure.get_atoms():
        for btom in structure.get_atoms():
            distvec=list(np.array(btom.get_coord()) - np.array(atom.get_coord()))
            print(magnitude(distvec))
            print('---')
            if magnitude(distvec)<=4:
                clash=clash+1
                print('Clash Found')
                print(clash)
            else:
                pass
    return clash
structure=read_input()
print("Clash Score: "+str(calculate_contacts(structure)))
print('All Done')
