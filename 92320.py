#! /usr/bin/env python3

import sys
import os
import Bio.PDB
import numpy as np


def read_input(input_file):
    """ Read the PDB-formatted input file and return a structure """
    from Bio.PDB import PDBParser
    while True:
        try:
            parser=PDBParser()
            structure=parser.get_structure("sample",input_file)
            break
        except:
            print('Invalid Filename, Please Try Again')
            continue
    return structure

def calculate_contacts(structure):
    """ For every pair of residues in the structure, find if the residues
    or any of the alternate locations make contact.  If so add them to
    the list of contacts.
    [NOTE: not sure if this should be a list or a dictionary...] """
    clash=-2
    for atom in structure.get_atoms():
        if atom.get_name()=='CA' or atom.get_name()=='C' or atom.get_name=='N':
            continue
        else:
            for btom in structure.get_atoms():
                if btom.get_name()=='CA' or btom.get_name()=='C' or btom.get_name=='N':
                    continue
                else:
                    dist= atom - btom
                    if dist<=3:
                        clash+=1
                    else:
                        pass
    return clash
def main():
    try:
        os.chdir('/Users/kevin/KN_2020_Fundamentals/')
    except:
        pass
    input_file=input('Name of File, please:')
    structure=read_input(input_file)
    print("Clash Score: "+str(calculate_contacts(structure)))
main()
