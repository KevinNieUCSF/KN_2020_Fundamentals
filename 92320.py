#! /usr/bin/env python3
import sys
import os
import Bio.PDB
try:
    os.chdir('C:\\Users\\mailk\\py4e')
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

    return None

def calculate_contacts(structure):
    """ For every pair of residues in the structure, find if the residues
    or any of the alternate locations make contact.  If so add them to
    the list of contacts.
    [NOTE: not sure if this should be a list or a dictionary...] """
    return None
print('All Done')
