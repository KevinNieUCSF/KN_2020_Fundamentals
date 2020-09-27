#!/bin/env python3
#$ -S /bin/env python3            #-- the shell for the job
#$ -q gpu.q                #-- use the gpu queue
#$ -o sge_logs             #-- log output directory (MUST EXIST)
#$ -j y                    #-- tell the system that the STDERR and STDOUT should be joined
#$ -cwd                    #-- tell the job that it should start in your working directory
#$ -l mem_free=4G          #-- submits on nodes with enough free memory
#$ -l h_rt=2:00:00         #-- runtime limit 2 hours (opens up more gpu.q nodes)
#$ -R yes                  #-- SGE host reservation
#$ -N md_demo               #-- job name
#$ -binding linear:4       #-- bind the job to N adjacent cores - set to same val as GMX_CPUS

#import modules
import sys
import numpy as np
import scipy
import matplotlib as plt
import math
import pypdb as pd
"""
this code is made to recreate figure 1a and 1b from Laitaoja, et al. "Zinc Coordination Spheres in Protein Structures" 


a level of abstraction I am seeing here will be to have a parser that makes a list for each search term 
then we have to do some sorting of each of those resulting lists (eg. remove duplicates and realy similar structure)
and annotation maybe of the molecular weight for each of the structures. 
finally we plot them out 

"""

def parsepdb(inputquery1, inputquery2): #general pipeline, tentative
    #fetch files from pdb that contain zinc, were collected by NMR, x-ray  generate a method; this generates a list?
    query=pd.make_query(inputquery1, querytype='ExpTypeQuery')
    for hit in pd.do_search(query):
        print(hit)
    return None


def graph(pdblist):
    #generate a graph using the list generated from parsepdb


def pdbListFilter():
    #this function will take two search lists and merge for the final figure -- need to remove duplicate structures etc that "will impact further statistics"


def main():
    #currently trying to get the query function in parsepdb() to do two search parameters at once but the code only accepts one search term and one
    #query type, I think I will first look around to see any solutions and then email the coders
    inputquery1=input('Please input the experimental method term (use "NMR"):')
    inputquery2=input('Please input the advanced search term (use "zinc"):')
    pdblist=parsepdb(inputquery1, inputquery2)
    
    

    #graph(pdblist):



main()
#print(list(pd.get_info('1pvn'))
