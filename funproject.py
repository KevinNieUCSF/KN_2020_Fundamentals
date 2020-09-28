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
"""this code is made to recreate figure 1a and 1b from Laitaoja, et al. "Zinc Coordination Spheres in Protein Structures" """

def parsepdb(inputquery1, inputquery2,inputquery3): #general pipeline, tentative
    """fetch IDs from pdb that contain zinc, then compares the list to IDs which used NMR or xray, returns two
    lists that contain zinc containing structures that used NMR or xray"""
    querymain=pd.make_query(inputquery1, querytype='AdvancedKeywordQuery')
    queryexp1=pd.make_query(inputquery2, querytype='ExpTypeQuery')
    queryexp2=pd.make_query(inputquery3, querytype='ExpTypeQuery')
    h=0
    for hit in pd.do_search(querymain):
        h+=1
    print(h)
    return None
def graph(pdblist):
    """generate a graph using the list generated from parsepdb"""
def main():
    #currently trying to get the query function in parsepdb() to do two search parameters at once but the code only accepts one search term and one
    #query type, I think I will first look around to see any solutions and then email the coders
    inputquery1=input('Please input the advanced search term (use "zinc"):')
    inputquery2=input('Please input experimental method 1 (use "NMR"):')
    inputquery3=input('Please input experimental method 1 (use "X-RAY"):')
    pdblist=parsepdb(inputquery1, inputquery2,inputquery3)
    #graph(pdblist):
main()
