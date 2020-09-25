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
"""this code is made to recreate figure 1a and 1b from Laitaoja, et al. "Zinc Coordination Spheres in Protein Structures" """

def parsepdb(input_param): #general pipeline, tentative
"""fetch files from pdb but filter out non zinc finger domain possessing ones, generate a method
to count how many of a certain molmass appears, this generates a list. We are iterating over the entire
database so I think this is most efficient"""
    return None
def graph(pdblist):
"""generate a graph using the list generated from parsepdb"""
def main():
    pdblist=parsepdb(input_param)
