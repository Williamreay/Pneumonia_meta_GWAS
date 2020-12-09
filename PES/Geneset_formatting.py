#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 11:48:31 2019

#Automation of PES score generation

@author: williamreay
"""

import os
import sys
import argparse
from subprocess import Popen


##Add command line flags

parser = argparse.ArgumentParser('Generation of PES scores for candidate pathways')
parser.add_argument('--geneset_ID', metavar='[str]', type = str, help = 'Name of geneset [case sensitive]'),
parser.add_argument('--path_to_gmt_file',metavar='[str]',type= str, help="Path to .gmt file")
inputFlags = parser.parse_args()

print(inputFlags)

#Make directory to store output files

Popen(""" mkdir """ + inputFlags.geneset_ID, shell=True).wait()

#Extract genes within the candidate set
Popen(""" grep -w """ + inputFlags.geneset_ID + """ """ + inputFlags.path_to_gmt_file + """MsigDB_canonical_and_hallmark_Tclin_v_6.1.0.gmt | tr "\t" "\n" | awk 'NR > 2 { print }' > """ \
      + inputFlags.geneset_ID + """_genes.txt """, shell=True).wait()


#Move to directory
Popen(""" mv """ + inputFlags.geneset_ID + """_genes.txt """ + inputFlags.geneset_ID, shell=True).wait()
