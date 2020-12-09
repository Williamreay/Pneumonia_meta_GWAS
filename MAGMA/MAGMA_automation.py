#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 15:15:23 2019
Automation of MAGMA pipeline
@author: williamreay
"""

import os
import sys
import argparse
from subprocess import Popen

##Add command line flags for automation of MAGMA protocol, set argsDict variable to be a dictionary containing the command line arguments


parser = argparse.ArgumentParser('MAGMA gene level and geneset level association analyses')
parser.add_argument('--gwasfile', metavar=['filename'], type = str, help = 'GWAS summary statistics in the following format - SNP, CHR, BP, P')
parser.add_argument('--geneloc', metavar=['filename'], type = str, help = 'Gene location file - ID, BP1, BP2')
parser.add_argument('--phenotype', metavar='[str]', type = str, help = 'GWAS phenotype name')
parser.add_argument('--bfile', metavar= '{prefix}', type = str, help = 'Reference population plink binary fileset - .bed, .bim + .fam')
parser.add_argument('--samplesize', metavar='[str]', type = str, help = 'Sample size which GWAS performed on - cases + controls')
parser.add_argument('--genesets', metavar=['filename'], type = str, help = 'Input pathways in standard MSigDB format, see MAGMA manual for further clarification')
inputFlags = parser.parse_args()

print(inputFlags)


##Annotate GWAS SNPs to genes, set genic boundaries to encompass 5 kb upstream and 1.5 kb downstream

def geneAnnotation(gwasfile, geneloc, phenotype):
    print("Annotate SNPs from GWAS file to genes with 5kb upstream and 1.5 kb downstream boundaries")
    Popen("""./magma --annotate window=35,10 --snp-loc """ + gwasfile + """ --gene-loc """ + geneloc + """ --out """ + phenotype, shell=True).wait()

geneAnnotation(inputFlags.gwasfile, inputFlags.geneloc, inputFlags.phenotype)



##Gene level association analysis using SNPwise mean model at each P value threshold


def geneAssociation(bfile, gwasfile, samplesize, phenotype):
    print("MAGMA gene level association analysis at each p value threshold")
    Popen(""" ./magma --bfile """ + bfile + """ --pval """ + gwasfile + """ N=""" + samplesize + """ --gene-annot """ + phenotype + """.genes.annot --out """ + phenotype + """_commonSNPs""", shell=True).wait()

geneAssociation(inputFlags.bfile, inputFlags.gwasfile, inputFlags.samplesize, inputFlags.phenotype)

##Geneset association on MSigDB pathways

def geneSetAssociation(genesets, phenotype):
    print("MAGMA geneset association analysis at each p value threshold")
    Popen(""" ./magma --gene-results """ + phenotype + """_commonSNPs.genes.raw --set-annot """ + genesets + """ --out """ + phenotype + """commonSNPs""", shell=True).wait()

geneSetAssociation(inputFlags.genesets, inputFlags.phenotype)
