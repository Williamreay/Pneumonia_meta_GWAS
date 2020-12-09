#!/usr/bin/env python3

"""
Created on Mon Jan 06 2020

Pharmagenic enrichment score (PES) candidate pathway derivation - MAGMA v1.07 version

Author: William Reay - william.reay@uon.edu.au - github: https://github.com/Williamreay
"""

import os
import sys
import argparse
from subprocess import Popen
import pdb

##Add command line flags for automation of MAGMA protocol, set argsDict variable to be a dictionary containing the command line arguments


parser = argparse.ArgumentParser('MAGMA gene level and geneset level association analyses at different P value thresholds')
parser.add_argument('--gwasfile', metavar=['filename'], type = str, help = 'GWAS summary statistics in the following format - SNP, CHR, BP, P')
parser.add_argument('--upstream_boundary', metavar='[str]', type = str, help = "Upstream boundary for gene definition")
parser.add_argument('--downstream_boundary', metavar='[str]', type = str, help = "Downstream boundary for gene definition")
parser.add_argument('--geneloc', metavar=['filename'], type = str, help = 'Gene location file - ID, BP1, BP2')
parser.add_argument('--phenotype', metavar='[str]', type = str, help = 'GWAS phenotype name')
parser.add_argument('--bfile', metavar= '{prefix}', type = str, help = 'Reference population plink binary fileset - .bed, .bim + .fam')
parser.add_argument('--samplesize', metavar='[str]', type = str, help = 'Sample size which GWAS performed on - cases + controls')
parser.add_argument('--genesets', metavar=['filename'], type = str, help = 'Input pathways in standard MSigDB format, see MAGMA manual for further clarification')
inputFlags = parser.parse_args()

print(inputFlags)


##Annotate GWAS SNPs to genes

def geneAnnotation(upstream_boundary, downstream_boundary, gwasfile, geneloc, phenotype):
    print("Annotate SNPs from GWAS file to genes with specified boundaries")
    Popen("""./magma --annotate window=""" + upstream_boundary + """,""" + downstream_boundary + """ --snp-loc """ + gwasfile + """ --gene-loc """ + geneloc + """ --out """ + phenotype, shell=True).wait()

geneAnnotation(inputFlags.upstream_boundary, inputFlags.downstream_boundary, inputFlags.gwasfile, inputFlags.geneloc, inputFlags.phenotype)


##Generate different P value threshold files for input into the MAGMA model - all SNPs, P < 0.5, P < 0.05, P < 0.005, P < 0.00005

def pThresholds(gwasfile):
    print("Prunning SNP set at different P value thresholds for input into the MAGMA model")
    print("P < 0.5")
    Popen(""" awk ' $4 < 0.5 { print $1, $4} ' """ + gwasfile + """ > 0.5_threshold.loc""", shell=True).wait()
    print("P < 0.05")
    Popen(""" awk ' $4 < 0.05 { print $1, $4} ' """ + gwasfile + """ > 0.05_threshold.loc""", shell=True).wait()
    print("P < 0.005")
    Popen(""" awk ' $4 < 0.005 { print $1, $4} ' """ + gwasfile + """ > 0.005_threshold.loc""", shell=True).wait()


pThresholds(inputFlags.gwasfile)

##Gene level association analysis using SNPwise mean model at each P value threshold


def geneAssociation(bfile, gwasfile, samplesize, phenotype):
    print("MAGMA gene level association analysis at each p value threshold")
    print("All SNPs as input")
    Popen(""" ./magma --bfile """ + bfile + """ --pval """ + gwasfile + """ N=""" + samplesize + """ --gene-annot """ + phenotype + """.genes.annot --out """ + phenotype + """_allSNPs""", shell=True).wait()
    print("P < 0.5 SNPs as input")
    Popen(""" ./magma --bfile """ + bfile + """ --pval 0.5_threshold.loc N=""" + samplesize + """ --gene-annot """ + phenotype + """.genes.annot --out """ + phenotype + """_0.5_thresholdSNPs""", shell=True).wait()
    print("P < 0.05 SNPs as input")
    Popen(""" ./magma --bfile """ + bfile + """ --pval 0.05_threshold.loc N=""" + samplesize + """ --gene-annot """ + phenotype + """.genes.annot --out """ + phenotype + """_0.05_thresholdSNPs""", shell=True).wait()
    print("P < 0.005 SNPs as input")
    Popen(""" ./magma --bfile """ + bfile + """ --pval 0.005_threshold.loc N=""" + samplesize + """ --gene-annot """ + phenotype + """.genes.annot --out """ + phenotype + """_0.005_thresholdSNPs""", shell=True).wait()

geneAssociation(inputFlags.bfile, inputFlags.gwasfile, inputFlags.samplesize, inputFlags.phenotype)

##Geneset association on druggable pathways - at least one Tclin gene

def geneSetAssociation(genesets, phenotype):
    print("MAGMA geneset association analysis at each p value threshold")
    print("All SNPs as input")
    Popen(""" ./magma --gene-results """ + phenotype + """_allSNPs.genes.raw --set-annot """ + genesets + """ --out """ + phenotype + """allSNPs""", shell=True).wait()
    print("P < 0.5 SNPs as input")
    Popen(""" ./magma --gene-results """ + phenotype + """_0.5_thresholdSNPs.genes.raw  --set-annot """ + genesets + """ --out """ + phenotype + """_0.5_SNPs""", shell=True).wait()
    print("P < 0.05 SNPs as input")
    Popen(""" ./magma --gene-results """ + phenotype + """_0.05_thresholdSNPs.genes.raw --set-annot """ + genesets + """ --out """ + phenotype + """_0.05_SNPs""", shell=True).wait()
    print("P < 0.005 SNPs as input")
    Popen(""" ./magma --gene-results """ + phenotype + """_0.005_thresholdSNPs.genes.raw --set-annot """ + genesets + """ --out """ + phenotype + """_0.005_SNPs""", shell=True).wait()

geneSetAssociation(inputFlags.genesets, inputFlags.phenotype)
