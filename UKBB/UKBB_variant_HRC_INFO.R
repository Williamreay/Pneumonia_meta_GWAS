#############################################

## Selecting HRC sites and INFO > 0.8 to retain in the UKBB sample

## William Reay (2020)

##############################################

library(data.table)
library(dplyr)


setwd("~/data/UKBB/UKBB_variant_QC/")

## Read in HRC sites

HRC_raw <- fread("HRC.r1-1.GRCh37.wgs.mac5.sites.tab", header = T)

## Keep unique HRC variants with AF > 1e-4

HRC_keep <- HRC_raw %>% filter(AF_EXCLUDING_1000G > 1e-4) %>% select(ID)

HRC_keep <- unique(HRC_keep)

## Read in imputation metrics

INFO_raw <- fread()
