#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

gVCF_file <- args[1]
RNA_VCF_file <- args[2]
chromosome <- args[3]
start_pos <- args[4]
end_pos <- args[5]

write.table(args[1],"testing")
