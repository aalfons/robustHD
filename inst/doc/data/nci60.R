# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------


## This is a modified version of the script 'process_data.R' in the following
## GitHub repository: https://github.com/aalfons/nci60
##
## The aforementioned script was used to preprocess the data for the article
## Alfons, Croux & Gelper (2013, The Annals of Applied Statistics), where
## information on the cancer cell lines, proteins, and genes was not relevant
## and therefore discarded.  This script adds this information since it may
## be relevant to users.
##
## The raw data were downloaded as zip-files from the following website:
## https://discover.nci.nih.gov/cellminer/


## load and process protein expression data
# load data
fileProtein <- "https://raw.githubusercontent.com/aalfons/nci60/main/nci60_Protein__Lysate_Array_log2.txt/nci60_Protein__Lysate_Array_log2.txt"
protein <- read.table(fileProtein, sep = "\t", skip = 1)
# extract meta information
proteinInfo <- protein[, 1:4]
names(proteinInfo) <- c("Experiment", "Probe", "Symbol", "ID")
# remove those descriptive columns, then convert to matrix and transpose
# because the cancer cell lines are stored in columns instead of rows
protein <- t(as.matrix(protein[, -(1:4)]))
rownames(protein) <- 1:nrow(protein)


## load and prepare gene expression data
# load data
fileGene <- "https://raw.githubusercontent.com/aalfons/nci60/main/nci60_RNA__Affy_HG_U133(A_B)_GCRMA.txt/nci60_RNA__Affy_HG_U133(A_B)_GCRMA.txt"
gene <- read.table(fileGene, sep = "\t", na.strings = "-", skip = 1)
# only keep rows that correspond to the Affymetrix HG-U133A chip
keep <- gene[, 1] == "RNA: Affy HG-U133(A) GCRMA"
gene <- gene[keep, ]
# extract meta information
geneInfo <- gene[, 1:4]
names(geneInfo) <- c("Experiment", "Probe", "Symbol", "ID")
# remove those descriptive columns, then convert to matrix and transpose
# because the cancer cell lines are stored in columns instead of rows
gene <- t(as.matrix(gene[, -(1:4)]))
rownames(gene) <- 1:nrow(gene)


## load and prepare metadata for cancer cell lines
# load data
fileCellLine <- "https://github.com/aalfons/nci60/raw/main/nci60_RNA__Affy_HG_U133(A_B)_GCRMA.txt/nci60_cellline_metadata.xls"
tempFile <- tempfile(fileext = ".xls")
download.file(fileCellLine, destfile = tempFile, mode = "wb", quiet = TRUE)
cellLine <- readxl::read_excel(tempFile, na = "?", skip = 1)
# extract information on the cell lines
cellLineInfo <- as.data.frame(cellLine[, 1:15])
# change some column names
names(cellLineInfo)[1] <- "Name"
names(cellLineInfo)[2] <- "Tissue"
names(cellLineInfo)[3] <- "Age"
names(cellLineInfo)[4] <- "Sex"
names(cellLineInfo)[5] <- "PriorTreatment"
names(cellLineInfo)[7] <- "Histology"
names(cellLineInfo)[8] <- "Source"
names(cellLineInfo)[9] <- "Ploidy"
names(cellLineInfo)[11] <- "MDR"
names(cellLineInfo)[12] <- "DoublingTime"


## gene expression data contains one observation with only NAs, so remove it
remove <- which(apply(is.na(gene), 1, all))
protein <- protein[-remove, ]
gene <- gene[-remove, ]
cellLineInfo <- cellLineInfo[-remove, ]

## save data to RData file
save(protein, gene, proteinInfo, geneInfo, cellLineInfo,
     file = "data/nci60.RData")
