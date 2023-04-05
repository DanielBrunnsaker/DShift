
library(stringr)
library(readxl)
library(data.table)
library(curl)
library(FELLA)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/supp.R")

graph <- buildGraphFromKEGGREST(organism = "sce")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph, databaseDir = tmpdir,internalDir = FALSE,matrices = c("diffusion","hypergeom"),normality = c("diffusion"),niter = 100)
fella.data <- loadKEGGdata(databaseDir = tmpdir,internalDir = FALSE,loadMatrix = "diffusion")

#### strains in ethanol phase, select the one if interest

#mdat = read.csv('../results/strain/ethanol/TDA1.csv')
#mdat = read.csv('../results/strain/ethanol/YGR067C.csv')
#mdat = read.csv('../results/strain/ethanol/RME1.csv')
mdat = read.csv('../results/strain/ethanol/DLD3.csv')
#mdat = read.csv('../results/strain/ethanol/GAL11.csv')
#mdat = read.csv('../results/strain/ethanol/OCA1.csv')
#mdat = read.csv('../results/strain/ethanol/FAA1.csv')
#mdat = read.csv('../results/strain/ethanol/RTS3.csv')
#mdat = read.csv('../results/strain/ethanol/PCL1.csv')
#mdat = read.csv('../results/strain/ethanol/MEK1.csv')

#### strains in gflucose phase

#mdat = read.csv('../results/strain/glucose/TDA1.csv')
#mdat = read.csv('../results/strain/glucose/YGR067C.csv')
#mdat = read.csv('../results/strain/glucose/RME1.csv')
#mdat = read.csv('../results/strain/glucose/DLD3.csv')
#mdat = read.csv('../results/strain/glucose/GAL11.csv')
#mdat = read.csv('../results/strain/glucose/OCA1.csv')
#mdat = read.csv('../results/strain/glucose/FAA1.csv')
#mdat = read.csv('../results/strain/glucose/RTS3.csv')
#mdat = read.csv('../results/strain/glucose/PCL1.csv')
#mdat = read.csv('../results/strain/glucose/MEK1.csv')

####
#Remove duplicates (keep the ones with the highest certainty of identification, e.g. MS2 spectra)
to_remove = c('243.05', '132.10147','156.07756','160.0429',
              '160.04425','160.04459','136.06218','136.06087',
              '296.06558','376.1325','396.1585','250.17918')
mdat = remove_unwanted_features(mdat, to_remove)
mdat$inchikey = NA
mdat$KEGG = NA

annomdat = read_excel('../data/metabolomics/B1_filtered.xlsx', sheet = 'Sheet1')[c(4,5,14)]
colnames(annomdat) = c('mz','name','inchi')
kegg_inchi = read.csv('../other/kegg_map.csv')

#for each metabolite
mdat = add_inchi_kegg(mdat, annomdat, kegg_inchi)
threshold = 0.1
tempdat = p_threshold (mdat, threshold)
compounds.ds = tempdat$KEGG

analysis.ds <- defineCompounds(compounds = compounds.ds,data = fella.data)





analysis.ds <- runDiffusion(object = analysis.ds,data = fella.data,approx = "normality")

nlimit <- 300
vertex.label.cex <- 1
plot(analysis.ds,method = "diffusion", data = fella.data, nlimit = nlimit, vertex.label.cex = vertex.label.cex, threshold = 0.01,vertex.size=30)
tab.all <- generateResultsTable(method = "diffusion",nlimit = nlimit,object = analysis.ds,data = fella.data,threshold = 0.01)
tab.all[tab.all$Entry.type == 'pathway',]



table = mdat[c('KEGG','logFC','AveExpr','t','P.Value')]
table = table[table$KEGG != 'No result',]
table = na.omit(table)
library(xlsx)
write.xlsx(table, "/Volumes/Samsung_T5/DShift_repo/DShift/tables/RME1_E.xlsx", row.names = FALSE)


