#load in functions
library(stringr)
library(readxl)
library(data.table)
library(curl)
library(FELLA)
library(rstudioapi)

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

source("supp/figure2C_supp.R")

graph <- buildGraphFromKEGGREST(organism = "sce")
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph, databaseDir = tmpdir,
                   internalDir = FALSE,matrices = c("diffusion","hypergeom"),
                   normality = c("diffusion"),niter = 100)
fella.data <- loadKEGGdata(databaseDir = tmpdir,internalDir = FALSE,
                           loadMatrix = "diffusion")

mdat = read.csv("../results/phase/phase.csv")
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

plot(analysis.ds,method = "diffusion", data = fella.data, nlimit = nlimit, vertex.label.cex = vertex.label.cex, threshold = 0.05,vertex.size=30)
tab.all <- generateResultsTable(method = "diffusion",nlimit = nlimit,object = analysis.ds,data = fella.data,threshold = 0.05)
tab.all[tab.all$Entry.type == 'pathway',]


