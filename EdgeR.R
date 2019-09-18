rm(list = ls())
pacman::p_load("edgeR", "SummarizedExperiment")

setwd("/home/antonio/Escritorio/Antonio/Bioinformatica/Cordoba/Oscar/paraAntonio/pulmon_adenocarcinoma")
data <- read.csv("luad_gene_read_paired-normal-tumor.txt", sep = "\t")
head(data)

a <- colnames(data)
condicion <- factor(substr(a, nchar(a), nchar(a)))

dge <- DGEList(data, group = condicion)
keep <- rowSums(cpm(dge)>1) >= 57 # Filtro para aquellos genes que tengan al menos cpm >= 1 en la mitad de las muestras
dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts) # Recalculo el tamaño de libreria
dge.n <- calcNormFactors(dge, method = "TMM")
dge.c <- estimateCommonDisp(dge.n)
dge.t <- estimateTagwiseDisp(dge.c)
et.c <- exactTest(dge.c)
et.t <- exactTest(dge.t)

topTags(et.c)
topTags(et.t)

# Determinamos los genes significativos

de.c <- decideTestsDGE(et.c, adjust.method = "BH", p.value = .05)
de.t <- decideTestsDGE(et.t, adjust.method = "BH", p.value = .05)

summary(de.t) # Resumen de los genes
summary(de.c)

# Genes up y down para cada clase
# Dispersión común

# Upregulados

upre.c <- which(de.c@.Data == 1)
upre.c <- et.c[upre.c,]
upre.c <- rownames(topTags(upre.c, n = 100))

# Downregulados

downre.c = which(de.c@.Data == -1)
downre.c = et.c[downre.c,]
downre.c = rownames(topTags(downre.c, n = 100))

# Cada gen tiene su dispersión

# Upregulados

upre.t = which(de.t@.Data == 1)
upre.t = et.t[upre.t,]
pvup = upre.t$table$PValue
upre.t = rownames(topTags(upre.t, n = "all")) 

# Downregulados

downre.t = which(de.t@.Data == -1)
downre.t = et.t[downre.t,]
pvdown = downre.t$table$PValue
downre.t = rownames(topTags(downre.t, n = "all")) 

# Hasta aquí hemos obtenido los nombres de los distintos genes up y down
# por los métodos de dispersión común y entre genes.
# Agrupamos los genes (Grupos de genes con expresión diferencial significativa)

forenrich <- data.frame(genes = c(upre.t, downre.t), pval = c(pvup, pvdown))
ordenado <- forenrich[order(-forenrich$pval),]
ordenado$genes <- gsub(".*\\|", "", ordenado$genes)
final <- ordenado$pval
names(final) <- ordenado$genes

save(final, file = "significativosenrichmentLUADpvalores.rda")

genes.c = c(upre.c, downre.c) # Dispersión común
genes.t = c(upre.t, downre.t) # Dispersión por gen, vamos a usar estos genes para Emerging Patterns.

save(genes.t, file = "GenessigLUAD.rda") # Guardamos los genes significativos 


# ### Análisis de grupos mediante Fisher (Category+GOStats)
# # Construcción del universo de genes
# 
# ## Aqui tengo que quedarme con los entrezid, nose si vale que me quede yo antes con ellos
# ## haciendole un split al nombre
# 
# nombres <- row.names(data)
# head(nombres)
# a <- sapply(strsplit(nombres,"\\|"), tail, 1)
# table(a)
# remove <- c("?")
# hgnc <- a[! a %in% remove]
# table(hgnc)
# a <- sapply(strsplit(genes.t, "\\|"), tail, 1)
# 
# pacman::p_load(GOstats, org.Hs.eg.db) # Cargamos más librerías necesarias
# 
# G2.entrezid <- unlist(hgnc)
# G2.entrezid
# G1.entrezid <- unique(G2.entrezid)
# G1.entrezid
# anyDuplicated(G1.entrezid) # 0
# 
# params.t <- new("GOHyperGParams", geneIds = a, universeGeneIds = G1.entrezid, ontology = "BP", pvalueCutoff = .05, conditional = F, 
#                 testDirection = "over", annotation = "hgu95av2.db")
# overRepresented.t = hyperGTest(params.t)
# summary(overRepresented.t)
# 
# 
# ## Para mirar la densidad de los Top 100 genes Up y Down
# conteos <- data.frame(dge.n$counts)
# head(conteos)
# conteos <- conteos[genes.t,]
# head(conteos)
# save(conteos, file = "Top100SigedgeR.r")
