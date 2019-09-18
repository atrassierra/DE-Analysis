rm(list = ls())
library(EBSeq)
setwd("/home/antonio/Escritorio/isoformas/")

# Leemos los datos
data <- read.csv("brca_iso_read_paired.txt", sep = "\t")

# Asignamos los grupos experimentales, en este caso Normales vs Tumorales
a <- colnames(data)
condicion <- as.factor(substr(a, nchar(a), nchar(a)))

# EBSeq necesita por un lado de la matriz de conteos, que puede estar en formato RSEM,
# además de dos vectores con los nombres de las isoformas y los nombres de los genes a los que pertenecen las isoformas.
# Hay que tener en cuenta que en el vector de nombres de genes, estos genes se van a repetir si tienen más de una isoforma.
# Es decir, gen 1, gen 2, gen 2, gen 3 se corresponde con iso 1-1, iso 2-1, iso 2-2, iso 3-1, etc. Es importante que ocupen
# las posiciones que corresponden en el vector.

# La matriz de conteos la tenemos directamente en data
# Para el nombre de genes y de isoformas:

genesiso <- rownames(data)
# quitar <- grep("\\?", genesiso) # Quitar genes hipotéticos 
genes <- unlist(sapply(strsplit(genesiso, ","), "[", 1))
isoformas <- unlist(sapply(strsplit(genesiso, ","), "[", 2))

data.matriz <- data.matrix(data)
rownames(data.matriz) <- isoformas
data.size <- MedianNorm(data.matriz) # Isoform-level library size factors

# En el siguiente paso tenemos que elegir en cuantos grupos queremos que nos separe los genes dependiendo del número
# de isoformas que estos tengan. Por ejemplo, si nuestros genes tienen entre 1 y 3 isoformas pondremos 3 grupos.
# Para saber cuántos tenemos, podemos ver cuantos genes tenemos duplicados n veces. O el gen con más duplicados.

duplicados <- table(genes)
NgList <- GetNg(isoformas, genes, TrunThre = 3)
IsoNgTrun <- NgList$IsoformNgTrun

# Análisis de expresión diferencial

IsoEBOut <- EBTest(Data = data.matriz, NgVector = IsoNgTrun, Conditions = condicion, sizeFactors = data.size, maxround = 5)

# Isoformas diferencialmente expresadas

deiso <- GetDEResults(IsoEBOut, FDR = .05)              
deiso
IsoFC <- PostFC(IsoEBOut)
IsoFC$PostFC