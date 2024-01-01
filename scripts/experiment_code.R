 setwd("D:/master/TFM_GEN/scripts")
# options(timeout = max(300, getOption("timeout")))
## CARGAR SERIES Y DATOS DE PLATAFORMA DESDE GEO
library(GEOquery)
gset <- getGEO("GSE45255", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)
varLabels(gset)
## SELECCION DE VARIABLES FENOTIPICAS DE INTERES
histological_grade <- gset$`histological grade:ch1`
er_status <- gset$`er status:ch1`
her2_status <- gset$`her2 status:ch1`
ln_status <- gset$`ln status:ch1`
size <- gset$`size (mm):ch1`

#GENERAR DATAFRAME DE VARIABLES PARA LUEGO ANEXAR
fenot <- data.frame(histological_grade, er_status, her2_status, ln_status, size)
rownames(fenot) = colnames(ex)

# CARGA DE LIBRERIAS Y DESCARGA DE DATOS

gcel = getGEOSuppFiles("GSE45255")

# DESCOMPRIMIR DATOS Y LECTURA
setwd("GSE45255/")
system("tar xvf GSE45255_RAW.tar")
library(affy)
gse45255raw = ReadAffy()
#GUARDANDO DATOS DE EXPRESION
#save(gse45255raw,file = "gse45255raw.rda")
#CARGANDO DATOS Y PREPROCESANDO MEDIANTE METODO MAS5
# Metodo que realiza una correción de fondo y el cálculo del valor de expresión.
load("gse45255raw.rda")
gse45255_mas5 = affy::mas5(gse45255raw)
pData(gse45255_mas5) <- fenot
# GUARDANDO DATOS PREPROCESADOS
#save(gse45255_mas5,file = "gse45255_mas5.rda")
load("gse45255_mas5.rda")

#REMOVIENDO NAs presentes
filter <- which(gse45255_mas5@phenoData@data$histological_grade == 'NA' | gse45255_mas5@phenoData@data$her2_status == "NA")
gse45255_mas5 <- gse45255_mas5[,-c(filter)]

# FILTRANDO GENES
#Se utiliza el rango intercuartilico y la desviación típica. Para remover genes con poca relevancia
library(genefilter)
gse45255.filt1 = nsFilter(gse45255_mas5,var.func=IQR,var.cutoff=0.5, require.GOBP=TRUE)
gse45255.filt2 = nsFilter(gse45255_mas5,var.func=sd,var.cutoff=0.5, require.GOBP=TRUE)
sel = intersect(featureNames(gse45255.filt1), featureNames(gse45255.filt2))
gse45255_mas5_filt = gse45255_mas5[sel,]


# ANALISIS DE COMPONENTES PRINCIPALES
pca = prcomp(t(exprs(gse45255_mas5_filt)),scale=TRUE,center=TRUE)
summary(pca)
# Representación gráfica de las dos primeras componentes para las variables ln_status, er_status y her2_status
pacman::p_load(ggfortify)
df0 = t(exprs(gse45255_mas5_filt))
ln_status = pData(gse45255_mas5_filt)[, "ln_status"]
df = data.frame(ln_status,df0)
my_plot <- autoplot(pca, data = df, colour="ln_status")
plot(my_plot)

er_status = pData(gse45255_mas5_filt)[, "er_status"]
df = data.frame(er_status,df0)
my_plot <- autoplot(pca, data = df, color="er_status")
plot(my_plot)

her2_status = pData(gse45255_mas5_filt)[, "her2_status"]
df = data.frame(her2_status,df0)
my_plot <- autoplot(pca, data = df, color="her2_status")
plot(my_plot)

histological_grade = pData(gse45255_mas5_filt)[, "histological_grade"]
df = data.frame(histological_grade,df0)
my_plot <- autoplot(pca, data = df, color="histological_grade")
plot(my_plot)

#MAPA DE CALOR
library(pheatmap)
corMatrix <- cor(exprs(gse45255_mas5_filt),use="c")
pheatmap(corMatrix)

## Test de Fisher
# Division de genes en grupos significativos en base a la información de Gene Ontology para luego aplicar un test de Fisher.

pacman::p_load(genefilter,multtest)

y=pData(gse45255_mas5_filt)[,"ln_status"]
# 1 Determinar grupos de genes significativos:
tt = rowttests(gse45255_mas5_filt,factor(y))
p0 = tt$p.value
p1 = mt.rawp2adjp(p0, "BH")
orden.original = order(p1$index)
p.BH = p1$adjp[orden.original,2]
significativos = which(p.BH < 0.05)
pacman::p_load(hgu133a.db)
# 2 Construccion de universo de genes verificando que no exista duplicidades
G1.entreizd = unlist(mget(featureNames(gse45255_mas5_filt), hgu133aENTREZID))
anyDuplicated(G1.entreizd)
# 3 Cambiando la identificación para que relacionar con la identificación de la base de datos de los genes identificativos:
seleccionados = unlist(mget(featureNames(gse45255_mas5_filt[significativos,]),hgu133aENTREZID))
# 4 Aplicación de Test de Fisher (SE PLICAN ESTOS 4 PASOS CON EL RESTO DE VARIABLES)
pacman::p_load(GO.db,Category,GOstats)
params = new("GOHyperGParams", geneIds = seleccionados,universeGeneIds = G1.entreizd,annotation = annotation(gse45255_mas5_filt), ontology = "BP",pvalueCutoff = 0.001,conditional = FALSE,testDirection = "over")
overRepresented = hyperGTest(params)
htmlReport(overRepresented, file = "GSE45255_ln_status_overRepresented.html")
head(summary(overRepresented))

y=pData(gse45255_mas5_filt)[,"er_status"]
tt = rowttests(gse45255_mas5_filt,factor(y))
p0 = tt$p.value
p1 = mt.rawp2adjp(p0, "BH")
orden.original = order(p1$index)
p.BH = p1$adjp[orden.original,2]
significativos = which(p.BH < 0.05)
G1.entreizd = unlist(mget(featureNames(gse45255_mas5_filt), hgu133aENTREZID))
anyDuplicated(G1.entreizd)
seleccionados = unlist(mget(featureNames(gse45255_mas5_filt[significativos,]),hgu133aENTREZID))
params = new("GOHyperGParams", geneIds = seleccionados,universeGeneIds = G1.entreizd,annotation = annotation(gse45255_mas5_filt), ontology = "BP",pvalueCutoff = 0.001,conditional = FALSE,testDirection = "over")
overRepresented = hyperGTest(params)
htmlReport(overRepresented, file = "GSE45255_er_status_overRepresented.html")
head(summary(overRepresented))

y=pData(gse45255_mas5_filt)[,"her2_status"] 
tt = rowttests(gse45255_mas5_filt,factor(y))
p0 = tt$p.value
p1 = mt.rawp2adjp(p0, "BH")
orden.original = order(p1$index)
p.BH = p1$adjp[orden.original,2]
significativos = which(p.BH < 0.05)
G1.entreizd = unlist(mget(featureNames(gse45255_mas5_filt), hgu133aENTREZID))
anyDuplicated(G1.entreizd)
seleccionados = unlist(mget(featureNames(gse45255_mas5_filt[significativos,]),hgu133aENTREZID))
params = new("GOHyperGParams", geneIds = seleccionados,universeGeneIds = G1.entreizd,annotation = annotation(gse45255_mas5_filt), ontology = "BP",pvalueCutoff = 0.001,conditional = FALSE,testDirection = "over")
overRepresented = hyperGTest(params)
htmlReport(overRepresented, file = "GSE45255_her2_status_overRepresented.html")
head(summary(overRepresented))

#GSA
# Análisis de cojutnos de genes para determinar cuáles tienen asociación positiva con la variable y comprobar si el resultado es similar al obtenido en el análisis anterior. 
library(annotate)
annotation(gse45255_mas5)
library(hgu133a.db)
library(GSEABase)
# Dividir los datos de los genes en conjuntos de acuerdo a los grupos definidos en Gene Ontology con la anotación correspondiente de los genes. 
gse45255.gsc = GeneSetCollection(gse45255_mas5_filt,setType=GOCollection())
names(gse45255.gsc) = unlist(lapply(gse45255.gsc,setName))
gsc = gse45255.gsc

gruposGrandes = gsc[which(sapply(geneIds(gsc),length) > 50)]
gse = gse45255_mas5_filt
# se realiza el procedimiento con las variables "ln_status", "er_status" y "her2_status".
node.num = as.numeric(factor(pData(gse)[,"ln_status"]))
pacman::p_load(GSA)
gse45255.gsa = GSA(exprs(gse),node.num, genenames=featureNames(gse),genesets=geneIds(gsc),resp.type="Two class unpaired", nperms=1000)
GSA.plot(gse45255.gsa)
#negativa
(ind.lo = which(gse45255.gsa$pvalues.lo <.05))
head(names(gsc[ind.lo]))
#positiva
(ind.hi = which(gse45255.gsa$pvalues.hi <.05))
head(names(gsc[ind.hi]))
#guardado en archivos
neg=names(gsc[ind.lo])
pos=names(gsc[ind.hi])
write.csv(neg, file = "neg_ln_status.csv")
write.csv(pos, file = "pos_ln_status.csv")
#er_status
node.num = as.numeric(factor(pData(gse)[,"er_status"]))
pacman::p_load(GSA)
gse45255.gsa = GSA(exprs(gse),node.num, genenames=featureNames(gse),genesets=geneIds(gsc),resp.type="Two class unpaired", nperms=1000)
GSA.plot(gse45255.gsa)

(ind.lo = which(gse45255.gsa$pvalues.lo <.05))
head(names(gsc[ind.lo]))
(ind.hi = which(gse45255.gsa$pvalues.hi <.05))
head(names(gsc[ind.hi]))
neg=names(gsc[ind.lo])
pos=names(gsc[ind.hi])
write.csv(neg, file = "neg_er_status.csv")
write.csv(pos, file = "pos_er_status.csv")
#her2_status
node.num = as.numeric(factor(pData(gse)[,"her2_status"]))
pacman::p_load(GSA)
gse45255.gsa = GSA(exprs(gse),node.num, genenames=featureNames(gse),genesets=geneIds(gsc),resp.type="Two class unpaired", nperms=1000)
GSA.plot(gse45255.gsa)

(ind.lo = which(gse45255.gsa$pvalues.lo <.05))
head(names(gsc[ind.lo]))
(ind.hi = which(gse45255.gsa$pvalues.hi <.05))
head(names(gsc[ind.hi]))
neg=names(gsc[ind.lo])
pos=names(gsc[ind.hi])
write.csv(neg, file = "neg_her2_status.csv")
write.csv(pos, file = "pos_her2_status.csv")

