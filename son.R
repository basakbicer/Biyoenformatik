# Gerekli paketlerin y??klenmesi
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor paketlerini y??klemek i??in
BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("hgu133plus2.db")
BiocManager::install("sva")

# Gerekli paketlerin y??klenmesi
library(GEOquery)
library(affy)
library(hgu133plus2.db)
library(sva)

# GSE23402, GSE26428 ve GSE28688 C'alD1EmalarD1nD1n verilerini C'ekme
GSE23402 <- getGEO("GSE23402")
GSE26428 <- getGEO("GSE26428")
GSE28688 <- getGEO("GSE28688")

# Verilerin uygun formatta alD1nmasD1
GSE23402 <- GSE23402[[1]][, 1:24]  # D0lk 24 C6rneDi al
GSE26428 <- GSE26428[[1]]
GSE28688 <- GSE28688[[1]]

# Verilerin log2 dC6nC<EC<mC<ne tabi tutulmasD1 ve 16 bit hassasiyete dC6nC<EtC<rC<lmesi
exprs(GSE23402) <- log2(exprs(GSE23402))
exprs(GSE26428) <- log2(exprs(GSE26428))
exprs(GSE28688) <- log2(exprs(GSE28688))

# TC<m veri setlerinin aynD1 hassasiyette olduDundan emin olunmasD1
exprs(GSE23402) <- round(exprs(GSE23402), 16)
exprs(GSE26428) <- round(exprs(GSE26428), 16)
exprs(GSE28688) <- round(exprs(GSE28688), 16)

# GSE23402, GSE26428 ve GSE28688 veri setlerinin deDiEken sayD1larD1nD1 kontrol etme
ncol(GSE23402)
ncol(GSE26428)
ncol(GSE28688)

# Veri setlerini birleEtirme
birlesmis_veri <- cbind(GSE23402, GSE26428, GSE28688)


# Veri iEleme ve makine C6Drenimi adD1mlarD1
set.seed(111)
sinir <- floor(0.80 * nrow(birlesmis_veri)) # Makine C6Drenmesi yC6ntemlerinde,
print(sinir) # EDitim iC'in ayrD1lan veri sayD1sD1nD1 gC6sterir.

ind <- sample.int(n = nrow(birlesmis_veri), size = sinir, replace = FALSE) # Crneklerimin indis sayD1sD1nD1 bulmak istiyorum
veriegitim <- birlesmis_veri[ind,]
veriegitim_df <- as.data.frame(veriegitim)
sinifegitim <- durum[ind]
veritest <- birlesmis_veri[-ind,]
siniftest <- durum[-ind]


veriegitim <- as.data.frame(veriegitim)

install.packages("rJava")
library(e1071)
library(RWeka)
library(rJava)
library(C50)
library(SMO)
modelj48 <- J48(sinifegitim ~ ., data = veriegitim) # Karar aDacD1 algoritmasD1
summary(modelj48)
plot(modelj48, type = "simple")
ongoru <- predict(modelj48, veritest)
con <- table(ongoru, siniftest)

dogruluk <- sum(diag(con)) / sum(con)
print(dogruluk)

# Cnemli genler
head(C5imp(modelj48), 10)
library(caret)
data_frame1 <- as.data.frame(sonveri)
data_frame1$durum <- durum
data_frame1$durum <- as.factor(data_frame1$durum)
# https://www.rdocumentation.org/packages/caret/versions/4.47/topics/train
model <- train(durum ~ ., data = data_frame1, method = "rf")
# Get variable importance
feature_importance <- varImp(model)
print(feature_importance)[10]

# Random Forest, C5.0 algoritmasD1 ile denemeler yapmanD1z.
# ArdD1ndan diDer yC6nteme gC6re GEO GDS4794 akciDer kanseri ve sizin bulduDunuz veriler C<zerine deneme yapD1nD1z.

library(randomForest)
modelrf <- randomForest(sinifegitim ~ ., data = veriegitim)
show(modelrf)

# Test iElemi ve modelin performansD1 iC'in
ongoru <- predict(modelrf, veritest)
con <- table(ongoru, siniftest)
dogruluk <- sum(diag(con)) / sum(con)
print(dogruluk)

library(pROC)
rocveri <- roc(as.numeric(siniftest), as.numeric(ongoru))
plot(rocveri, col = "blue")

auc(rocveri)

# Cnemli genleri ortaya C'D1kartmamD1z lazD1m
head(varImp(modelrf, scale = TRUE), 10)