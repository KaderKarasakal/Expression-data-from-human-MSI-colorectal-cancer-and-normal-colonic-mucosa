if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)
# get dataset
gse <- getGEO("GSE24514", GSEMatrix =TRUE)

pd <- pData(phenoData(gse[[1]]))

expression_matrix <- exprs(gse[[1]])

sample_names <- colnames(exprs(gse[[1]]))

tumor_samples <- sample_names[sample_names >= "GSM604484" & sample_names <= "GSM604517"]
normal_samples <- sample_names[sample_names >= "GSM604518" & sample_names <= "GSM604532"]
#separate tumor and normal samples
tumor_expression_matrix <- exprs(gse[[1]])[, tumor_samples]
normal_expression_matrix <- exprs(gse[[1]])[, normal_samples]

# created combined matrix with labels for PCA
combinedMatrix <- cbind(tumor_expression_matrix, normal_expression_matrix)
sampleLabels <- c(rep("T", ncol(tumor_expression_matrix)), rep("N", ncol(normal_expression_matrix)))

pcaResult <- prcomp(t(combinedMatrix), center = TRUE, scale. = TRUE)

colors <- ifelse(sampleLabels == "Tumor", "red", "blue")
# PCA plot tumor->red, normal->blue.
plot(pcaResult$x[,1], pcaResult$x[,2], col=colors, xlab="PC1", ylab="PC2", main="PCA of GSE24514")
legend("topright", legend=c("Tumor", "Normal"), col=c("red", "blue"), pch=1)

#applied t-test
pValues <- apply(combinedMatrix, 1, function(x) {
  t_test_result <- t.test(x ~ sampleLabels)
  return(t_test_result$p.value)
})

library(stats)
# for benjamini Hochgberg correction
adjusted_pValues <- p.adjust(pValues, method = "BH")

significant_BH <- rownames(combinedMatrix)[adjusted_pValues < 0.05]
print(significant_BH)
length(significant_BH) # 8183 significant gen varmış.

#distance matrix for hierarchial clustering
distMatrix <- dist(t(combinedMatrix), method = "euclidean")

hc <- hclust(distMatrix, method = "complete")

plot(hc, labels = sampleLabels, main = "Hierarchical Clustering ")

# created my trainData and testData
test_indices <- c(sample(which(sampleLabels == "T"), 1), sample(which(sampleLabels == "N"), 1))

train_indices <- setdiff(1:ncol(combinedMatrix), test_indices)

trainData <- t(combinedMatrix[, train_indices])
trainLabels <- sampleLabels[train_indices]
testData <- t(combinedMatrix[, test_indices])
testLabels <- sampleLabels[test_indices]

# for SVM model
install.packages("e1071")

library(e1071)

svmModel <- svm(trainData, as.factor(trainLabels), kernel = "linear")

predictions <- predict(svmModel, testData)

result <- data.frame(Actual = testLabels, Predicted = predictions)

print(result)
#            Actual Predicted
#GSM604507      T         T
#GSM604526      N         N

# for apply Gaussian Graphical Model (GGM)
topGenes <- head(significant_BH, 5)

topGenesData <- combinedMatrix[topGenes, ]


install.packages("glasso")

library(glasso)

# Estimation of the sparse inverse covariance matrix
sparseCovariance <- glasso(cov(topGenesData), rho = 0.01)  

pearsonCorrelation <- cor(topGenesData)

print(sparseCovariance$wi)

print(pearsonCorrelation)

