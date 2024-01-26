#Question 1: Downloaded the GSE data
#GSE24514

library(GEOquery)
data_set <- getGEO("GSE24514", AnnotGPL = TRUE)
data_set <- data_set[[1]]
m_data <- exprs(data_set) #49 sample: first 34 sample tumor, last 15 sample normal
log_data <- log2(m_data) #data mızı normal dağılıma daha çok yaklaştırmak için bu aşamayı gerçekleştirdik
hist(log_data)

#create the transpose data, samples in rows
tr_data <- t(log_data)

#Question 2: Performed the Principal Component Analysis (PCA)

install.packages("ggplot2")  #for created the plots
library(ggplot2)

sample_names <- colnames(m_data)
tumor_samples <- sample_names[sample_names >= "GSM604484" & sample_names <= "GSM604517"]
normal_samples <- sample_names[sample_names >= "GSM604518" & sample_names <= "GSM604532"]
#separate tumor and normal samples
tumor_expression_matrix <- m_data[, tumor_samples]
normal_expression_matrix <- m_data[, normal_samples]
# created combined matrix with labels for PCA
combinedMatrix <- cbind(tumor_expression_matrix, normal_expression_matrix)
sampleLabels <- c(rep("T", ncol(tumor_expression_matrix)), rep("N", ncol(normal_expression_matrix)))

pca_result <- prcomp(tr_data, center = TRUE, scale = TRUE, retx = TRUE)
summary(pca_result)

plot(pca_result$x[, 1], pca_result$x[, 2], col = c(rep("red", ncol(tumor_expression_matrix)), rep("green", ncol(normal_expression_matrix))),  
     xlab = "PC1", ylab = "PC2", main = "PCA of GSE24514")

#plot(pca_result$x[, 1], pca_result$x[, 2], col = c("red", rep("blue", ncol(normal_expression_matrix))),xlab = "PC1", ylab = "PC2", main="PCA of GSE24514")
text(pca_result$x[, 1], pca_result$x[, 2], labels = rownames(tr_data), col = c(rep("red", 34), rep("green", ncol(normal_expression_matrix))), pos = 3)
legend("bottomleft", legend=c("Tumor", "Normal"), col=c("red", "green"), pch=1)


#Question 3: Benjamini-Hochberg Correction - t-test
p_val = NULL
for (i in 1:nrow(log_data)) {
  p_val[i] = t.test(log_data[i,1:34], log_data[i,35:49])$p.value
}

#Used the "p.adjust (p, method = “BH”)" 
#Benjamini-Hochberg correction cut-off of 0.05 to identify significantly changed genes.
corrected_p_values = p.adjust(p_val, method = "BH")
significant_genes <- which(corrected_p_values < 0.05)
sig_genes <- length(significant_genes) #8329

#Symbols of these genes
data_feature = data_set@featureData@data$'Gene Symbol'
BH_Sign_Genes_symbol = data_feature[significant_genes, 3] ###???

#3 genes-most significant p-values, Find the indices of the smallest three values
most_3 <- order(corrected_p_values)[1:3]
print(most_3) #3442  5773 17769 index of the Gene list
a <- data_set@featureData@data$`Gene symbol`
a[most_3] #most common 3 genes names; "CXCL9" "MICB"  "SNX10"

#Question 5: hclust

install.packages(c("gplots", "RColorBrewer"))
library(gplots)

data <- log_data[significant_genes,]

#We used to pearson correlation in as.dist function.
pg <- as.dist(1- cor(t(data), method = "pearson"))
ps <- as.dist(1-cor(data, method = "pearson"))

#After, used to hclust function and in dist_samples and centroid method
hc_pg <- hclust(pg, method = "centroid")
hc_ps <- hclust(ps, method = "centroid")

#we created the dendrogram using the as.dendrogram fucntion
#after, visualize the h_samples dendrogram
plot(as.dendrogram(hc_ps)) #compared to the samples

dg <- as.dendrogram(hc_pg)
ds <- as.dendrogram(hc_ps)

heatmap(data, Colv = ds, Rowv = dg, scale = "row")

#Question 7: Classification with SVM

install.packages("e1071")
library(e1071)

#tumor sample 0: 1:34, normal sample 1: 35:49: 
#train_data= 2:48, test data sample 1, and sample 49, tumor and normal samples

label <- c(rep(0, 34), rep(1 ,15))
frame = data.frame(label,t(data))

train_data = frame[2:48,]
test_data = frame[c(1,49),] 

svm_model <- svm(formula = label ~ ., 
                 data = train_data, 
                 type = 'C-classification', 
                 kernel = 'linear', scale = TRUE, cost=10)
summary(svm_model)

try_svm <- predict(svm_model, train_data) #control svm with train_data
print(table(try_svm))

predictions <- predict(svm_model, test_data) #test_data
predictions

#Question 8: GGM-Based Network Inference 

install.packages("ppcor")
library(ppcor)

sort_pvalue = sort(p_val)
most_Sg = which(p_val < sort_pvalue[6])
print(most_Sg) #3442  4866  5773  9038 17769
signgenes_symbol = data_set@featureData@data$`Gene symbol`[most_Sg] 
#"CXCL9" "STIL"  "MICB"  "RIPK2" "SNX10"

PPC_data = log_data[c(3442, 4866,  5773,  9038, 17769),] #PPC_data <- log_data[most_5,]

PCOR_Result <- pcor(t(PPC_data), "pearson")

#described the estimate and pvalue threshold our output
result_ggm <- ifelse((ifelse(abs(PCOR_Result$estimate) >= 0.35, 1, 0)) & (ifelse(abs(PCOR_Result$p.value) < 0.05, 1, 0))==TRUE,1,0) 


