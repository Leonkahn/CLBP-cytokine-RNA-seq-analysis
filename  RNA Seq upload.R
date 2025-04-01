install.packages("FactoMineR")
install.packages("factoextra")
install.packages("ggplot2")
library(factoextra)
library(FactoMineR)

head(Count_QC_rpkm_2)
nrow(Count_QC_rpkm_2)
ncol(Count_QC_rpkm_2)

#delete Individuals 111 and 115
class(Count_QC_rpkm_2)
Original_data_clean <- Count_QC_rpkm_2[, !colnames(Count_QC_rpkm_2) %in% c("B0034_111_1_BR", "B0034_115-1_BR")]
ncol(Original_data_clean)

#transpose data- switches individuals into rows and variables into columns
RNA_transposed <- t(Original_data_clean)
nrow(RNA_transposed)
ncol(RNA_transposed)

gene_names <- RNA_transposed[2,]
head(gene_names)
length(gene_names)


RNA_numeric <- RNA_transposed[-c(1, 2),] #get rid of row "id" and gene name bc they are not numeric  
nrow(RNA_numeric)
ncol(RNA_numeric)

all(sapply(RNA_numeric, is.numeric))
#convert data to numeric
RNA_numeric <- apply(RNA_numeric, 2, as.numeric)
str(RNA_numeric)
any(is.na(RNA_numeric))

#**Filter out low expressed genes:**
#calculate variance for each gene
gene_variances <- apply(RNA_numeric, 2, var)
summary(gene_variances)

#log transformation:
log_gene_variances <- log(gene_variances + 1)

#Plot a histogram 
hist(log_gene_variances, xlab = "Variance", ylab = "Frequency")
quantile(log_gene_variances, probs = seq(0, 1, 0.1))

#set a variance threshold
var_threshold <- 6.868996e-01  #is 50%, 0.5 variance 



#Filter genes with variance >= var_threshold
filtered_genes <- RNA_numeric[ ,log_gene_variances >= var_threshold]
nrow(filtered_genes)
ncol(filtered_genes)
colnames(RNA_numeric) <- gene_names  
colnames(filtered_genes) <- gene_names[log_gene_variances >= var_threshold]  

filtered_gene_names <- gene_names[log_gene_variances >= var_threshold]

filtered_gene_names <- gsub("\\.[0-9]+$", "", filtered_gene_names)  # Remove `.1`, `.2`, etc.
filtered_gene_names <- gsub("-[0-9]+$", "", filtered_gene_names)    # Remove `-1`, `-2`, etc.


#checking for missing values NA
any(is.na(filtered_genes))

#Centering and Scaling:
str(filtered_genes)
all(sapply(filtered_genes, is.numeric))

scaled_RNA <- scale(filtered_genes, center = TRUE, scale = TRUE)
is.data.frame(scaled_RNA)
scaled_RNA_df <- as.data.frame(scaled_RNA)
dim(scaled_RNA_df)

colnames(scaled_RNA_df) <- make.unique(gsub("\\.[0-9]+$", "", colnames(scaled_RNA_df)))
colnames(scaled_RNA_df) <- make.unique(gsub("-[0-9]+$", "", colnames(scaled_RNA_df)))



#Check for NA:
any(is.na(scaled_RNA_df))
sum(is.na(scaled_RNA_df))

#1.1 Histogram random variable 
DPM1 <- scaled_RNA_df$"DPM1"  
hist(DPM1, main = "Histogram of DPM1", xlab = "DPM1", col = "lightblue", border = "black")

GCLC <- scaled_RNA_df$GCLC
hist(GCLC, main = "Histogram of GCLC", xlab = "GCLC", col = "lightblue", border = "black")

#1.3
RNA_pca <- FactoMineR::PCA(scaled_RNA_df, graph = FALSE, scale.unit = FALSE) #FALSE bc i have already scaled 
print(RNA_pca)
summary(RNA_pca)

#extract PCA results:
pca_loadings <- RNA_pca$var$coord 
pca_loadings_df <- as.data.frame(pca_loadings)

rownames(pca_loadings_df) <- make.unique(gsub("\\.[0-9]+$", "", rownames(pca_loadings_df)))
rownames(pca_loadings_df) <- make.unique(gsub("-[0-9]+$", "", rownames(pca_loadings_df)))

any(duplicated(rownames(pca_loadings_df)))  #should say FALSE

head(pca_loadings_df)


#Visualization and Interpretation:#3
eig.val_RNA <- get_eigenvalue(RNA_pca)
eig.val_RNA 

sum(eig.val_RNA[,1]) 
ncol(scaled_RNA_df)

#check why sum(eig_value) not the same is as number columns 
#Calculate means of each variable
column_means <- colMeans(scaled_RNA_df)
head(column_means)


#scree plot:
fviz_eig(RNA_pca, addlabels = TRUE, ylim = c(0,65))

explained_variance <- eig.val_RNA[, 2]  
barplot(explained_variance, 
        main = "Proportion of Variance Explained", 
        xlab = "Principal Components", 
        ylab = "Variance (%)", 
        col = "skyblue")

cumulative_variance <- cumsum(explained_variance)
plot(cumulative_variance, pch = 19, col = "red", 
     main = "Cumulative Variance Explained", 
     xlab = "Principal Components",
     ylab = "Cumulative Variance (%)",font.main = 1)



#Graph of variables: #5
var <- get_pca_var(RNA_pca)
var

head(var$coord)
head(var$cor)
head(var$cos2) 
head(var$contrib)

dim(var$coord)    
nrow(var$contrib) 

#top 10 varibales contrib to Dimension 1
#extract contrib for dimension 1
contrib_dim1 <- var$contrib[, 1]  
contrib_df_dim1 <- data.frame(Variable = rownames(var$contrib), Contribution = contrib_dim1)
#sort data frame by contributions in descending order
contrib_df_sorted_dim1 <- contrib_df_dim1[order(-contrib_df_dim1$Contribution), ]
top_20_contrib_dim1 <- head(contrib_df_sorted_dim1, 20)
print(top_20_contrib_dim1)

contrib_dim2 <- var$contrib[, 2] 
contrib_df_dim2 <- data.frame(Variable = rownames(var$contrib), Contribution = contrib_dim2)
contrib_df_sorted_dim2 <- contrib_df_dim2[order(-contrib_df_dim2$Contribution), ]
top_20_contrib_dim2 <- head(contrib_df_sorted_dim2, 20)
print(top_20_contrib_dim2)

library(ggplot2)

#lollipop plot
ggplot(top_20_contrib_dim2, aes(x = reorder(Variable, Contribution), y = Contribution)) +
  geom_segment(aes(xend = Variable, yend = 0), color = "steelblue") +
  geom_point(color = "steelblue", size = 4) +
  coord_flip() +
  labs(title = "Top 20 genes contributing to Dim2", x = "Gene", y = "Contribution (%)") +
  theme_minimal()


fviz_pca_var(RNA_pca, 
             col.var = "black",
             select.var = list(contrib = 20))


library(corrplot)
subset_vars <- var$contrib[order(-var$contrib[, 1]), ][1:30, ]  
corrplot(subset_vars, is.corr = FALSE)

subset_vars <- var$contrib[order(-rowSums(var$contrib[, 1:5])), ][1:30, ]
corrplot(subset_vars, is.corr = FALSE)


install.packages("pheatmap")
library(pheatmap)
pheatmap(var$contrib[1:50, ], cluster_rows = TRUE, cluster_cols = TRUE,main = "Top 50 genes contribution to PCs", repel = TRUE)

fviz_pca_var(RNA_pca, col.var = "cos2",
             gradient.cols = c("yellow2", "orange", "red4"), select.var = list(cos2 = 20),
             repel = TRUE)


#Dimension description:
#dimdesc 
res.desc <- dimdesc(RNA_pca, axes = c(1, 2), proba = 0.05)

#top contributor for Dim1 and Dim2
top_contributor_dim1 <- head(res.desc$Dim.1$quanti, 20)  
top_contributor_dim2 <- head(res.desc$Dim.2$quanti, 20)  

print(top_contributor_dim1)
print(top_contributor_dim2)

#Extract gene names of top contributors
top_genes_dim1 <- rownames(top_contributor_dim1)
cat(top_genes_dim1, sep = "\n")

top_genes_dim2 <- rownames(top_contributor_dim2)
cat(top_genes_dim2, sep = "\n")

RNA_pca_scores <- as.data.frame(RNA_pca$ind$coord)
#Analysis of individuals:
ind_ID <- colnames(Original_data_clean)[-c(1, 2)] #exclude gene names and feature ID
head(ind_ID)

#add ind_id to RNA_pca_score
RNA_pca_scores$ind_ID <- ind_ID
head(RNA_pca_scores)
print(RNA_pca_scores)

#Graph of individuals:
ind <- get_pca_ind(RNA_pca)
ind

head(ind$coord)

#quality of individuals 
head(ind$cos2) 
print(ind$cos2)
#contributions of individuals 
head(ind$contrib)
print(ind$contrib)


#extract contributions of individuals
ind_contrib <- as.data.frame(RNA_pca$ind$contrib)
top_contrib_ind <- ind_contrib[order(-ind_contrib$Dim.1), ]  
top_20_contrib_ind <- head(top_contrib_ind, 20)
print(top_20_contrib_ind[, c("Dim.1", "Dim.2")])

top_20_contrib_table <- data.frame(Individual = rownames(top_20_contrib_ind), Dim.1 = top_20_contrib_ind$Dim.1)
print(top_20_contrib_table, row.names = FALSE)



library(ggplot2)
# visualize samples in PCA space:
ggplot(RNA_pca_scores, aes(x = Dim.1, y = Dim.2)) +
  geom_point(alpha = 0.5, color = "black") +  
  geom_density_2d(color = "blue") +
  theme_minimal()


library(factoextra)
fviz_pca_ind(RNA_pca, labelsize = 3, repel = TRUE) 
fviz_pca_ind(RNA_pca, select.ind = list(cos2 = 0.5), repel = TRUE, labelsize = 3) 
fviz_pca_ind(RNA_pca, select.ind = list(cos2 = 20)) 
fviz_pca_ind(RNA_pca, select.ind = list(contrib = 20))

fviz_cos2(RNA_pca, choice = "ind")  
fviz_contrib(RNA_pca, choice = "ind", axes = 1:2)  


library(ggplot2)
ggplot(pca_loadings_df, aes(x = Dim.1, y = Dim.2)) +
  geom_point(alpha = 0.3, color = "black") +  
  geom_density_2d(color = "blue") +  
  labs(title = "2D Density Estimation on PCA Results", x = "PC1", y = "PC2") +
  theme_minimal()

ggplot(pca_loadings_df, aes(x = Dim.1, y = Dim.2)) +
  geom_point(alpha = 0.3, color = "black") +
  geom_density_2d_filled(alpha = 0.6) +  # Filled density with color scale
  labs(title = "2D Density Estimation on PCA Results", x = "PC1", y = "PC2", fill = "Density") +
  theme_minimal()





#only numeric columns for clustering
numeric_data <- RNA_pca_scores[, sapply(RNA_pca_scores, is.numeric)]

library(factoextra)

#determine optimal number of clusters
fviz_nbclust(numeric_data, kmeans, method = "wss") +
  labs(title = "Elbow Method for Optimal Clusters")

fviz_nbclust(numeric_data, kmeans, method = "silhouette") +
  labs(title = "Silhouette Analysis for Optimal Clusters")



#k-means clustering
set.seed(123)  
kmeans_result <- kmeans(numeric_data, centers = 3)

RNA_pca_scores_temp <- RNA_pca_scores
RNA_pca_scores_temp$cluster <- as.factor(kmeans_result$cluster)

ggplot(RNA_pca_scores_temp, aes(x = Dim.1, y = Dim.2, color = cluster)) +
  geom_point(alpha = 0.5) +
  labs(title = "K-means Clustering on PCA Results", x = "PC1", y = "PC2") +
  theme_minimal()



#ACHTUNG do first the cohort comparisoon part to add cohort to PCA result!!!!
RNA_pca_scores_with_cohort$cluster <- as.factor(kmeans_result$cluster)
print(table(RNA_pca_scores_with_cohort$cluster, RNA_pca_scores_with_cohort$cohort))

#plot with clusters and cohort labels
ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.1, y = Dim.2, color = as.factor(cluster), shape = as.factor(cohort))) +
  geom_point(alpha = 0.7, size = 3) +
  labs(
    title = "K-means Clustering with Cohort Labels",
    x = "PC1",
    y = "PC2",
    color = "Cluster",
    shape = "Cohort ") + theme_minimal()



#Chi-Square test
chi_square_test <- chisq.test(contingency_table)
print(chi_square_test)

boxplot(Dim.1 ~ cluster, data = RNA_pca_scores_with_cohort, main = "PC1 Distribution by Cluster")
boxplot(Dim.2 ~ cluster, data = RNA_pca_scores_with_cohort, main = "PC2 Distribution by Cluster")

boxplot(Dim.1 ~ cluster,
         data = RNA_pca_scores_with_cohort,
         main = "",  # suppress default title
         xlab = "Cluster",
         ylab = "PC1",
         col = "lightblue")

title(main = "PC1 Distribution by Cluster", font.main = 1)



fviz_contrib(RNA_pca, choice = "ind", axes = 1) +
  labs(title = "Contributions to PC1", x = "Individuals", y = "Contribution")
fviz_contrib(RNA_pca, choice = "ind", axes = 2) +
  labs(title = "Contributions to PC2", x = "Individuals", y = "Contribution")


#Biplot:
fviz_pca_biplot(RNA_pca, select.var = list(contrib = 10), repel = TRUE)

fviz_pca_biplot(RNA_pca, 
                select.var = list(contrib = 5),  
                repel = TRUE,
                col.var = "black",    
                col.ind = "blue",     
                label = "var")       


#Cohort comparison:
install.packages("readxl")  
library("readxl")
RNAseq_cohort <- read_excel("/Users/leonkahn/Documents/Masterthesis/Data Excel/RNA Seq PBMC/RNAseq Cohort.xlsx")
head(RNAseq_cohort)

cohort_column <- RNAseq_cohort$cohort
nrow(RNA_pca_scores)
length(cohort_column)
str(cohort_column)
RNA_pca_scores_with_cohort <- cbind(cohort = cohort_column, RNA_pca_scores )
head(RNA_pca_scores_with_cohort)

#Visualization:
library(factoextra)
fviz_pca_biplot(RNA_pca, 
                repel = TRUE,
                select.var = list(contrib = 10),
                col.var = "green",  
                col.ind = as.factor(RNA_pca_scores_with_cohort$cohort),  
                palette = c("blue", "red"),  
                legend.title = list(ind = "Cohort", var = "Variables"),
                label = "var") +  
  scale_shape_manual(values = c(16, 16))  

fviz_pca_biplot(RNA_pca, 
                repel = TRUE,
                select.var = list(contrib = 5),
                col.var = "darkgrey",   
                col.ind = as.factor(RNA_pca_scores_with_cohort$cohort),  
                palette = c("blue", "red"),  
                legend.title = list(ind = "Cohort", var = "Variables"),
                label = "var",
                pointsize = 2            
) + 
  scale_shape_manual(values = c(16, 16)) +    
  guides(shape = "none") +                    
  ggtitle("PCA-Biplot Cohort Comparison") +  
  labs(color = "Cohort")                      




RNA_pca_scores_with_cohort$cohort <- as.factor(RNA_pca_scores_with_cohort$cohort)

ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.1, y = Dim.2, color = cohort)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim1 vs Dim2",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.1, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim1 vs Dim3",
       x = "PC1",
       y = "PC2",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()


ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.2, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim2 vs Dim3",
       x = "PC2",
       y = "PC3",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

#statistical testing:
t.test(Dim.1 ~ cohort, data = RNA_pca_scores_with_cohort)
t.test(Dim.2 ~ cohort, data = RNA_pca_scores_with_cohort)
t.test(Dim.3 ~ cohort, data = RNA_pca_scores_with_cohort)


#labeled with ID 
ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.2, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = ind_ID), size = 2, vjust = 2, hjust = 1) +
  labs(title = "PCA Plot: Dim2 vs Dim3",
       x = "Principal Component 2",
       y = "Principal Component 3",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

#Pearson correlation coefficient:  
#Dim1
cor(RNA_pca_scores_with_cohort$Dim.1, as.numeric(RNA_pca_scores_with_cohort$cohort))
#Dim2
cor(RNA_pca_scores_with_cohort$Dim.2, as.numeric(RNA_pca_scores_with_cohort$cohort))
#Dim3
cor(RNA_pca_scores_with_cohort$Dim.3, as.numeric(RNA_pca_scores_with_cohort$cohort))#moderate correlation 
#Dim4
cor(RNA_pca_scores_with_cohort$Dim.4, as.numeric(RNA_pca_scores_with_cohort$cohort))
#Dim5
cor(RNA_pca_scores_with_cohort$Dim.5, as.numeric(RNA_pca_scores_with_cohort$cohort))

#Pearson's product-moment correlation 
cor.test(RNA_pca_scores_with_cohort$Dim.1, as.numeric(RNA_pca_scores_with_cohort$cohort))

#example using MANOVA
manova_result_cohort<- manova(cbind(Dim.1, Dim.2, Dim.3) ~ cohort, data = RNA_pca_scores_with_cohort)
summary(manova_result_cohort)

#boxplot of distribution of Dim1 by cohort
ggplot(RNA_pca_scores_with_cohort, aes(x = cohort, y = Dim.1, color = cohort)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.1 by Cohort",
       x = "Cohort",
       y = "Principal Component 1") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.1, color = cohort)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.1 by Cohort",
       x = "Principal Component 1",
       y = "Density",
       color = "Cohort") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()



#Dim2 across cohort:
ggplot(RNA_pca_scores_with_cohort, aes(x = cohort, y = Dim.2, color = cohort)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.2 by Cohort",
       x = "Cohort",
       y = "Principal Component 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_cohort, aes(x = Dim.2, color = cohort)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.2 by Cohort",
       x = "Principal Component 2",
       y = "Density",
       color = "Cohort") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()




#Adding Gender to data:

gender_column <- RNAseq_cohort$gender
RNA_pca_scores_with_gender <- cbind(gender = gender_column, RNA_pca_scores )
head(RNA_pca_scores_with_gender)

#both added
RNA_pca_scores_with_cohort_gender <- cbind(gender = gender_column, RNA_pca_scores_with_cohort )
head(RNA_pca_scores_with_cohort_gender)
print(RNA_pca_scores_with_cohort_gender)



#analysis of PCA compared to gender:
RNA_pca_scores_with_gender$gender <- as.factor(RNA_pca_scores_with_gender$gender)

str(RNA_pca_scores_with_gender)
ggplot(RNA_pca_scores_with_gender, aes(x = Dim.1, y = gender, color = gender)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim.1 by Gender",
       x = "Principal Component 1",
       y = "Gender",
       color = "Gender") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_gender, aes(x = gender, y = Dim.1, color = gender)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.1 by Sex",
       x = "Sex",
       y = "Principal Component 1") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_gender, aes(x = Dim.1, color = gender)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.1 by Gender",
       x = "Principal Component 1",
       y = "Density",
       color = "Gender") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()


ggplot(RNA_pca_scores_with_gender, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = gender), alpha = 0.7) +
  facet_wrap(~ gender) +
  labs(title = "PCA Plot by Gender",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()



#t-test:
t.test(Dim.1 ~ gender, data = RNA_pca_scores_with_gender)
t.test(Dim.2 ~ gender, data = RNA_pca_scores_with_gender)
t.test(Dim.3 ~ gender, data = RNA_pca_scores_with_gender)




#example using MANOVA
manova_result_gender <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ gender, data = RNA_pca_scores_with_gender)
summary(manova_result_gender)



#adding AGE to pca result
age_column <- RNAseq_cohort$Alter

RNA_pca_scores_with_cohort_gender_age <- cbind(gender = gender_column, age = age_column, RNA_pca_scores_with_cohort)
head(RNA_pca_scores_with_cohort_gender_age)
str(RNA_pca_scores_with_cohort_gender_age)

RNA_pca_scores_with_cohort_gender_age$age <- as.numeric(RNA_pca_scores_with_cohort_gender_age$age)

head(RNA_pca_scores_with_cohort_gender_age)


#correlation coefficeint
cor(RNA_pca_scores_with_cohort_gender_age$Dim.1, RNA_pca_scores_with_cohort_gender_age$age)
cor(RNA_pca_scores_with_cohort_gender_age$Dim.2, RNA_pca_scores_with_cohort_gender_age$age)
cor(RNA_pca_scores_with_cohort_gender_age$Dim.3, RNA_pca_scores_with_cohort_gender_age$age)


#MANOVA 
manova_result_age <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ age, data = RNA_pca_scores_with_cohort_gender_age)
summary(manova_result_age)


ggplot(RNA_pca_scores_with_cohort_gender_age, aes(x = age, y = Dim.1, color = cohort)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", se = FALSE, aes(color = cohort), linetype = "dashed") +  
  labs(title = "Association Between Age and PC1 by Cohort",
       x = "Age (years)",
       y = "PC1",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_cohort_gender_age, aes(x = age, y = Dim.2, color = cohort)) +
  geom_point(alpha = 0.7) +  
  geom_smooth(method = "lm", se = FALSE, aes(color = cohort), linetype = "dashed") +  
  labs(title = "Association Between Age and PC2 by Cohort",
       x = "Age (years)",
       y = "PC2",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

ggplot(RNA_pca_scores_with_cohort_gender_age, aes(x = age, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +  
  geom_smooth(method = "lm", se = FALSE, aes(color = cohort), linetype = "dashed") +  
  labs(title = "Association Between Age and PC3 by Cohort",
       x = "Age (years)",
       y = "PC3",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()


#Creating a table for overview:
summary_table <- data.frame(
  id = RNA_pca_scores_with_cohort_gender$ind_ID,  
  Dim.1 = RNA_pca_scores_with_cohort_gender$Dim.1,
  Cohort = RNA_pca_scores_with_cohort_gender$cohort)

# View the table
print(summary_table)





