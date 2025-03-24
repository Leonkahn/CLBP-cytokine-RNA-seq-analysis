install.packages("FactoMineR")
install.packages("factoextra")
install.packages("ggplot2")
library(factoextra)
library(FactoMineR)

head(cytokine_data)
nrow(cytokine_data)
ncol(cytokine_data)


install.packages("dplyr")
library(dplyr)

#Count cohort members
cytokine_data %>%
  count(cohort)


cytokine_data_NA <- cytokine_data[,-c(1, 2)] #get rid of column "id" bc is not numeric 
head(cytokine_data_NA)
nrow(cytokine_data_NA)

#1.2 log transformation 
cytokine.log <- log(cytokine_data_NA + 1 ) #+1 to avoid taking log of 1

#handling missing values NA
dim(na.omit(cytokine_data_NA)) 
install.packages("VIM")
install.packages("naniar")
install.packages("missMDA")

print(sum(is.na(cytokine.log)))

library(naniar)
gg_miss_var(cytokine.log) #plot of NA

library(VIM)
resVIM<-summary(aggr(cytokine.log, sortVar=TRUE))$combinations

matrixplot(cytokine.log, sortby = 2)


library(ggplot2)
vis_miss(cytokine.log) +
  labs(x = "Variables", y = "Individuals") +
  scale_fill_manual(values = c("gray80", "red"),
                    labels = c("Present", "Missing")) +
  guides(fill = guide_legend(title = "Missing Data")) 



#PCA with missing values
library(missMDA)
nb <- estim_ncpPCA(cytokine.log,method.cv = "Kfold", verbose = FALSE) 
print(nb$ncp)

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")
res.comp <- imputePCA(cytokine.log, ncp = nb$ncp) 
cytokine_data_imputed <- res.comp$completeObs
print(head(cytokine_data_imputed))
print(cytokine_data_imputed[42:53])

class(cytokine_data_imputed)
cytokine_data_imputed_df <- as.data.frame(cytokine_data_imputed)
class(cytokine_data_imputed_df)

#check why/if negative values after imputation
has_negative_values <- sum(cytokine_data_imputed_df < 0, na.rm = TRUE)
print(has_negative_values)


#adding cohort back into data frame
if (nrow(cytokine_data) == nrow(cytokine_data_imputed_df)) {
  cytokine_data_with_cohort <- cbind(cytokine_data$cohort, cytokine_data_imputed_df)
  colnames(cytokine_data_with_cohort)[1] <- "Cohort"
  
  print(head(cytokine_data_with_cohort))}


#Descriptive Statistics: 

#data separated into cohorts
df_cohort_0 <- subset(cytokine_data_with_cohort, Cohort == 0)
df_cohort_1 <- subset(cytokine_data_with_cohort, Cohort == 1)
head(df_cohort_0)
head(df_cohort_1)


#values calculated separately per Cohort:
calculate_descriptive_stats <- function(data) {
  numeric_vars <- data[, sapply(data, is.numeric)]
  
  mean_values <- colMeans(numeric_vars, na.rm = TRUE)
  median_values <- sapply(numeric_vars, median, na.rm = TRUE)
  
  get_mode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  mode_values <- sapply(numeric_vars, get_mode)
  
  sd_values <- sapply(numeric_vars, sd, na.rm = TRUE)
  variance_values <- sapply(numeric_vars, var, na.rm = TRUE)
  min_values <- sapply(numeric_vars, min, na.rm = TRUE)
  max_values <- sapply(numeric_vars, max, na.rm = TRUE)
  range_values <- max_values - min_values
  iqr_values <- sapply(numeric_vars, IQR, na.rm = TRUE)
  
  absolute_deviation_values <- sapply(numeric_vars, function(x) {
    mean(abs(x - mean(x, na.rm = TRUE)), na.rm = TRUE)
  })
  
  data.frame(
    Mean = mean_values,
    Median = median_values,
    Mode = mode_values,
    `Standard Deviation` = sd_values,
    Variance = variance_values,
    Minimum = min_values,
    Maximum = max_values,
    Range = range_values,
    IQR = iqr_values,
    `Absolute Deviation` = absolute_deviation_values
  )
}

#statistic for complete dataset 
descriptive_stats_comb <- calculate_descriptive_stats(cytokine_data_with_cohort)
round(descriptive_stats_comb, 3)

#statistics for Cohort 0
descriptive_stats_cohort_0 <- calculate_descriptive_stats(df_cohort_0)
round(descriptive_stats_cohort_0,3)

#statistics for Cohort 1
descriptive_stats_cohort_1 <- calculate_descriptive_stats(df_cohort_1)
round(descriptive_stats_cohort_1, 3)


#1.1 einzelne Histogram
MMP1 <- cytokine_data_imputed_df$MMP1
hist(MMP1)

NPY <- cytokine_data_imputed_df$NPY
hist(NPY)

PDGF_BB <- cytokine_data_imputed_df$`PDGF-BB`
hist(PDGF_BB, main = "Histogram of PDGF-BB", col = "lightblue", border = "black", xlab = "PDGF-BB")

#using lapply to create histograms
par(mfrow = c(3, 6), mar = c(2, 2, 2, 2)) 

invisible(lapply(names(cytokine_data_imputed_df), function(var_name) {
  hist(cytokine_data_imputed_df[[var_name]], main = paste(" ", var_name), col = "lightblue", border = "black", xlab = var_name)}))

par(mfrow = c(1, 1)) #reset plotting parameter

# PCA
cytokine.pca <- FactoMineR::PCA(cytokine_data_imputed_df, graph = FALSE, scale.unit = TRUE) 
print(cytokine.pca)#2
summary(cytokine.pca)

#Visualization and Interpretation:#3
library(factoextra)
eig.val_cytokine <- get_eigenvalue(cytokine.pca)
eig.val_cytokine #ask 
sum(eig.val_cytokine[,1]) 

#scree plot:#4
fviz_eig(cytokine.pca, addlabels = TRUE, ylim = c(0,30))
fviz_screeplot(cytokine.pca, addlabels= TRUE, ylim = c(0,30)) #same as above

#Graph of variables:
var <- get_pca_var(cytokine.pca)#5
var

head(var$coord)
head(var$cor)
head(var$cos2) 
head(var$contrib)

print(var$coord)
print(var$cor)
print(var$cos2)
print(var$contrib)

#6
fviz_pca_var(cytokine.pca, col.var = "black",
             repel = TRUE)

library("corrplot") #quality of representation 
plot.new()
dev.off()
corrplot(var$cos2, is.corr = FALSE) #visualization of cos2 , is.corr = FALSE means no correlation matrix

#7
fviz_cos2(cytokine.pca, choice = "var", axes = 1:2) #bar plot
#total cos2 of variables on Dim.1 and Dim.2 
fviz_pca_var(cytokine.pca, col.var = "cos2",
             gradient.cols = c("yellow2", "orange", "red4"),
             repel = TRUE) #avoid text overlapping

fviz_pca_var(cytokine.pca, col.var = "cos2",
             gradient.cols = c("yellow2", "orange", "red4"), select.var = list(cos2 = 4), #top 4
             repel = TRUE)

fviz_pca_var(cytokine.pca, alpha.var = "cos2") #just transparency of arrows according to cos2 value 


#Contributions of variables to PCs:#8
head(var$contrib, 4)#The larger the value of the contribution, the more the variable contributes to the component.
library("corrplot")
plot.new()
dev.off()
corrplot(var$contrib, is.corr = FALSE)

fviz_contrib(cytokine.pca, choice = "var", axes = 1, top = 10) #Dim 1
fviz_contrib(cytokine.pca, choice = "var", axes = 2, top = 10) #Dim 2
fviz_contrib(cytokine.pca, choice = "var", axes = 1:2, top = 10) #both 
fviz_pca_var(cytokine.pca, col.var = "contrib",              
             gradient.cols = c("yellow","orange", "blue")) #highlighting the important ones

fviz_pca_var(cytokine.pca, alpha.var = "contrib") #display in transparency



#kmeans clustering
#create 3 groups 
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25) #three clusters, and 25 tries
grp <- as.factor(res.km$cluster)


fviz_pca_var(cytokine.pca, col.var = grp,           
             palette = c("darkblue", "yellow", "brown"),    
             legend.title = "Cluster")

#11
#dimension description: 
res.desc <- dimdesc(cytokine.pca, axes = c(1,2), proba = 0.05)
#description of dimension 1
res.desc$Dim.1
res.desc$Dim.2

#12
#graph of individuals:
ind <- get_pca_ind(cytokine.pca)
ind
#coordinates of individuals 
head(ind$coord)

#quality of individuals 
head(ind$cos2) 

#contributions of individuals 
head(ind$contrib)


#Plots: quality and contribution
fviz_pca_ind(cytokine.pca)
fviz_pca_ind(cytokine.pca, select.ind = list(cos2 = 0.5)) 
fviz_pca_ind(cytokine.pca, select.ind = list(cos2 = 20)) 
fviz_pca_ind(cytokine.pca, select.ind = list(contrib = 20)) 
fviz_pca_ind(cytokine.pca, select.ind = list(name = c("3","25","29"))) 
  
fviz_pca_ind(cytokine.pca, col.ind = "cos2",
             gradient.cols = c("green","yellow","darkviolet"),
             repel = TRUE)

fviz_pca_ind(cytokine.pca, pointsize = "cos2",     
             pointshape = 21, fill = "lightblue",      
             repel = TRUE)

fviz_pca_ind(cytokine.pca, col.ind = "cos2", pointsize = "contrib",     
             gradient.cols = c("red", "green", "violet"),       
             repel = TRUE)

fviz_cos2(cytokine.pca, choice = "ind")  
fviz_contrib(cytokine.pca, choice = "ind", axes = 1:2)  


#Graph customization:#13
#Variables on dimesnions 2 and 3
fviz_pca_var(cytokine.pca, axes = c(2,3))
#Individuals on dimesnions 2 and 3
fviz_pca_ind(cytokine.pca, axes = c(2,3))

#show variable point and text labels
fviz_pca_var(cytokine.pca, geom.var = c("point", "text"))
#show individuals text label only
fviz_pca_var(cytokine.pca, geom.ind = "text")

#customize labelsize 
fviz_pca_var(cytokine.pca, arrowsize = 0.8, labelsize = 5,
             repel = TRUE)
fviz_pca_ind(cytokine.pca,
             pointsize= 3, pointshape=21, fill="lightblue",
             labelsize=5, repel= TRUE)

#Biplot:#14
fviz_pca_biplot(cytokine.pca, repel = TRUE,
                col.var = "lightgreen",
                col.ind = "orange")





#Add ID back to PCA-scores:
pca_scores <- as.data.frame(cytokine.pca$ind$coord)

id_column <- cytokine_data$id
pca_results_with_id <- cbind(id = id_column, pca_scores)
head(pca_results_with_id)


#Add Cohort back to PCA-scores:
cohort_column <- cytokine_data$cohort
pca_results_with_id_cohort <- cbind(cohort = cohort_column, pca_results_with_id )
head(pca_results_with_id_cohort)


#Visualization:
library(factoextra)
fviz_pca_biplot(cytokine.pca, 
                repel = TRUE,
                axes = c(1,3),
                col.var = "darkgrey",  
                col.ind = pca_results_with_id_cohort$cohort,  
                legend.title = list(ind = "Cohort", var = "Variables"),
                palette = c("blue", "red"))  




pca_results_with_id_cohort$cohort <- as.factor(pca_results_with_id_cohort$cohort)

ggplot(pca_results_with_id_cohort, aes(x = Dim.1, y = Dim.2, color = cohort)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim1 vs Dim2",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()


ggplot(pca_results_with_id_cohort, aes(x = Dim.1, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim1 vs Dim3",
       x = "Principal Component 1",
       y = "Principal Component 3",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()


ggplot(pca_results_with_id_cohort, aes(x = Dim.2, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim2 vs Dim3",
       x = "Principal Component 2",
       y = "Principal Component 3",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

#statistical testing:
t.test(Dim.1 ~ cohort, data = pca_results_with_id_cohort)
t.test(Dim.2 ~ cohort, data = pca_results_with_id_cohort)
t.test(Dim.3 ~ cohort, data = pca_results_with_id_cohort)

#further investigation Dim3
ggplot(pca_results_with_id_cohort, aes(x = Dim.3, color = cohort)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.3 by Cohort",
       x = "Principal Component 3",
       y = "Density",
       color = "Cohort") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(pca_results_with_id_cohort, aes(x = cohort, y = Dim.3, color = cohort)) +
  geom_boxplot() +
  labs(title = "Boxplot of Dim.3 by Cohort",
       x = "Cohort",
       y = "Principal Component 3",
       color = "Cohort") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()




#Exclude the outlier
filtered_data <- pca_results_with_id_cohort %>% 
  filter(!(cohort == 1 & Dim.3 > 4))  
#Create the box plot without the outlier
ggplot(filtered_data, aes(x = cohort, y = Dim.3)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.3 by Cohort (Outlier Excluded)",
       x = "Cohort",
       y = "Principal Component 3") +
  theme_minimal()



#labeled with ID 

ggplot(pca_results_with_id_cohort, aes(x = Dim.2, y = Dim.3, color = cohort)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = id), size = 2, vjust = 1, hjust = 1) +
  labs(title = "PCA Plot: Dim2 vs Dim3",
       x = "Principal Component 2",
       y = "Principal Component 3",
       color = "Cohort") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

  
  
#Pearson correlation coefficient:  
#Dim1
cor(pca_results_with_id_cohort$Dim.1, as.numeric(pca_results_with_id_cohort$cohort))
#Dim2
cor(pca_results_with_id_cohort$Dim.2, as.numeric(pca_results_with_id_cohort$cohort))
#Dim3
cor(pca_results_with_id_cohort$Dim.3, as.numeric(pca_results_with_id_cohort$cohort))#moderate correlation 
#Dim4
cor(pca_results_with_id_cohort$Dim.4, as.numeric(pca_results_with_id_cohort$cohort))
#Dim5
cor(pca_results_with_id_cohort$Dim.5, as.numeric(pca_results_with_id_cohort$cohort))

#correlation for Dim3
cor.test(pca_results_with_id_cohort$Dim.3, as.numeric(pca_results_with_id_cohort$cohort))

#MANOVA example
manova_result_cohort<- manova(cbind(Dim.1, Dim.2, Dim.3) ~ cohort, data = pca_results_with_id_cohort)
summary(manova_result_cohort)



#k-means clsutering: determine hwo many clusters 
set.seed(123)
numeric_data <- pca_results_with_id_cohort[, sapply(pca_results_with_id_cohort, is.numeric)]

fviz_nbclust(numeric_data, kmeans, method = "wss") +
  labs(title = "Elbow Method for Optimal Clusters")

fviz_nbclust(numeric_data, kmeans, method = "silhouette") +
  labs(title = "Silhouette Analysis for Optimal Clusters")

#K-means clustering
set.seed(123) 
kmeans_result <- kmeans(numeric_data, centers = 2)
pca_results_with_id_cohort$cluster <- as.factor(kmeans_result$cluster)


ggplot(pca_results_with_id_cohort, aes(x = Dim.1, y = Dim.2, color = cluster)) +
  geom_point(alpha = 0.5) +
  labs(title = "K-means clustering on PCA results of cytokines", x = "PC1", y = "PC2") +
  theme_minimal()

#with cohort Labels 
ggplot(pca_results_with_id_cohort, aes(x = Dim.1, y = Dim.2, color = as.factor(cluster), shape = as.factor(cohort))) +
  geom_point(alpha = 0.7, size = 3) +
  scale_shape_manual(values = c(16, 17)) +  # 16 = circle (Cohort 0), 17 = triangle (Cohort 1)
  scale_color_manual(values = c("1" = "red", "2" = "green")) +  # Cluster 1 = Red, Cluster 2 = Green
  labs(
    title = "K-means Clustering of Cytokines with Cohort Labels ",
    x = "PC1",
    y = "PC2",
    color = "Cluster",
    shape = "Cohort"
  ) +
  theme_minimal()


pca_results_with_id_cohort$cluster <- as.factor(kmeans_result$cluster)

#contingency table
contingency_table_cytokine <- table(pca_results_with_id_cohort$cluster, pca_results_with_id_cohort$cohort)
print(contingency_table_cytokine)


#Chi-Square test
chi_square_test_cytokine <- chisq.test(contingency_table_cytokine)
print(chi_square_test_cytokine)

#Fisher test:
fisher.test(contingency_table_cytokine)



#Adding Gender to data: ADD EXCEL
class(Cytokine_Age_Gender)
Gender_cytokine<- Cytokine_Age_Gender[ ,-c(2,3)]

install.packages("dplyr")
library(dplyr)
pca_results_gender<- pca_results_with_id_cohort %>%
  left_join(Gender_cytokine, by = "id")
head(pca_results_gender)

library(ggplot2)
pca_results_gender$gender <- as.factor(pca_results_gender$gender)


#count gender 0 & gender 1
table(pca_results_gender$gender)
table(pca_results_gender$cohort)

#cohort1
gender_in_cohort1 <- pca_results_gender %>%
  filter(cohort == 1) %>%  
  count(gender)  

print(gender_in_cohort1)

#cohort 0
gender_in_cohort0 <- pca_results_gender %>%
  filter(cohort == 0) %>%  
  count(gender) 

print(gender_in_cohort0)


#scatter plot of Dim.2 vs. gender 
str(pca_results_gender)
ggplot(pca_results_gender, aes(x = Dim.2, y = gender, color = gender)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA Plot: Dim.2 by Gender",
       x = "Principal Component 2",
       y = "Gender",
       color = "Gender") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(pca_results_gender, aes(x = gender, y = Dim.2, color = gender)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.2 by Gender",
       x = "Gender",
       y = "Principal Component 2") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()



ggplot(pca_results_gender, aes(x = Dim.1, color = gender)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.1 by Gender",
       x = "Principal Component 1",
       y = "Density",
       color = "Gender") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(pca_results_gender, aes(x = Dim.2, color = gender)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.2 by Gender",
       x = "Principal Component 2",
       y = "Density",
       color = "Gender") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()

ggplot(pca_results_gender, aes(x = Dim.3, color = gender)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Dim.3 by Gender",
       x = "Principal Component 3",
       y = "Density",
       color = "Gender") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()


library(factoextra)
fviz_pca_biplot(cytokine.pca, 
                repel = TRUE,
                col.var = "orange",  
                col.ind = pca_results_gender$gender,
                legend.title = list(ind = "gender", var = "Variables"),
                palette = c("blue", "red"))  


ggplot(pca_results_gender, aes(x = Dim.2, y = Dim.3)) +
  geom_point(aes(color = gender), alpha = 0.7) +
  facet_wrap(~ gender) +
  labs(title = "PCA Plot by Gender",
       x = "Principal Component 2",
       y = "Principal Component 3") +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  theme_minimal()


#t-test

t.test(Dim.1 ~ gender, data = pca_results_gender)
ggplot(pca_results_gender, aes(x = gender, y = Dim.1)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.1 by Gender",
       x = "Gender",
       y = "Principal Component 1") +
  theme_minimal()


t.test(Dim.2 ~ gender, data = pca_results_gender)
ggplot(pca_results_gender, aes(x = gender, y = Dim.2)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.2 by Gender",
       x = "Gender",
       y = "Principal Component 2") +
  theme_minimal()


t.test(Dim.3 ~ gender, data = pca_results_gender)
ggplot(pca_results_gender, aes(x = gender, y = Dim.3)) +
  geom_boxplot() +
  labs(title = "Distribution of Dim.3 by Gender",
       x = "Gender",
       y = "Principal Component 3") +
  theme_minimal()

#MANOVA
manova_result_gender <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ gender, data = pca_results_gender)
summary(manova_result_gender)


#ANOVA 
anova_dim1 <- aov(Dim.1 ~ gender, data = pca_results_gender)
anova_dim2 <- aov(Dim.2 ~ gender, data = pca_results_gender)
anova_dim3 <- aov(Dim.3 ~ gender, data = pca_results_gender)


summary(anova_dim1)
summary(anova_dim2)
summary(anova_dim3)
print(anova_dim1)



#adding age to PCA results
Age_cytokine <- Cytokine_Age_Gender[, c("id", "age")]

pca_results_age <- pca_results_with_id_cohort %>%
  left_join(Age_cytokine, by = "id")

# Rename 'age' column to 'Age'
colnames(pca_results_age)[colnames(pca_results_age) == "age"] <- "Age"

# Remove redundant cohort columns after the join
pca_results_age <- pca_results_age[ , !grepl("cohort", colnames(pca_results_age))]

head(pca_results_age)
print(pca_results_age)


#Age distribution
pca_results_age <- pca_results_with_id_cohort %>%
  left_join(Age_cytokine, by = "id") %>%
  rename(Age = age) %>%
  select(id, cohort, Age, starts_with("Dim."))

library(dplyr)
colnames(pca_results_age)

age_summary <- pca_results_age %>%
  group_by(cohort) %>%
  summarise(
    mean_age = round(mean(Age, na.rm = TRUE), 1),
    sd_age = round(sd(Age, na.rm = TRUE), 1),
    n = n())

print(age_summary)


#correlation coefficient 
cor(pca_results_age$Dim.1, pca_results_age$Age)
cor(pca_results_age$Dim.2, pca_results_age$Age)
cor(pca_results_age$Dim.3, pca_results_age$Age)

#plots of age vs Dimensions
ggplot(pca_results_age, aes(x = Age, y = Dim.1)) + geom_point() + ggtitle("Age vs Dim.1")
ggplot(pca_results_age, aes(x = Age, y = Dim.2)) + geom_point() + ggtitle("Age vs Dim.2")
ggplot(pca_results_age, aes(x = Age, y = Dim.3)) + geom_point() + ggtitle("Age vs Dim.3")




#summary 
summary_table <- data.frame(
  id = pca_results_with_id_cohort$id,  
  Dim.1 = pca_results_with_id_cohort$Dim.1,
  Cohort = pca_results_with_id_cohort$cohort)

print(summary_table)
print(head(pca_results_with_id_cohort, 20))



install.packages("sessioninfo")  
library(sessioninfo)
session_info()






