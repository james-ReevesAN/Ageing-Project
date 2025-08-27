####Call Ageing genes in nSolver genenames
Agestr <- c("BACH2","TIGIT","KLRF1","PRF1","SESN3","LEF1","CCR7","B3GAT1","COX16","MYC", "GZMB" , "CD79A" , "FOXP3" , "CXCR4" , 'HLA-DPA1', "CX3CR1" , "CD27" , "CX3CR1" , "NT5E" , "CD163" , "TNFRSF1B" , "IL15", "NLRP3", "BTLA", "IFNG" , "IL7R" , "CD160",
            "RORC", "CD97", "CD40LG" , "IRAK1", "ITGA6" , "IRF1" , "CCL13" , "BATF")

Age_df <- gene_data[gene_data$Name %in% Agestr, ]
remove(Agestr)

#AGEing genes not found #######
Agestr[ !Agestr %in% rownames(Age_df)] 

#Read in csv 
gene_data <- read_csv("~/R_Projects/Relative gene Expresssion .csv")

gene_data <- as.data.frame(gene_data)

NormCounts.PCA <- as.data.frame(gene_data[3:12])
rownames(NormCounts.PCA) <- as.character(gene_data$NanoID)
# To investigate samples for ten-Gene panel
NormCounts.PCA <- NormCounts.PCA[18:47,]

Age_df %>% 
  filter(colnames(Age_df) rownames(NormCounts.PCA))

Age_df <- Age_df[,colnames(Age_df) %in% rownames(NormCounts.PCA) ]

Age_df1 <- left_join(Age_df_t, NormCounts.PCA, by = NanoID)

NormCounts.PCA <- NormCounts.PCA[rownames(NormCounts.PCA) != "20211103_AN Set9_D068_10.RCC", ]

                                               
Age_df_t <- as.data.frame(t(Age_df))
NormCounts.PCA$NanoID <- rownames(NormCounts.PCA)
Age_df_t$NanoID <- rownames(Age_df_t)
  
Age_df_t <- Age_df_t[,1:5]  

Age_df_t[, 1:5] <- scale(Age_df_t[, 1:5])

######Plot to compare Nanostring vs qPCR dqta
# Select relevant columns + NanoID
df_plot <- Age_df_t %>%
  select(NanoID, LEF1, KLRF1, TIGIT, CCR7, PRF1)

# Reshape to long format for ggplot
df_long <- df_plot %>%
  pivot_longer(cols = c(LEF1, KLRF1, TIGIT, CCR7, PRF1),
               names_to = "Gene",
               values_to = "Expression")

# Bar plot
ggplot(df_long, aes(x = NanoID, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  #scale_y_log10() +
  labs(title = "Gene Expression per Sample",
       x = "Sample (NanoID)",
       y = "Expression",
       fill = "Gene") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(df_long, aes(x = NanoID, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity") +
  #scale_y_log10() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Gene Expression per Sample (Log Scale)",
       x = "Sample (NanoID)",
       y = "Log10(Expression)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

NormCounts.PCA_t <- NormCounts.PCA[, c("CCR7", "KLRF1", "LEF1", "PRF1", "TIGIT")]



df_plot <- NormCounts.PCA %>% 
  select(NanoID, LEF1, BACH2, TIGIT, KLRF1 , PRF1 , SESN3 ,LEF1 , CCR7 ,B3GAT1 ,COX16)

# Reshape to long format for ggplot
df_long <- df_plot %>%
  pivot_longer(cols = c(LEF1, BACH2, TIGIT, KLRF1 , PRF1 , SESN3 ,LEF1 , CCR7 ,B3GAT1 ,COX16),
               names_to = "Gene",
               values_to = "Expression")


#Run UMAP Algorithmn
NormCounts.PCA_t  |>
  uwot::umap(
    n_neighbors = 6, n_threads = (parallel::detectCores()-2),
    min_dist = 0.01,
    n_sgd_threads = (parallel::detectCores()-2), 
    verbose = TRUE
  ) -> umap_exprn

umap_exprn <- data.table::as.data.table(umap_exprn)
colnames(umap_exprn) <- c("UMAP_X", "UMAP_Y")

#cbind umap values to exprn data.table
NormCounts.PCA_t <- cbind(NormCounts.PCA_t,umap_exprn)

remove(umap_exprn)

###########UMAP Plots

ggplot(gene_data_NS, aes(x = UMAP_X, y = UMAP_Y, color = TIGIT)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient2(
    low = "#FFA46B",       # Low scores
    mid = "#F1EAE4",      # Midpoint (score = 0)
    high = "#5830F7",       # High scores
    midpoint = 0        # This sets 0 as the central color
  ) +
  theme_minimal() 

ggplot(gene_data_QP, aes(x = UMAP_X, y = UMAP_Y, color = SurvivalFive)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Yes" = "green3", "No" = "orangered")) +
  theme_minimal()

ggplot(gene_data_QP, aes(x = UMAP_X, y = UMAP_Y, color = RelapseFive)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("No" = "green3", "Yes" = "orangered")) +
  theme_minimal()

ggplot(NormCounts.PCA_t, aes(x = UMAP_X, y = UMAP_Y, color = NanoID)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()


ggplot(join, aes(x = UMAP_X, y = UMAP_Y, color = Age)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient2(
    low = "#FFA46B",       # Low scores
    mid = "#F1EAE4",      # Midpoint (score = 0)
    high = "#5830F7",       # High scores
    midpoint = 30        # This sets 0 as the central color
  ) +
  theme_minimal() 


ggplot(gene_data_NS, aes(x = UMAP_X, y = UMAP_Y, color = DonorSex)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal()


gene_dataNS <-  gene_data[gene_data$NanoID %in% rownames(Age_df_t), ]
gene_data_QP <- gene_data[gene_data$NanoID %in% rownames(NormCounts.PCA_t), ]

gene_data_NS <- cbind(gene_dataNS, Age_df_t[,7:8])
  
gene_data_QP <- cbind(gene_data_QP, NormCounts.PCA_t[,6:7])








