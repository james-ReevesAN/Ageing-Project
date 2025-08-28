#Filepaths

nsfilepath <- "/Users/jamesreeves/R_Projects/NanoData/Donors"
sampleTabpath <- "//Users/jamesreeves/R_Projects/NanoData/Groups-Donors.csv"


# Define the directory containing .rcc files
folder_path <- "/Users/jamesreeves/R_Projects/NanoData/Donors"


#Normalize by Pos. and HK scale factors and remove all genes under background noise level (mean +2sd of neg.)
DataDefault.nSolver <- processNanostringData(nsFiles = nsfilepath,
                                             sampleTab = sampleTabpath,
                                             idCol = "SAMPLE", groupCol = "SurvivalFive", normalization = "nSolver", bgType = "t.test")

#Split dataset into young & old
young <- clinical_data[clinical_data$DonorAge < 30,]
old <- clinical_data[clinical_data$DonorAge > 40,]

split <- rbind(young,old)

#Subset gene data based on young/old splits
young_genes <- cbind(gene_data$Name,gene_data[colnames(gene_data) %in% young$SAMPLE | old$SAMPLE])
old_genes <- cbind(gene_data$Name,gene_data[colnames(gene_data) %in% old$SAMPLE])
split_genes <- cbind(gene_data$Name,gene_data[colnames(gene_data) %in% split$SAMPLE])

young_genes <- DataDefault.nSolver[DataDefault.nSolver$SAMPLE %in% young$SAMPLE]
old_genes <- DataDefault.nSolver[DataDefault.nSolver$SAMPLE %in% old$SAMPLE]
split_genes <- DataDefault.nSolver[DataDefault.nSolver$SAMPLE %in% split$SAMPLE]

res.limma.nsolver <- NanoTube::runLimmaAnalysis(split_genes , base.group = "Old" )
DFexpression.nSolver <- topTable(res.limma.nsolver, number = nrow(res.limma.nsolver), coef = 2)

#Get the lis of DEGs
DF<- DFexpression.nSolver %>% filter (P.Value < 0.05)
DF<-DF[,2]

DFexpression.nSolver$significance <- NA_integer_

# Add significance and -log10(p-value)
significance <- DFexpression.nSolver %>%
  mutate(
    significance = case_when(
      P.Value < 0.05 & logFC > 0.3  ~ "Up",
      P.Value < 0.05 & logFC < -0.3 ~ "Down",
      TRUE ~ "Not Significant"
    )
  )

DFexpression.nSolver$significance <- significance$significance
# Pick top 10 most significant genes (Up or Down)
significance <- significance$significance

top_labels <- filter_all(DFexpression.nSolver$significance != "Not Significant")
  slice_head(order_by = P.Value, n = 10)
top_labels$Name <- rownames(top_labels)

ggplot(DFexpression.nSolver, aes(x = logFC, y = P.Value)) +
  scale_y_reverse() +
  geom_point(alpha = 0.7) +
 # ylim(-1.25,0.25)+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot with Top Genes",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~P~value)
  ) +
  theme_void()

# Plot
ggplot(DFexpression.nSolver, aes(x = logFC, y = P.Value, color = significance)) +
  scale_y_reverse() +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Up" = "#5830F7", "Down" = "#FFA46B", "Not Significant" = "gray")) +
 # geom_text(data = top_labels, aes(label = Name), vjust = -0.5, size = 3, show.legend = FALSE) +
  labs(
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~P~value),
    color = "Significance"
  ) +
  theme_minimal()


# Extract phenotype data
pdata <- pData(split_genes)

# Mutate the group column
pdata <- pdata %>%
  mutate(
    groups = case_when(
      DonorAge <= 30 ~ "Young",
      DonorAge >= 40 ~ "Old",
      TRUE ~ NA_character_
    )
  )

# Assign back to ExpressionSet
pData(split_genes) <- pdata


split_genes <- split_genes[, !sampleNames(split_genes) %in% "20211103_AN Set9_D066_09.RCC"]



# Extract phenotype data
pdata <- pData(split_genes)

# Identify valid sample names by DonorAge condition
valid_samples <- rownames(pdata)[pdata$DonorAge <= 30 | pdata$DonorAge >= 40]

# Subset the ExpressionSet using those sample names
split_genes <- split_genes[, valid_samples]

