library(TCGAbiolinks)
library(dplyr)
library(edgeR)
library(ggplot2)
library(factoextra)
library(gplots)

# Download the expression data: SKCM
query.exp <- GDCquery(project = "TCGA-SKCM",
                      data.category = "Transcriptome Profiling",
                      data.type ="Gene Expression Quantification" ,
                      workflow.type = 'HTSeq - Counts',
                  )

GDCdownload(query.exp)

data <- GDCprepare(query.exp, save = TRUE, 
                   save.filename = "ProteinExpression.rda",
                   remove.files.prepared = TRUE)


countdata <- TCGAanalyze_Preprocessing(data) %>%
  as.data.frame()


##### CLEANING DATA ### 

myCPM <- cpm(countdata)
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
clean_counts <- countdata[keep,]


rownames(clean_counts) <- data@colData@listData$patient


#plot
plot(myCPM[,1],HTSeq_Counts[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5)

#### SAVING DATA ######
genesNames <- as.data.frame(data@rowRanges@elementMetadata@listData[["external_gene_name"]])
rownames(genesNames)<- data@rowRanges@elementMetadata@listData[["ensembl_gene_id"]]
names(genesNames) <- "gene_name"


write.csv(clean_counts,'/home/sergio/Projects/IVJornadas.Bioinformatica/melanoma-datathon/melanoma.datathon/clean.counts.csv')
write.csv(genesNames,'/home/sergio/Projects/IVJornadas.Bioinformatica/melanoma-datathon/melanoma.datathon/genesNames.csv')
write.csv(countdata,'/home/sergio/Projects/IVJornadas.Bioinformatica/melanoma-datathon/melanoma.datathon/countdata.csv')
info_df <- as.data.frame(data@colData)
newdf <- na.omit(df_tsne_model)

###################################################################################################################################

####### T-SNE #############

library(Rtsne)

tsne_model <- Rtsne(clean_counts)

df_tsne_model <- as.data.frame(tsne_model$Y)
rownames(df_tsne_model) <- data@colData@rownames

ggplot(df_tsne_model, aes(x=V1, y=V2)) +  
  geom_point(size=1) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

# determine number of clusters: TODOS LOS METODOS NOS DICEN QUE EL NUMERO ÓPTIMO DE CLUSTERS ES DE 3 (kmeans)
fviz_nbclust(df_tsne_model, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df_tsne_model, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
set.seed(123)
fviz_nbclust(df_tsne_model, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

###### kmeans clustering on tsne reduction ###########
fit_cluster_kmeans <- kmeans(scale(df_tsne_model), 3)

ggplot(df_tsne_model, aes(x=V1, y=V2)) +  
  geom_point(size=1, col = fit_cluster_kmeans$cluster) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")


###### hieratical clustering with heatmaps ###########

logcounts <- cpm(clean_counts,log=TRUE)
var_genes <- apply(logcounts, 1, var)

# Get the gene names for the top 500 most variable genes
selected_patients <- names(sort(var_genes, decreasing=TRUE))[1:1500]

# Subset logcounts matrix
highly_variable_patients<- logcounts[selected_patients,] 


#heatmap clustering
library(pheatmap)
out2 <- pheatmap(highly_variable_patients, scale = "row",
                 show_rownames = F, show_colnames = F,
                 color = greenred(5))

results <- cutree(out2$tree_col, 3)


View(results)


#Select genes "expresion" from clustered patients
results_df <- as.data.frame(results)
names(results_df) <- "hieratical_cluster_patient"
results_df$k_cluster_patient<- fit_cluster_kmeans$cluster

library(plyr)
results_df$hieratical_cluster_patient<- mapvalues(results_df$hieratical_cluster_patient, from = c(1,2,3), to = c("Group1", "Group2", "Group3"))
results_df$k_cluster_patient<- mapvalues(results_df$k_cluster_patient, from = c(1,2,3), to = c("Group1", "Group2", "Group3"))

write.csv(results_df,'/home/sergio/Projects/IVJornadas.Bioinformatica/melanoma-datathon/melanoma.datathon/results_cluster.csv')

#estudio de expresion diferencial entre los clusters del heatmap

labels <- paste(results_df$hieratical_cluster_patient)
group <- paste(results_df$hieratical_cluster_patient,sep=".")
group <- factor(group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# voom normaliza los conteos con log(CPMs)
v <- voom(clean_counts,design,plot = TRUE) 

fit <- lmFit(v)
names(fit)

cont.matrix <- makeContrasts(results_df$hieratical_cluster_patient,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
topgenes_hieratical <- topTable(fit.cont,coef=1,sort.by="p")

common <- intersect(rownames(topgenes_hieratical), rownames(genesNames))  
selected_names <- genesNames[common,] %>%
  as.data.frame()

topgenes_hieratical <- cbind(topgenes_hieratical, selected_names)
topgenes_hieratical$method <- "Hierarchical clustering"



#ES DECIR, independientemente del cluster obtenemos los mismos gnes resultado en cuanto a la importancia genética ¿cambia la expresion segun los clusters?
common_df <- rbind(topgenes_clusters_kmeans,topgenes_hieratical)
ggplot(common_df,aes(y= AveExpr, x= ., fill = method)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_light(base_size=10)


####### CLINICAL DATA FROM ARTICLE ######## 

NIHMS698912_supplement_3 <- read_excel("Projects/IVJornadas.Bioinformatica/melanoma-datathon/melanoma.datathon/NIHMS698912-supplement-3.xlsx", 
                                       +     sheet = "Supplemental Table S1D")

clinical_data <- as.data.frame(NIHMS698912_supplement_3)
rownames(clinical_data) <- clinical_data$Name
clinical_data$Name <- NULL

rownames(results_df) <- substr(rownames(results_df), 1, 15)


common_patients <- intersect(rownames(results_df), rownames(clinical_data))  

selected_patients <- clinical_data[common_patients,] %>%
  as.data.frame()

selected_patients1 <- results_df[common_patients,] %>%
  as.data.frame()

results_df_clinical <- cbind(selected_patients1, selected_patients)

write.csv(results_df_clinical, "/home/sergio/Projects/IVJornadas.Bioinformatica/melanoma-datathon/melanoma.datathon/results_df_clinical.csv")


# Plot final map the heatmap
## data prep ## 
highly_variable_patients2 <- highly_variable_patients %>%
  as.data.frame()
colnames(highly_variable_patients2)<- substr(colnames(highly_variable_patients2),1,15)

highly_variable_patients2 <- highly_variable_patients2[,common_patients] %>%
  as.matrix()
  
## complex heatmap ## https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.html#toc_0
library(ComplexHeatmap)
plotha <- results_df_clinical[,1:5]
plotha$ALL_SAMPLES <- NULL
plotha$Tissue_origin <- results_df_clinical$REGIONAL_VS_PRIMARY


ha = HeatmapAnnotation(df = plotha,col = list(hieratical_cluster_patient = c("Group1" = "orange", "Group2" = "purple", "Group3" = "red"),
                                              k_cluster_patient = c("Group1" = "blue", "Group2" = "green", "Group3" = "darkorange"),
                                              MUTATIONSUBTYPES = c("-" = "white", "Triple_WT"= "grey", "RAS_Hotspot_Mutants" = "dodgerblue1", "NF1_Any_Mutants" = "darkgreen", "BRAF_Hotspot_Mutants" = "gold")))

Heatmap(highly_variable_patients2, name = "Expression level", show_column_names = F, show_row_names = F, 
        column_title = "Patients", row_title = "Genes",
        col=rev(morecols(1000)),
        top_annotation = ha)

