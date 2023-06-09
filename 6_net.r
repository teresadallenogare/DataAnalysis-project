# ============================ NETWORK BASED ANALYSIS ==============================
# Author: Teresa Dalle Nogare
# Version : 19-04-2023
# ORA : DAVID tool
# GSEA : GSEA tool
# NBA : expand list of genes including genes functionally related with original ones -> through non interaction of list
# then, do functional analysis with extended list
# STRING : for each gene find neighbors that are genes which are functionally related. Need a network of interactions 
# to find functionally related genes
# given a gene interaction network, clusters correspond to genes that are functionally related
# XD score : estimate of how close the gene set is close to the original list, both on distance and number of 
# connections, without accounting for other features
# p-value : to capture the statistical significance
# try to figure out the function both considering the gene ontology bout also the interactions - collections of interactions
# of genes that are not in the list but that are neignours

# STRING : first shell are the neighbours
# ENRICHNET : XD-score anything below 1 is not interesting
# regression plot : x : Fisher t-test , y :XD- significance of overlap 
# PATHFINDR
# ==============================================================================
# ====================== Install packages  ===================================== 
install.packages("BiocManager")
install.packages("igraph") 
install.packages("ggplot2")
install.packages('caret')
install.packages("glmnet") 
install.packages('Matrix')
install.packages('igraph')
BiocManager::install("GEOquery")
BiocManager::install("useful")
BiocManager::install("genefilter")
BiocManager::install("rScudo")
BiocManager::install("enrichplot")

library("GEOquery")
library("useful")
library(ggplot2)
library("RColorBrewer")
library("genefilter")
library("MASS")
library("gplots")
library('ggplot2')
library('lattice')
library("caret")
library("glmnet")
library("rScudo")
library("enrichplot")
library("XML")
library(annotate)

BiocManager::install("KEGGREST") 
BiocManager::install("KEGGgraph") 
BiocManager::install("AnnotationDbi") 
BiocManager::install("org.Hs.eg.db") 
BiocManager::install("HsAgilentDesign026652.db")
library("KEGGREST") 
library("KEGGgraph") 
library("AnnotationDbi") 
library("org.Hs.eg.db")
library('HsAgilentDesign026652.db') # agilent conversion id probes
install.packages("pathfindR")
library("pathfindR")
library('limma')
library("randomForest") 
library('DESeq2')




load(file = "~/Desktop/DA-project/0_data/Y.Rdata")
# order dataset based on treatment type and group
Y <- Y[, c(rep(1:3), rep(10:12), rep(19:21), rep(28:30), rep(4:9), rep(13:18), rep(22:27), rep(31:35))]
group <- c(rep('FBS',12), rep('GCRACOMB',23))
f <- factor(group)
# determination of logFC
design <- model.matrix(~ f+0 , data.frame(Y))

fit <- lmFit(Y, design)
cts <- paste('fFBS', 'fGCRACOMB', sep="-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", number = length(Y))

# redefinition of group
group <- c(rep('FBS',12), rep('GCRA-COMB',23))
f <- factor(group)

# Random forest : selection of 100 most relevant genes
set.seed(1234) # deterministic algorithm
rf <- randomForest(x=t(Y), y=as.factor(group), ntree=1000)
probe.names  <- rownames(rf$importance)
top100index  <- order(rf$importance,decreasing = TRUE)[1:100]
top100       <- probe.names[order(rf$importance,decreasing = TRUE)[1:100]]
top100df     <- data.frame(top100)
# select only the 100 genes in top100
tT <- tT[row.names(tT) %in% top100, ]
# sort data in the same order of ttY
tT <- tT[match(top100, rownames(tT)), ]

gene.entrez  <- lookUp(top100, "HsAgilentDesign026652", "ENTREZID")
gene.ensembl <- lookUp(top100, "HsAgilentDesign026652", "ENSEMBL")
gene.symbols <- lookUp(top100, "HsAgilentDesign026652", "SYMBOL")

top100df["ENZYME"]  <- list(as.character(gene.entrez))
top100df["SYMBOL"]  <- list(as.character(gene.symbols))
top100df["ENSEMBL"] <- list(as.character(gene.ensembl))

# row ttest
set.seed(1234)
ttY <- rowttests(Y,f)

top100df["p.value"] <- ttY$p.value[top100index]
top100df["logFC"]   <- tT$logFC

#top100df_entrez <- top100df[c(2,5)]

top100df.pathfindR  <- top100df[,c("SYMBOL","logFC","p.value")]
PF_input <- top100df.pathfindR
input_testing(PF_input, p_val_threshold = 0.05)

## -------------------------- Enrichment analysis -----------------------------
PF_output <- run_pathfindR(PF_input, 
                           gene_sets = "GO-BP", # GO-ALL, KEGG
                           pin_name_path = "STRING", # Biogrid, STRING, IntAct
                           iterations = 1,
                          # output_dir = '~/Desktop/DA-project/gsea_data/FBS-GCRA/pathfindR/',
                           list_active_snw_genes = TRUE,
                          plot_enrichment_chart = TRUE)  # report non significant active subnetwork genes
# returns a dotplot : displays the most signficant enriched terms

# ------------------------ Clustering of enriched terms ------------------------
#par(mfcol=c(1.5,3.5))
#PF_clustering <- cluster_enriched_terms(PF_output, 
#                                        plot_hmap = TRUE, # heatmap of 
#                                        plot_dend = TRUE,
#                                        clu_method = 'average')      
# ------------------------ Visualization enrich results  -----------------------
#input_processed <- input_processing(PF_input)
#visualize_terms(result_df = PF_output,
#                input_processed = input_processed,
#                hsa_KEGG = TRUE)
# pathways saved in folder : DA_project/term_visualizations  
#visualize_terms(result_df = PF_output,
#                input_processed = input_processed,
#                hsa_KEGG = FALSE,
#                pin_name_path = 'Biogrid')
# ------------------------ Term genes  -----------------------------------------
term_gene_heatmap(result_df = PF_output, PF_input)
term_gene_graph(result_df = PF_output,
                use_description = FALSE)

UpSet_plot(PF_output, PF_input,
           method = 'heatmap',
           num_terms = 100
             )




# sort top100df based on pvalue
top100df_sort <- top100df[order(top100df$p.value,decreasing = FALSE),]


library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de)

