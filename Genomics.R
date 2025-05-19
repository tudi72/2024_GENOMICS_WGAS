 # if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")


# BiocManager::install("pbapply")
 # BiocManager::install("xfun")
 # BiocManager::install("regioneR")
 # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19.masked")
 # BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
 # BiocManager::install("karyoploteR",force=TRUE)
 # BiocManager::install("biomaRt")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("vcfR")

# install.packages("httpgd")
# install.packages("languageserver")

library("mice")
library("httr")
library("jsonlite")
library("dplyr")
library("vcfR")
library("pbapply")
library("karyoploteR")
library("BSgenome.Hsapiens.UCSC.hg19.masked")
library("BSgenome.Hsapiens.UCSC.hg38")
library("regioneR")
library("biomaRt")
library("caret")



#----------------TASK 1.WHOLE-GENOME DATASET -----------------------------------
set.seed(322)
set_config(timeout(360000))  # Set timeout to 600 seconds

createDataset <- function(num.snps = 20000, max.peaks = 5) {
 hg19.genome  <- filterChromosomes(getGenome("hg19")) 
 snps         <- sort(createRandomRegions(nregions    = num.snps, 
                                          length.mean = 1, 
                                          length.sd   = 0, 
                                          genome       = filterChromosomes(getGenome("hg19"))))
 names(snps)              <- paste0("rs", seq_len(num.snps))
 snps$pval                <- rnorm(n = num.snps, mean = 0.5, sd = 1)
 snps$pval[snps$pval < 0] <- -1*snps$pval[snps$pval < 0]
 #---------------- SIGNIFICANT PEAKS -------------------------------------------

 peaks <- createRandomRegions(runif(1, 1, max.peaks), 8e6, 4e6)

 for(npeak in seq_along(peaks)) {
  snps.in.peak            <- which(overlapsAny(snps, peaks[npeak]))
  snps$pval[snps.in.peak] <- runif(n = length(snps.in.peak), 
                                  min = 0.1, 
                                  max = runif(1, 6, 8))
 }
 snps$pval <- 10^(-1*snps$pval)
 return(list(peaks = peaks, snps = snps))
}


ds <- createDataset()
ds$snps

#-----------------------Visualizing genotyping data-----------------------------
png("plot1_basic_manhattan.png", width=1200, height=600)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data = ds$snps)
dev.off()
png("plot2_highlight_chr3.png", width=1200, height=600)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = "chr3:1-30000000")
dev.off()
png("plot3_highlight_peaks.png", width=1200, height=600)
kp <- plotKaryotype(plot.type = 4)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = ds$peaks, points.cex = 0.8)
dev.off()
png("plot4_multicolor_tracks.png", width=1200, height=800)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set1", r0=autotrack(1,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "2blues", r0=autotrack(2,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "greengray", r0=autotrack(3,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "rainbow", r0=autotrack(4,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = c("orchid", "gold", "orange"), r0=autotrack(5,5))
dev.off()
#-----------------Creating a gradient for point visualization-------------------
transf.pval <- -log10(ds$snps$pval)
points.col  <- colByValue(transf.pval, colors=c("#BBBBBB", "orange"))
kp          <- plotKaryotype(plot.type=4)
kp          <- kpPlotManhattan(kp, data=ds$snps, points.col = points.col)

#-----------------Modifying line and highlight colors---------------------------
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, 
  data            = ds$snps,
  highlight       = ds$peaks, 
  highlight.col   = "orchid",
  suggestive.col  ="orange", 
  suggestive.lwd  = 3,
  genomewide.col  = "red", 
  genomewide.lwd  = 6)
#-----------------------Controlling axis look-----------------------------------
library(karyoploteR)

transf.pval <- -log10(ds$snps$pval) # Manhattan -log10(pval)
png("plot5_log10pval.png", width=1200, height=600)
kp          <- plotKaryotype(plot.type=4)
kp          <- kpPlotManhattan(kp, data=ds$snps, pval = transf.pval, logp = FALSE)
dev.off()

png("plot6_axis_autoymax.png", width=1200, height=600) # Manhattan with axis
kp        <- plotKaryotype(plot.type=4)
kp        <- kpPlotManhattan(kp, data=ds$snps)
kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
dev.off()

png("plot7_custom_ticks.png", width=1200, height=600)
kp      <- plotKaryotype(plot.type=4)
kp      <- kpPlotManhattan(kp, data=ds$snps)
ymax    <- kp$latest.plot$computed.values$ymax
ticks   <- c(0, seq_len(floor(ymax)))
kpAxis(kp, ymin=0, ymax=ymax, tick.pos = ticks)
dev.off()

png("plot8_ymax_10_ticks.png", width=1200, height=600)
kp      <- plotKaryotype(plot.type=4)
kp      <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)
dev.off()

png("plot9_manhattan_ymax10.png", width=1200, height=600)
kp      <- plotKaryotype(plot.type=4)
kp      <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)
dev.off()
#---------------------Labelling top significant SNP per chromosome--------------
snps            <- kp$latest.plot$computed.values$data
suggestive.thr  <- kp$latest.plot$computed.values$suggestiveline

#---------------Get the names of the top SNP per chr----------------------------
top.snps  <- tapply(seq_along(snps), seqnames(snps), function(x) {
 in.chr   <- snps[x]
 top.snp  <- in.chr[which.max(in.chr$y)]
 return(names(top.snp))
})
#--------------------Filter by suggestive line----------------------------------
top.snps <- top.snps[snps[top.snps]$y>suggestive.thr]
top.snps <- snps[top.snps]

png("plot10_manhattan_plot_peak.png", width = 1200, height = 800)  # You can adjust width/height
kp      <- plotKaryotype(plot.type=4)
kp      <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)
kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3)
dev.off()
#-------------------SIGNIFICANT VARIANTS INFO-----------------------------------
ensembl   <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

nt.biomart <- getBM(
  attributes = c("refsnp_id","allele","chr_name","chrom_start","chrom_strand","associated_gene","ensembl_gene_stable_id"),
  values     = names(top.snps),
  mart       = ensembl,
  verbose    = TRUE,
  uniqueRows = TRUE,
  filters    = "snp_filter"
  )
# ---------TASK2. PANCREATIC CANCER RISK GENOTYPE------------------------------
genes_risk <- c("PDX1" , "KLF5"     , "BCAR1", 
                "HNF1B", "LINC00673", "GRP", 
                "WNT2B", "NOC2L"    , "ZNRF3",
                "ETAA1", "TP63"     , "TERT", 
                "TNS3" , "SUGCT"    , "HNF4G", 
                "PVT1")

genes_risk_rownames <- do.call(paste, c(expand.grid(genes_risk, c("low", "intermediate", "high")), sep = "_"))
prognosis_table     <- matrix(nrow = length(genes_risk_rownames))
outcome_vector      <- c()

set.seed(322)
for (patient in 1:1000) {
  outcome_seed <- sample(1:100, 1)
  
  if (outcome_seed <= 80) {
    prob <- c(0.6, 0.3, 0.1)
    outcome_vector <- c(outcome_vector, "low")
  } else if (outcome_seed >= 95) {
    prob <- c(0.2, 0.3, 0.5)
    outcome_vector <- c(outcome_vector, "high")
  } else {
    prob <- c(0.5, 0.3, 0.2)
    outcome_vector <- c(outcome_vector, "intermediate")
  }

  gene_outcome <- c()
  for (gene in genes_risk) {
    total_snps       <- sample(10:30, 1)
    snp_entry        <- table(sample(x = c("low", "intermediate", "high"), size = total_snps, replace = TRUE, prob = prob))
    names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
    gene_outcome     <- c(gene_outcome, snp_entry)
  }
  
  prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])
}

prognosis_table <- prognosis_table[, -1]
rownames(prognosis_table) <- genes_risk_rownames
colnames(prognosis_table) <- paste0("patient", 1:1000)
prognosis_table[is.na(prognosis_table)] <- 0

# Introduce missing data
missing_data <- sample(1:1000, 250)
for (patient in missing_data) {
  sim_missing <- prognosis_table[, patient]
  num_missing <- sample(genes_risk, sample(1:7, 1))
  sim_missing[grepl(paste(num_missing, collapse = "|"), names(sim_missing))] <- NA
  prognosis_table[, patient] <- sim_missing
}

cat("Shape of the matrix:", dim(prognosis_table)[1], "x", dim(prognosis_table)[2], "\n")

# ----------------- IMPUTATION ------------------------------------------------
library(parallel)
library(doParallel)
library(mice)

cl <- makeCluster(4)
registerDoParallel(cl)

prognosis_table_transposed <- t(prognosis_table)

imputed_data <- mice(prognosis_table_transposed, m = 1, method = 'pmm', seed = 123, maxit = 1, parallel = TRUE)

stopCluster(cl)

# Get complete imputed data (patients x features)
imputed_patients <- complete(imputed_data)

# ------------------ K-MEANS CLUSTERING -----------------------------
# Optional scaling
imputed_scaled <- scale(imputed_patients)

# Perform K-means clustering
set.seed(123)
k_clusters          <- 3
kmeans_result       <- kmeans(imputed_scaled, centers = k_clusters, nstart = 25)
cluster_assignments <- kmeans_result$cluster
cat("Cluster sizes:\n")
print(table(cluster_assignments))

# ------------------ PCA VISUALIZATION --------------------
library(ggplot2)

pca_result         <- prcomp(imputed_scaled)
pca_df             <- as.data.frame(pca_result$x[, 1:2])
pca_df$Cluster     <- factor(cluster_assignments)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "K-Means Clustering of Pancreatic Cancer Genotype Profiles",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()
