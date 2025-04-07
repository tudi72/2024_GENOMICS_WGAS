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
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data = ds$snps)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = "chr3:1-30000000")
kp <- plotKaryotype(plot.type = 4)
kp <- kpPlotManhattan(kp, data=ds$snps, highlight = ds$peaks, points.cex = 0.8)

kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "brewer.set1", r0=autotrack(1,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "2blues", r0=autotrack(2,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "greengray", r0=autotrack(3,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = "rainbow", r0=autotrack(4,5))
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = c("orchid", "gold", "orange"), r0=autotrack(5,5))

#-----------------Creating a gradient for point visualization-------------------
transf.pval <- -log10(ds$snps$pval)
points.col <- colByValue(transf.pval, colors=c("#BBBBBB", "orange"))
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, points.col = points.col)

#-----------------Modifying line and highlight colors---------------------------
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps,
                      highlight = ds$peaks, highlight.col = "orchid",
                      suggestive.col="orange", suggestive.lwd = 3,
                      genomewide.col = "red", genomewide.lwd = 6)

#-----------------------Controlling axis look-----------------------------------
transf.pval <- -log10(ds$snps$pval)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, pval = transf.pval, logp = FALSE )
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)
kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps)
ymax <- kp$latest.plot$computed.values$ymax
ticks <- c(0, seq_len(floor(ymax)))
kpAxis(kp, ymin=0, ymax=ymax, tick.pos = ticks)
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)

#---------------------Labelling top significant SNP per chromosome--------------
kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)


snps <- kp$latest.plot$computed.values$data
suggestive.thr <- kp$latest.plot$computed.values$suggestiveline

#---------------Get the names of the top SNP per chr----------------------------
top.snps <- tapply(seq_along(snps), seqnames(snps), function(x) {
 in.chr <- snps[x]
 top.snp <- in.chr[which.max(in.chr$y)]
 return(names(top.snp))
})
#--------------------Filter by suggestive line----------------------------------
top.snps <- top.snps[snps[top.snps]$y>suggestive.thr]
top.snps <- snps[top.snps]

kp <- plotKaryotype(plot.type=4)
kp <- kpPlotManhattan(kp, data=ds$snps, ymax=10)
kpAxis(kp, ymin = 0, ymax=10, numticks = 11)
kpText(kp, data = top.snps, labels = names(top.snps), ymax=10, pos=3)

#-------------------SIGNIFICANT VARIANTS INFO-----------------------------------
ensembl <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

nt.biomart <- getBM(
  attributes = c("refsnp_id","allele","chr_name","chrom_start","chrom_strand","associated_gene","ensembl_gene_stable_id"),
  values     = names(top.snps),
  mart       = ensembl,
  verbose    = TRUE,
  uniqueRows = TRUE,
  filters    = "snp_filter"
  )

#---------TASK2. PANCREATIC CANCER RISK GENOTYPE--------------------------------
#Generate a processed data set based on risk genes for pancreatic cancer
genes_risk <- c("PDX1" ,        "KLF5",   "BCAR1", 
                "HNF1B",   "LINC00673",     "GRP", 
                "WNT2B",       "NOC2L",   "ZNRF3",
                "ETAA1",        "TP63",    "TERT", 
                "TNS3" ,       "SUGCT",   "HNF4G", 
                "PVT1")
genes_risk_rownames <- do.call(paste, c(expand.grid(genes_risk, c("low", "intermediate", "high")), sep = "_"))

prognosis_table <- matrix(nrow = length(genes_risk_rownames))
outcome_vector <- c()

set.seed(322)
# iterate each patient 
for (patient in 1:1000) {
  # create a seed for display high or low
  outcome_seed <- sample(1:100, 1)
  if (outcome_seed <= 80) {
    gene_outcome <-  c()
    for (gene in genes_risk) {
      total_snps       <- sample(10:30, 1)
      snp_entry        <- table(sample(x = c("low", "intermediate", "high"), size = total_snps,  replace = TRUE, prob = c(0.6, 0.3, 0.1)))
      names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
      gene_outcome     <- c(gene_outcome, snp_entry)
    }
    outcome_vector  <- c(outcome_vector, "low")
    prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])

    
  } else if (outcome_seed >= 95) {
    gene_outcome <-  c()
    for (gene in genes_risk) {
      total_snps <- sample(10:30, 1)
      snp_entry <- table(sample(x = c("low", "intermediate", "high"), size = total_snps,  replace = TRUE, prob = c(0.2, 0.3, 0.5)))
      names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
      gene_outcome <- c(gene_outcome, snp_entry)
    }
    outcome_vector <- c(outcome_vector, "high")
    prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])
    
  } else {
    gene_outcome <-  c()
    for (gene in genes_risk) {
      total_snps <- sample(10:30, 1)
      snp_entry <- table(sample(x = c("low", "intermediate", "high"), size = total_snps,  replace = TRUE, prob = c(0.5, 0.3, 0.2)))
      names(snp_entry) <- paste(gene, names(snp_entry), sep = "_")
      gene_outcome <- c(gene_outcome, snp_entry)
    }
    outcome_vector <- c(outcome_vector, "intermediate")
    prognosis_table <- cbind(prognosis_table, patient = gene_outcome[match(genes_risk_rownames, names(gene_outcome))])
  }
}

#Clean up the results data set
prognosis_table                         <- prognosis_table[, -1]
rownames(prognosis_table)               <- genes_risk_rownames
colnames(prognosis_table)               <- names(outcome_vector) <- paste0("patient", 1:1000)
prognosis_table[is.na(prognosis_table)] <- 0

prognosis_table[1:5, 1:5]
outcome_vector[1:5]

#For added "realism" we can include missing data
missing_data <- sample(1:1000, 250)

for (patient in missing_data) {
  sim_missing                <- prognosis_table[, patient]
  num_missing                <- sample(genes_risk, sample(1:7, 1))
  sim_missing[grepl(paste(num_missing,collapse="|"), names(sim_missing))] <- NA
  prognosis_table[, patient] <- sim_missing
}

shape2 <- dim(prognosis_table)
cat("Shape of the matrix:", shape2[1], "x", shape2[2], "\n")

# ----------------- MISSING DATA-----------------------------------------------
num_missing        <- sum(is.na(prognosis_table))
missing_percentage <- sum( num_missing / (nrow(prognosis_table) * ncol(prognosis_table))) * 100
# there are 6.18% missing values 

# ----------------- IMPUTATION ------------------------------------------------
library(parallel)
library(doParallel)


cl <- makeCluster(4)
registerDoParallel(cl)

# (feature x patients) -> (patients x features)
prognosis_table_transposed <- t(prognosis_table)

# Predictive Mean ~ long time 
imputed_data <- mice(prognosis_table_transposed, m = 1, method = 'pmm', seed = 123,maxit=1,parallel=TRUE)

stopCluster(cl)

prognosis_table  <- complete(imputed_data)

#----------------DATASET--------------------------------------------------------
target_gene  <- "TERT_high" 
y            <- prognosis_table[[target_gene]]
y            <- make.names(y)

X            <- prognosis_table[, !(colnames(prognosis_table) %in% target_gene), drop = FALSE]

# -------------TRAIN_TEST_SPLIT ------------------------------------------------
library(ranger)
library(caret)
trainIndex    <- createDataPartition(prognosis_table[, ncol(prognosis_table)], p = .8, 
                                  list = FALSE, 
                                  times = 1)
train_data    <- prognosis_table[trainIndex, ]
test_data     <- prognosis_table[-trainIndex, ]

#Specify training parameters
train_control <- trainControl(method = "cv", number = 10,classProbs=TRUE)

if (any(is.na(train_data))) {
  train_data <- na.omit(train_data)
  warning("NA values detected and removed.")
}

#----------------PANCREATIC CANCER RISK GENE PREDICTION-------------------------

#Train a selected model with specified parameters
model <- train(y           = y_train,  # Target (now a factor)
               x           = X_train, # Assuming last column is the target
               method      = "ranger",
               metric      = "ROC",
               trControl   = train_control,
               tuneLength  = 5)

#Apply it to the testing dataset
y_predictions <- predict(model, newdata = as.data.frame(test_data))
y_predictions <- factor(y_predictions)
y_true        <- as.factor(test_data[[target_gene]])
all_levels    <- union(levels(y_true), levels(y_predictions))

y_predictions <- factor(y_predictions, levels = all_levels)
y_true        <- factor(y_true, levels = all_levels)
#------------------CONFUSION MATRIX---------------------------------------------
confusion     <- confusionMatrix(y_predictions,y_true)
conf_table    <- as.table(confusion)
conf_df       <- as.data.frame(conf_table)

head(conf_df)

cat("\n[OVERALL ACCURACY]: ", confusion$overall['Accuracy'], "\n")


cat("\n[PRECISION FOR EACH CLASS]:\n")
precision    <- confusion$byClass[, "Precision"]
print(precision)

