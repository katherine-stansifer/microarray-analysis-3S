library(oligo)
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)
library(affy)
library(mgu74av2.db)
library(mgu74bv2.db)
library(mgu74cv2.db)
library(mouse4302.db)


# Process GSE4408(synch data) 
starting_dir <- getwd()
cel.files <- list.celfiles('data/GSE54408_RAW/total/', full.names=TRUE)
GSE54408_raw <- read.celfiles(cel.files)
GSE54408_rma <- oligo::rma(GSE54408_raw, target = "core")
pnames<-featureNames(GSE54408_rma)
en_GSE54408 <- unlist(mget(pnames, mogene10sttranscriptclusterENTREZID, ifnotfound=NA))


# Process GSE12769 (one of the unsynch datasets)
cel.files <- list.celfiles('data/GSE12769_raw', full.names=TRUE)
GSE12769_raw <- read.celfiles(cel.files)
GSE12769_pma <- paCalls(GSE12769_raw)$calls
GSE12769_rma <- oligo::rma(GSE12769_raw)
en_GSE12769 <- unlist(mget(featureNames(GSE12769_rma),envir=get("mouse4302ENTREZID"),ifnotfound=NA))

# Process GSE926 (other unsynch dataset)
GSE926_raw_a <- read.celfiles(list.celfiles("data/GSE926_RAW/MGU74v2_A/", full.names=TRUE))
GSE926_pma_a <- paCalls(GSE926_raw_a)
GSE926_rma_a <- oligo::rma(GSE926_raw_a)
GSE926_raw_b <- read.celfiles(list.celfiles("data/GSE926_RAW/MGU74v2_B/", full.names=TRUE))
GSE926_pma_b <- paCalls(GSE926_raw_b)
GSE926_rma_b <- oligo::rma(GSE926_raw_b)
GSE926_raw_c <- read.celfiles(list.celfiles("data/GSE926_RAW/MGU74v2_C/", full.names=TRUE))
GSE926_pma_c <- paCalls(GSE926_raw_c)
GSE926_rma_c <- oligo::rma(GSE926_raw_c)

GSE926_rma_rough <- rbind(exprs(GSE926_rma_a),exprs(GSE926_rma_b),exprs(GSE926_rma_c))
GSE926_rma_all <-new("ExpressionSet")
exprs(GSE926_rma_all) <- GSE926_rma_rough
GSE926_pma_all <- rbind(GSE926_pma_a$calls,GSE926_pma_b$calls,GSE926_pma_c$calls)
entrezA<-mget(featureNames(GSE926_rma_a),envir=get("mgu74av2ENTREZID"),ifnotfound=NA)
entrezB<-mget(featureNames(GSE926_rma_b),envir=get("mgu74bv2ENTREZID"),ifnotfound=NA)
entrezC<-mget(featureNames(GSE926_rma_c),envir=get("mgu74cv2ENTREZID"),ifnotfound=NA)
en_GSE926 <- unlist(c(entrezA,entrezB,entrezC))

# Save work
save.image("results/basic_proc.RData")
