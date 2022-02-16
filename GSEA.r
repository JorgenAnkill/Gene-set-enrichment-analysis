######################################################################################################################################################################################################################################################################
# -- GSEA using the genesets from MSigDB
######################################################################################################################################################################################################################################################################

GSEA <- function(genes,geneset.collection,geneset_subset,pcut,geneset.threshold,version,minGenes){
if(missing(minGenes)){minGenes <- 0} 
if(missing(pcut)){pcut <- 0.05}
if(missing(geneset.threshold)){geneset.threshold <- 0}
if(missing(geneset.collection)){geneset.collection <- c("H","C1","C2","C4","C5","C6","C7","C8")}
if(missing(version)){version <- "v7.5.1"}

# Load packages and install if missing
my_packages <- c("fst")                                   
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed) #Install necessary packages if missing

library(fst)
if(version=="v7.5.1"){genesets <- read_fst("/open/work/Jorgen/Data/GSEA/Molecular signatures database/v7.5.1/msigdb_v7.5.1_GMTs/MSigDB_GSVA_data.fst")}

# -- With unspecified background
total.genes <- length(unique(genesets$Gene))
genes.length <- sum(unique(genesets$Gene)%in%genes)
genesets <- genesets[genesets$GeneSetCollection%in%geneset.collection,]
if(!missing(geneset_subset)){genesets <- genesets[genesets$GeneSetName%in%genesets[genesets$GeneSetName%in%unique(grep(paste(geneset_subset,collapse="|"),genesets$GeneSetName,value=TRUE)),]$GeneSetName,]} #HP, GOBP, GOMF, GOCC

res <- data.frame(matrix(nrow=length(unique(genesets$GeneSetName)),ncol=5))
colnames(res) <- c("GeneSetName","GenesInGeneset","GenesInOverlap","p.value","FDR.q.value")
res$GeneSetName <- unique(genesets$GeneSetName)

dup <- genesets[!duplicated(genesets$GeneSetName),]
res$GenesInGeneset <- dup$GenesInGeneset

out <- split(genesets,f=genesets$GeneSetName); out <- out[match(res$GeneSetName,names(out))]
vec <- vector()
for(i in 1:nrow(res)){vec[i] <- length(intersect(out[[i]]$Gene,genes))}

res$GenesInOverlap <- vec

pvals <- vector()
fe <- vector()
N <- total.genes
k <- genes.length
for(i in 1:nrow(res)){
q <- res$GenesInOverlap[i]
m <- res$GenesInGeneset[i]
n <- N-m
pvals[i] <- phyper(q-1,m,n,k,lower.tail=FALSE)
fe[i] <- q*N/(m*k)}

res$p.value <- pvals
res$FE <- fe
res$FDR.q.value <- p.adjust(pvals,method="fdr")
res <- res[order(res$p.value,decreasing=FALSE),]
res <- res[res$GenesInGeneset>=geneset.threshold,]
result <- res[res$FDR.q.value<=pcut & res$GenesInOverlap>=minGenes,]
return(result)
print(noquote("GSEA completed"))
}