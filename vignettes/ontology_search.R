# !diagnostics off
library(qusage)
###### Load data

brca <- readRDS("/home/mtr/rfactory/brca_data.RDS")
brca <- readRDS("C:/RFactory/Rdata_brca/brca_data.RDS")
cvfit <- readRDS("C:/RFactory/predsurv/performance_results/ridge2_model_brca_pooled.RDS")

hallmark_all <- qusage::read.gmt("/home/mtr/rfactory/predsurv/vignettes/Gene sets/h.all.v6.1.symbols.gmt")

hallmark_all <- qusage::read.gmt("C:/RFactory/predsurv/vignettes/Gene sets/c6.all.v6.1.symbols.gmt")

results <- read.csv("/home/mtr/rfactory/predsurv/vignettes/iclust2-vignette/summary_iclust2.csv")

results <- read.csv("C:/RFactory/predsurv/performance_results/summary_iclust2.csv")

onco_signature <- ontology_search(hallmark_list = hallmark_all, predictor = results,
                                  coef = median_coeff)

onco_signature_weight <- ontology_search(hallmark_list = hallmark_all, predictor = results,
                                         coef = weight)

capture.output(print(onco_signature), file = "onco_signature.txt")

########################### GSEA
# https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/



iclust2 <- brca[brca$intclust == 2, ]
#rm(brca)


library(glmnet)
library(dplyr)
optimal <- as.matrix(coef(cvfit, s = "lambda.min"))
optimal <- as.data.frame(optimal)
colnames(optimal) <- "mod"
optimal$Hugo_Symbol <- rownames(optimal)
optimal <- optimal %>% filter(mod!=0)



library(glmnet)
cvfit <- readRDS("C:/RFactory/predsurv/performance_results/enet_model_iclust2_filter.RDS")
genedata <- readr::read_tsv("C:/RFactory/Downloads/brca_metabric/data_expression.txt", col_names = TRUE)

genedata$Entrez_Gene_Id[match(names(geneList), genedata$Hugo_Symbol)]

if("age_std" %in% optimal$Hugo_Symbol){
  geneTable <- optimal[(-match(c("age_std", "npi"),optimal$Hugo_Symbol)),]
}else{
  geneTable <- optimal
}


geneList <- geneTable$mod
names(geneList) <- geneTable$Hugo_Symbol

gene <- genedata$Entrez_Gene_Id[match(names(geneList)[abs(geneList)>quantile(abs(optimal$mod),.5)], genedata$Hugo_Symbol)]

gene.df <- clusterProfiler::bitr(gene, fromType = "ENTREZID",
                                 toType = c("ENSEMBL", "SYMBOL"),
                                 OrgDb = org.Hs.eg.db)

head(gene.df)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego))

sort(geneList, decreasing = TRUE)

library(DOSE)


############give ENTREZID NAME
names(geneList) <- genedata$Entrez_Gene_Id[match(names(geneList), genedata$Hugo_Symbol)]

ncg <- gseNCG(sort(geneList, decreasing = TRUE),
           nPerm         = 10,
           minGSSize     =  2,
           maxGSSize =  100,
           pvalueCutoff  = 1,
           pAdjustMethod = "BH",
           verbose       = TRUE)

ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3)


dotplot(ncg, showCategory=30)
enrichMap(ncg, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
cnetplot(ncg, foldChange=geneList)


gene2 <- names(geneList)
ncg2 <- enrichNCG(gene2)
head(ncg2)

gsecc <- gseGO(geneList=sort(geneList, decreasing = TRUE),
               ont="BP",
               OrgDb=org.Hs.eg.db,
               verbose=F)
head(summary(gsecc))

gseaplot(gsecc,  geneSetID="GO:0006414")
plotGOgraph(gsecc)
cnetplot(gsecc, foldChange=geneList)
enrichMap(gsecc, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dotplot(gsecc, showCategory=30)
