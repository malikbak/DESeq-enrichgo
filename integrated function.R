library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
enrichment<-function(x,y){
  plot=enrichGO(
    x,
    org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = y,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  dotplot(plot)
}
data("geneList")
eid <- readxl::read_xlsx("Entrez_ID.xlsx")
DEgs <- read.csv('EnID.csv')
gene = names(geneList)[geneList > 1]
enrichment(eid$To, "CC")
Kegg <- enrichKEGG(eid$To)
Kegg@result
dotplot(Kegg)
kegg_Re <- Kegg@result
library(ggplot2)
top_20 <- kegg_Re[1:10,]
ggplot(top_20, aes(pvalue, Description, size= pvalue, color= Count))+
  geom_point()
MF <-enrichGO(
  DEgs$To,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = 'MF',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
barplot(MF)
MF_res <- MF@result
MF_res <- MF_res[1:20,]
library(ggplot2)
ggplot(MF_res, aes(pvalue, Description, size= pvalue, color= Count))+
  geom_point()
##################BP#####################
BP <-enrichGO(
  DEgs$To,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
barplot(BP)
BP_res <- BP@result
BP_res <- BP_res[1:20,]
library(ggplot2)
ggplot(BP_res, aes(pvalue, Description, size= pvalue, color= Count))+
  geom_point()
#################CC############
CC <-enrichGO(
  ttg$entrezgene_id,
  org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)
dotplot(CC)
CC_res <- CC@result
CC_res <- CC_res[1:10,]
library(ggplot2)
ggplot(CC_res, aes(pvalue, Description, size= pvalue, color= Count))+
  geom_point()
######################
library(enrichR)
dbs_pw <- "KEGG_2019_Human"
upEnriched_pw <- enrichr(genes = eid$From, databases = dbs_pw)
dnEnriched_pw <- enrichr(genes = dn.genes, databases = dbs_pw)

plotEnrich(upEnriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Genes - KEGG Pathway Analysis")
