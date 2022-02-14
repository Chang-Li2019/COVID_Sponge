library(DESeq2)
library(biomaRt)
setwd("~/Downloads/targetscan_70")
exp_table = read.table('covid19-transcriptomics-pathogenesis-diagnostics-results-master/data/swab_gene_counts.csv',
              sep = ',',header=TRUE)
meta_table <-  read.table('covid19-transcriptomics-pathogenesis-diagnostics-results-master/data/metatable_with_viral_status.csv',
                        sep = ',',header = TRUE)
names(meta_table)
# compare all three groups
dds <- DESeqDataSetFromMatrix(countData=exp_table, 
                              colData=meta_table, 
                              design=~viral_status, tidy = TRUE)
# compare CV2 with no infection (NF)
meta_table.cv2 <- meta_table[(meta_table$viral_status=='SC2') & (meta_table$SC2_rpm>=10),]
meta_table.ctr <- meta_table[(meta_table$viral_status=='no_virus') & (meta_table$SC2_rpm<1),]

meta_table.cv2.NF <- rbind(meta_table.cv2,meta_table.ctr)
keep = c('X')
keep <-  append(keep,meta_table.cv2.NF$CZB_ID)
exp_table.cv2.NF <- exp_table[keep]
names(exp_table.cv2.NF)[1] <- 'Genes'
dds <- DESeqDataSetFromMatrix(countData=exp_table.cv2.NF, 
                              colData=meta_table.cv2.NF, 
                              design=~viral_status, tidy = TRUE)
# Run DESeq
contrast <- c("viral_status", "SC2", "no_virus")
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

# sort by fold change
res <- res[order(res$log2FoldChange,decreasing = TRUE),]
screened = res[(res$log2FoldChange>0 & res$log2FoldChange<1 & res$padj<=0.000001),]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","refseq_mrna"),values=rownames(screened),mart= mart)
write.table(G_list$hgnc_symbol,'Dif_exp_1e-6',quote = FALSE,row.names = FALSE)
# get significantly up-regulated genes.
screened = res[(res$log2FoldChange>0 & res$log2FoldChange<1 & (res$padj>=0.05)),]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","refseq_mrna"),values=rownames(screened),mart= mart)
write.table(G_list$hgnc_symbol,'Dif_exp_ctr_nosig',quote = FALSE,row.names = FALSE)

# get first control
screened = res[(res$log2FoldChange<0 | res$log2FoldChange>1 & (res$padj<=0.000001)),]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","refseq_mrna"),values=rownames(screened),mart= mart)
write.table(G_list$hgnc_symbol,'Dif_exp_ctr_1e-6_0_1',quote = FALSE,row.names = FALSE)

# get control
screened = res[(( abs(res$log2FoldChange)>=1.5) & res$padj<=0.000001),]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","refseq_mrna"),values=rownames(screened),mart= mart)
write.table(G_list$hgnc_symbol,'Dif_exp_ctr_1e-6',quote = FALSE,row.names = FALSE)

rownames(screened)
# convert ensembl ID to gene symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","refseq_mrna"),values=rownames(screened),mart= mart)

write.table(G_list$hgnc_symbol,'Dif_exp_1e-8',quote = FALSE,row.names = FALSE)

# conver miRNA IDs 
library(miRBaseConverter)
acc <- read.table('mir_list',header=FALSE)

tmp <- miRNA_NameToAccession(acc$V1,version = "v22")
write.table(tmp,'converted_mir_list',quote=FALSE,sep='\t',row.names=FALSE)

# check correlation between mRNA expression and viral load
sc2_id <- meta_table$CZB_ID[meta_table$viral_status=='SC2']
sc2_viral <-  meta_table$SC2_rpm[meta_table$viral_status=='SC2']
sc2_cand <- as.data.frame(cbind(sc2_id,sc2_viral))
rownames(sc2_cand) = sc2_cand$sc2_id

gene_list <- c('AGRN','APOL6','CMTR1','CNP','EIF4A2','GBP3','IRF9','JADE2','NUB1','OPTN','PARP10','PARP11','PNPT1','PRMT7','PSMA2','PSMA6','PSMB8','RABGAP1L','RBCK1','SLC25A28','TDRD7','TGM2','TMSB10','TRAFD1','TRIM14','ZNFX1')
l1 <- c()
l2 <- c()
l3 <- c()
for (i in gene_list){
  cand_gene <- unique(G_list$ensembl_gene_id[G_list$hgnc_symbol==i])
  cand_gene_exp <- exp_table.cv2.NF[exp_table.cv2.NF$Genes==cand_gene,]
  t_exp <- t(cand_gene_exp)
  merged <-  merge(sc2_cand, t_exp,by=0)
  merged$sc2_viral <- as.numeric(merged$sc2_viral)
  merged[,4] <- as.numeric(merged[,4])
  corr <- cor.test(x=merged$sc2_viral, y=merged[,4], method = 'spearman')
  l1 <- append(l1,i)
  l2 <- append(l2,corr$estimate)
  l3 <- append(l3,corr$p.value)
  print('---------------------------------')
  print(i)
  print(corr)
}

# quality control

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd,intgroup="viral_status")

