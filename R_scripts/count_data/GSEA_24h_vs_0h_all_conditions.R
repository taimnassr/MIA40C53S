###this is here to do a GSEA of comparisons between 24h and 4h 
##across al strains

count_data = read.csv("count_data/valid_count_data.csv")
head(count_data)
count_data = count_data[,c(1,4:ncol(count_data))]
count_matrix = as.matrix(count_data[,2:ncol(count_data)]) #use ensembl ids
rownames(count_matrix) = count_data$Geneid

col_names = colnames(count_matrix)
condition = stringr::str_replace(col_names, "\\..*", "")
time = stringr::str_extract(col_names, "(\\d)+h")
coldata = data.frame(condition = condition,
                      time = time,
                      row.names = col_names)
all(rownames(coldata)==colnames(count_matrix))
coldata$time = factor(coldata$time,
                      levels = c("4h", "0h","12h","24h"))
coldata_c53s = coldata
coldata_c53s$condition = factor(coldata_c53s$condition,
                                levels = c("C53S", "WT", "Mock"))

coldata_wt = coldata
coldata_wt$condition = factor(coldata_wt$condition,
                                levels = c("WT", "C53S", "Mock"))

coldata_mock = coldata
coldata_mock$condition = factor(coldata_mock$condition,
                                levels = c("Mock", "C53S", "WT"))

##now we run deseq analysis

dds_c53s = DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,
                                          colData = coldata_c53s,
                                            design =~ condition + time + condition:time)
dds_c53s = DESeq2::DESeq(dds_c53s)

dds_mock = DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,
                                          colData = coldata_mock,
                                          design =~ condition + time + condition:time)
dds_mock = DESeq2::DESeq(dds_mock)

dds_wt = DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,
                                          colData = coldata_wt,
                                          design =~ condition + time + condition:time)
dds_wt = DESeq2::DESeq(dds_wt)

##here are the list of the DESEQ log2fc                     

c53s_24_4 = DESeq2::results(dds_c53s,name = "time_24h_vs_4h")
wt_24_4 = DESeq2::results(dds_wt, name = "time_24h_vs_4h")
mock_24_4 = DESeq2::results(dds_mock, name = "time_24h_vs_4h")

###now we arrange them according to the stats argument
c53s_24_4 = c53s_24_4[order(c53s_24_4$stat, decreasing = TRUE),]
wt_24_4 = wt_24_4[order(wt_24_4$stat, decreasing = TRUE),]
mock_24_4 = mock_24_4[order(mock_24_4$stat, decreasing = TRUE),]

###now we do the GSEA
##convert gene names to ensembl ids
gsea_c53s = as.data.frame(c53s_24_4[,"stat"])
rownames(gsea_c53s) = rownames(c53s_24_4)

gsea_wt = as.data.frame(wt_24_4[,"stat"])
rownames(gsea_wt) = rownames(wt_24_4)

gsea_mock = as.data.frame(mock_24_4[,"stat"])
rownames(gsea_mock) = rownames(mock_24_4)



library(clusterProfiler)
library(org.Hs.eg.db)

#convert into names vectors
gsea_c53s = setNames(gsea_c53s[[1]], rownames(gsea_c53s))
gsea_wt = setNames(gsea_wt[[1]], rownames(gsea_wt))
gsea_mock = setNames(gsea_mock[[1]], rownames(gsea_mock))

names(gsea_c53s) = sub("\\.\\d+$","", names(gsea_c53s))
names(gsea_wt) = sub("\\.\\d+$","", names(gsea_wt))
names(gsea_mock) = sub("\\.\\d+$","", names(gsea_mock))
gsea_result_c53s = gseGO(geneList = gsea_c53s,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",keyType = "ENSEMBL",
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05)
gsea_result_c53s = as.data.frame(gsea_result_c53s)

gsea_result_wt = gseGO(geneList = gsea_wt,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",keyType = "ENSEMBL",
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05)
gsea_result_wt = as.data.frame(gsea_result_wt)

gsea_result_mock = gseGO(geneList = gsea_mock,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",keyType = "ENSEMBL",
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05)
gsea_result_mock = as.data.frame(gsea_result_mock)
####now the graphs