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

GSEA_KEGG_df_to_dotplot <- function(GSEA_df, 
                                    num_entries =20, title){
  ##this is a function that takes a DATA frame of the
  #results of the GSEA or KEGG and then plots the dotplot
  #it all takes how many GO terms we need to graph,
  #the default is num_entries = 20
  
  ##make sure the input is converted to dataframe before
  
  
  
  ##here you have to make sure that the NES column, p.adj, 
  # setsize, and Description exists 
  
  ##title tell us what the title of the graph should be
  colnames_to_check = c("Description", "NES", "p.adjust",
                        "setSize")
  
  if (!all(colnames_to_check %in% colnames(GSEA_df))){
    stop("Required Columns are not present in the dataframe")
  }
  
  ###take the top num_entries entries based on their ABSOLUTE
  #NES value
  
  ##to make sure there is equal representations
  # we will take the top 10 upregulated and top 10 downregulated
  
  GSEA_pos = GSEA_df[GSEA_df$NES > 0,]
  GSEA_top = GSEA_pos[order(abs(GSEA_pos$NES),
                            decreasing = TRUE),]
  
  GSEA_top = GSEA_top[1:(num_entries/2),] ## here the top 10 upregulated terms
  
  GSEA_neg = GSEA_df[GSEA_df$NES < 0,]
  GSEA_btm = GSEA_neg[order(abs(GSEA_neg$NES),
                            decreasing = TRUE),]
  
  GSEA_btm = GSEA_btm[1:(num_entries/2),] ## here the top 10 downregulated terms
  GSEA_input = rbind(GSEA_top, GSEA_btm)
  ##then we order the input
  ## we order the description based on NES Value
  #this is important in the graph for factoring
  
  #remove rows containing NA
  GSEA_input = na.omit(GSEA_input)
  
  GSEA_input$Description = factor(GSEA_input$Description,
                                  levels = GSEA_input$Description[order(GSEA_input$NES)])
  
  
  suppressPackageStartupMessages(library(ggplot2))
  ##maybe we need to tranform NES to visualize better
  
  dotplot = ggplot(GSEA_input, aes(x = NES,
                                   y = Description,
                                   size = setSize,
                                   color = p.adjust))+
    scale_size_continuous(range = c(1,10),
                          limits = c(min(GSEA_input$setSize),
                                     max(GSEA_input$setSize)),
                          breaks = seq(min(GSEA_input$setSize),
                                       max(GSEA_input$setSize),
                                       by = 80))+
    geom_point() + 
    labs( x = "Normalized Enrichment Value (NES)",
          y = "GO Term",
          title = title,
          size = "Number of genes per Go Term",
          color = "Adjusted P.value")+ theme_minimal()+
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 14)) +
    scale_color_gradient(low = "blue",
                         high = "red",
                         limits = c(0,0.05)) 
  
  return(dotplot)
}

plot_c53s = GSEA_KEGG_df_to_dotplot(gsea_result_c53s, 20,
                                    "GSEA C53S 24h vs 4h")
pdf("plots/GSEA/GSEA_C53S_24h_vs_4h.pdf", width = 10, heigh=8)
plot_c53s
dev.off()

plot_wt = GSEA_KEGG_df_to_dotplot(gsea_result_wt, 20,
                                    "GSEA WT 24h vs 4h")
pdf("plots/GSEA/GSEA_WT_24h_vs_4h.pdf", width = 10, heigh=8)
plot_wt
dev.off()

plot_mock = GSEA_KEGG_df_to_dotplot(gsea_result_mock, 20,
                                  "GSEA Mock 24h vs 4h")
pdf("plots/GSEA/GSEA_MOck_24h_vs_4h.pdf", width = 10, heigh=8)
plot_mock
dev.off()

