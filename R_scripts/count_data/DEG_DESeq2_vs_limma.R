###### here we will try to compare
## the C53S 24 h vs C53S 0 h based on our RNA-seq data
#builidng  a model taking into account
# the conditions : wt, mock, C53S
#timepoints: 0h,4h,12h, and 24h

##here we will be comparing two methods
#DESeq2 and Limma + voom 
### let us load the data

setwd("/mnt/d/Data_Analysis/NGS_KW_MIA40C53S/count_data")
count_data = read.csv("valid_count_data.csv",
                      header = TRUE)
count_matrix = as.matrix(count_data[,
                          c(4:ncol(count_data))])
row.names(count_matrix) = count_data$Gene_name
#### from the count_matrix we want to make the
###col data file that is needed for DESeq2 input

col_names = colnames(count_matrix)
###first extract the conditions
condition = sub("\\..*", "", col_names)
time_point = sub("^[^.]*\\.(\\d{1,2}).*", "\\1", col_names)
time_point = as.numeric(time_point)
repitition = sub(".*\\.(\\d+)$", "\\1", col_names)
repitition = as.numeric(repitition)
#### so now we have our different criteria

col_data = data.frame(
  condition = condition,
  time_point = time_point,
  repitition = repitition,
  row.names = col_names
)
col_data$condition = factor(col_data$condition,
                            levels = c("C53S", "WT", "Mock"))
col_data$time_point = factor(col_data$time_point,
                                levels = c(0,4,12,24))

###check now if count matrix aligns with col_data
all(colnames(count_matrix) == rownames(col_data))#should be TRUE

##for count_matrix we need to remove duplicat names
#to do that we will make them unique

row_names = rownames(count_matrix)#check if they are unique
any(duplicated(row_names))##if TRUE then change them
row_names2 = make.unique(row_names)
rownames(count_matrix) = row_names2# now they are unique
# now we build or DESEQ models
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition + 
             time_point +
              condition:time_point) ## here the reference
#was C53S because we were interesed in that comparison
dds <- DESeq(dds) # now we have our results
resultsNames(dds)# the reference is C53S here 

#### now extract our comparison
DESeq_C53S_24_0 = results(dds, name = "time_point_24_vs_0")
DESeq_C53S_24_0 = as.data.frame(DESeq_C53S_24_0)
## now we want to remove NA and extract
#genes whose abs log2fc are more than 0.5 and
# adj p value less tham 
# we chose 0.5 instead of 1 because there were no results
#for when we chose 1
write.csv(DESeq_C53S_24_0, "DESeq2_results_C53S_24h_vs_0h.csv",
          row.names = TRUE)
### let us find now the DEG
DEG_DESeq_24_0 = DESeq_C53S_24_0[abs(DESeq_C53S_24_0$log2FoldChange) >= 0.5 &
                                 DESeq_C53S_24_0$padj < 0.05,]
DEG_DESeq_24_0 = na.omit(DEG_DESeq_24_0)
write.csv(DEG_DESeq_24_0, "DEG_DESeq2_C53S_24h_vs_0h.csv",
          row.names = TRUE)


####now let us do a GSEA and graph it 
### we first need to get the entrez ids

library(biomaRt)
mart <- useEnsembl(biomart = "genes",
                   dataset = "hsapiens_gene_ensembl",
                   mirror = "asia")
gene_symbols = rownames(DESeq_C53S_24_0)
gene_to_id = getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = mart
)
DESeq_C53S_24_0$gene = rownames(DESeq_C53S_24_0)
KEGG_DESeq = merge(DESeq_C53S_24_0,
                   gene_to_id,
                   by.x = "gene",
                   by.y = "hgnc_symbol")
KEGG_DESeq_input = KEGG_DESeq$stat ## arrange based on stat
names(KEGG_DESeq_input) = KEGG_DESeq$entrezgene_id
KEGG_DESeq_input = sort(KEGG_DESeq_input, decreasing = TRUE)
### now we do the KEGG analysis
KEGG_DESeq_input = KEGG_DESeq_input[!duplicated(names(KEGG_DESeq_input))]
KEGG_DESeq_input = na.omit(KEGG_DESeq_input)
KEGG_DESeq_input = KEGG_DESeq_input[!is.na(names(KEGG_DESeq_input))]

library(clusterProfiler)
KEGG_res_DESeq_24_0 = gseKEGG(geneList = KEGG_DESeq_input,
                           organism= "hsa",
                           keyType = "kegg",
                           minGSSize = 5,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           verbose = TRUE,
                           eps = 0) ##
KEGG_res_DESeq_24_0 = as.data.frame(KEGG_res_DESeq_24_0)
##### now we plot the results
DESeq2_KEGG = GSEA_KEGG_df_to_dotplot(GSEA_df = KEGG_res_DESeq_24_0,
                                      title = "KEGG C53S 24h vs 0H DESeq2") +
  xlim(-2.5,0)
pdf("Dotplot_KEGG_C53S_24h_vs_0h_DESeq2.pdf", width = 10, height = 8)
DESeq2_KEGG
dev.off()


##### after that we wanna save the KEGG in an excel
## we need to modify it
head(DEG_DESeq_24_0)
DEG_DESeq_24_0$X = rownames(DEG_DESeq_24_0)
new_KEGG_DESeq2 = task_entrez_to_DEG(KEGG_res_DESeq_24_0,
                                     gene_to_id,
                                     DEG_DESeq_24_0)
head(new_KEGG_DESeq2)
write.csv(new_KEGG_DESeq2, "KEGG_DESeq2_C53S_24h_vs_0h.csv",
          row.names = TRUE)


##### now first we want to the volcano plot for the
#DEG results for both C53S 24h vs 0h
#and I think also WT and MOCK 24h vs 0h

###let me load the count matrix and col_data
#from above



#####here we have to focus on the comparison
##for the col_data if we want to compare
#WT 24h vs 0H then the levels for the condition part should
#be WT
#and so on

#using the same matrix but different col_data
#with different levels
head(col_data)
col_data_c53s = col_data
col_data_wt = col_data
col_data_mock = col_data
col_data_c53s$condition = factor(col_data_c53s$condition,
                                 levels = c("C53S", "WT","Mock"))
#here now the comparison would be to C53S

col_data_wt$condition = factor(col_data_wt$condition,
                               levels = c( "WT","C53S","Mock"))
#here would be for wt 

col_data_mock$condition = factor(col_data_mock$condition,
                               levels = c("Mock" ,  "WT","C53S"))
#here would be for mock
#now we do the dds for all three
library(DESeq2)
dds_c53s = DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = col_data_c53s,
                                  design = ~ condition +
                                    time_point +
                                    condition:time_point)
dds_c53s = DESeq(dds_c53s)
dds_mock = DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = col_data_mock,
                                  design = ~ condition +
                                    time_point +
                                    condition:time_point)
dds_mock = DESeq(dds_mock)
dds_wt = DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = col_data_wt,
                                  design = ~ condition +
                                    time_point +
                                    condition:time_point)  
dds_wt =DESeq(dds_wt)##
###since I am interested in 24h vs 0 h comparison
dds_wt_24_vs_0 = results(dds_wt,
                         name = "time_point_24_vs_0" )
dds_wt_24_vs_0 = as.data.frame(dds_wt_24_vs_0)
write.csv(dds_wt_24_vs_0,
          "DESeq2_results_WT_24h_vs_0h.csv",
          row.names = TRUE)
dds_mock_24_vs_0 = results(dds_mock,
                           name = "time_point_24_vs_0" )
dds_mock_24_vs_0 = as.data.frame(dds_mock_24_vs_0)
write.csv(dds_mock_24_vs_0,
          "DESeq2_results_MOCK_24h_vs_0h.csv",
          row.names = TRUE)
dds_c53s_24_vs_0 = results(dds_c53s,
                           name = "time_point_24_vs_0",
                           independentFiltering = TRUE)
dds_c53s_24_vs_0 = as.data.frame(dds_c53s_24_vs_0)
head(dds_c53s_24_vs_0[order(dds_c53s_24_vs_0$log2FoldChange,
                           decreasing = TRUE),])


#############now let us plot the volcano plot
#using the function volcano plots
#we also wanna highlight the DEG
####now we change the volcano plots
## to account for significance


volcao_c53s = volcano_plots(df = dds_c53s_24_vs_0,
               title =  "C53S 24h vs 0h",
               log2fc_threshold = 0.5,
               padj_threshold =  0.05  ) +
  xlim(-1.5,3) + ylim(0,100)
volcao_c53s 


pdf("Volcano_plot_C53S_24h_vs_0h.pdf",
    width = 8, height = 6)
volcao_c53s
dev.off()

volcao_wt = volcano_plots(df = dds_wt_24_vs_0,
                            title =  "WT 24h vs 0h",
                            log2fc_threshold = 0.5,
                            padj_threshold =  0.05  ) +
  xlim(-1.5,3) + ylim(0,100)
volcao_wt 

pdf("Volcano_plot_WT_24h_vs_0h.pdf",
    width = 8, height = 6)
volcao_wt
dev.off()


volcao_mock = volcano_plots(df = dds_mock_24_vs_0,
                          title =  "Mock 24h vs 0h",
                          log2fc_threshold = 0.5,
                          padj_threshold =  0.05  ) +
  xlim(-1.5,1.5) + ylim(0,60)
volcao_mock

pdf("Volcano_plot_Mock_24h_vs_0h.pdf",
    width = 8, height = 6)
volcao_mock
dev.off()



head(dds_c53s_24_vs_0)
