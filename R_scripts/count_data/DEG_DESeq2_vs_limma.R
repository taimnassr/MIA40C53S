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



############so now after getting the DEG
#we want to see the RPKM of the top DEG accross the different
#conditions (WT, Mock, and C53S 24h vs 0h)
##first we need to load the DEG and extract them

list.files()
deg_wt = read.csv("count_data/DESeq2_results_WT_24h_vs_0h.csv")
deg_wt = deg_wt[abs(deg_wt$log2FoldChange) > 0.5&
                deg_wt$padj < 0.05,]
deg_wt = na.omit(deg_wt)
#now let us extract the top 5 DEG
deg_wt = deg_wt[order(abs(deg_wt$log2FoldChange),
                      decreasing = TRUE),]
genes_wt = deg_wt$X[1:5]#for the wt

deg_mock = read.csv("count_data/DESeq2_results_MOCK_24h_vs_0h.csv")
deg_mock = deg_mock[abs(deg_mock$log2FoldChange) > 0.5&
                  deg_mock$padj < 0.05,]
deg_mock = na.omit(deg_mock)
#now let us extract the top 5 DEG
deg_mock = deg_mock[order(abs(deg_mock$log2FoldChange),
                      decreasing = TRUE),]
genes_mock = deg_mock$X[1:5]#for the wt

deg_c53s = read.csv("count_data/DESeq2_results_C53S_24h_vs_0h.csv")
deg_c53s = deg_c53s[abs(deg_c53s$log2FoldChange) > 0.5&
                  deg_c53s$padj < 0.05,]
deg_c53s = na.omit(deg_c53s)
#now let us extract the top 5 DEG
deg_c53s = deg_c53s[order(abs(deg_c53s$log2FoldChange),
                      decreasing = TRUE),]
genes_c53s = deg_c53s$X[1:5]#for the wt

#now we combine the gene lists
genes_deg = unique(c(genes_wt, genes_c53s, genes_mock))


#now we want to plot the RPKM of these genes
#using the different conditions 

#first load the RPKM file
rpkm = read.csv("count_data/rpkm.csv")
names = names(rpkm)
#get comparisons of interest
conditions = names[grepl("0h|24h", names)]
gene_to_ensembl = read.csv("count_data/ensembl_to_gene_conversion.csv")
head(gene_to_ensembl)
rpkm = merge(rpkm,
             gene_to_ensembl,
             by.x = "X",
             by.y = "V1")
rpkm = rpkm[, c("V2", conditions)]
names(rpkm)[1] = "gene_name"
rpkm = rpkm[rpkm$gene_name %in% genes_deg,]
### so now we have our rpkm

##we change to long data format
library(tidyr)
library(stringr)
rpkm_long1 = pivot_longer( data = rpkm,
                           cols = contains("h"),
                           names_to = "condition_time",
                           values_to = "rpkm")

rpkm_long1 = as.data.frame(rpkm_long1)
rpkm_long1$condition = sub("\\..*$", "", rpkm_long1$condition_time)
rpkm_long1$time = str_extract(rpkm_long1$condition_time,
                              "\\d+h")
###now we can start graphing
library(ggplot2)
rpkm_graph = ggplot(rpkm_long1,
                    aes( x = gene_name,
                         y =log10(rpkm),
                         colour = time)) +
  geom_point() +
  labs(title = "RPKM for top DEG",
       x= "Gene",
       y = "Log10(RPKM)") +
  facet_grid( condition~.)
rpkm_graph

pdf("plots/RPKM_comparison_top_DEG_WT_C53S_Mock.pdf",
    width = 12, height = 8)
rpkm_graph
dev.off()


#what changes are specific to C53S between 24h vs 0h
#I kinda wanna see which genes     change differently

genes_c53s = deg_c53s$X

genes_wt = deg_wt$X
genes_wt
specific_c53s = setdiff(genes_c53s, genes_wt)
specific_c53s

###now we want to  only plot the specific c53s
head(rpkm)
rpkm_c53s = rpkm[rpkm$gene_name %in% specific_c53s,]
rpkm_c53s_l = data.frame(pivot_longer(rpkm_c53s,
                           cols = 2:ncol(rpkm_c53s),
                           names_to = "time_condition",
                           values_to = "rpkm"))
head(rpkm_c53s_l)
rpkm_c53s_l$condition = sub("\\..*$", 
                            "",
                            rpkm_c53s_l$time_condition)

rpkm_c53s_l$time = str_extract(rpkm_c53s_l$time_condition,
                               "\\d+h")

rpkm_graph_c53s = ggplot(rpkm_c53s_l,
                    aes( x = gene_name,
                         y =log10(rpkm),
                         colour = time)) +
  geom_point() +
  labs(title = "RPKM for top DEG",
       x= "Gene",
       y = "Log10(RPKM)") +
  facet_grid( condition~.)
##################
#### I wanna redo this according to 
## a tutorial that matthias sent
# which is found on the following link

## https://www.atakanekiz.com/technical/a-guide-to-designs-and-contrasts-in-DESeq2/

#so first I will only use the comparisons between
# wt and C53S since they are the most similar

count_data = read.csv("count_data/valid_count_data.csv")
##now I only wanna keep the data of C53S or WT
names = names(count_data)
names_keep = names[grepl("C53S|WT",names)]
count_wt_c53s = count_data[,c("Gene_name", 
                              names_keep)]##keep the relevant
count_matrix = as.matrix(count_wt_c53s[,
                          2:ncol(count_wt_c53s)])
rownames(count_matrix) = count_wt_c53s$Gene_name##
##now let us create to col_data
library(stringr)
col_names = colnames(count_matrix)
condition = str_replace_all(col_names,
                            "\\..*$", "")
time = str_extract(col_names,
                   "\\d+h")
coldata = data.frame(condition = factor(condition),
                     time = factor(time),
                     row.names = col_names)

#####now we can start with the dds.
###but here we want to use the function that
## is supplied by this website I am following

contraster <- function(dds,    # should contain colData and design
                       group1, # list of character vectors each with 2 or more items 
                       group2, # list of character vectors each with 2 or more items
                       weighted = F
){
  
  
  mod_mat <- model.matrix(design(dds), colData(dds))
  
  grp1_rows <- list()
  grp2_rows <- list()
  
  
  for(i in 1:length(group1)){
    
    grp1_rows[[i]] <- colData(dds)[[group1[[i]][1]]] %in% group1[[i]][2:length(group1[[i]])]
    
  }
  
  
  for(i in 1:length(group2)){
    
    grp2_rows[[i]] <- colData(dds)[[group2[[i]][1]]] %in% group2[[i]][2:length(group2[[i]])]
    
  }
  
  grp1_rows <- Reduce(function(x, y) x & y, grp1_rows)
  grp2_rows <- Reduce(function(x, y) x & y, grp2_rows)
  
  mod_mat1 <- mod_mat[grp1_rows, ,drop=F]
  mod_mat2 <- mod_mat[grp2_rows, ,drop=F]
  
  if(!weighted){
    
    mod_mat1 <- mod_mat1[!duplicated(mod_mat1),,drop=F]
    mod_mat2 <- mod_mat2[!duplicated(mod_mat2),,drop=F]
    
  }
  
  return(colMeans(mod_mat1)-colMeans(mod_mat2))
  
  
}
library(DESeq2)
coldata$condition = factor(coldata$condition,
                           levels = c("WT", "C53S"))
rownames(count_matrix) = make.unique(rownames(count_matrix))
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition + time + condition:time)
dds <- DESeq(dds)
resultsNames(dds)
res_c53s_24 = results(dds,
  contrast = list(c("time_24h_vs_0h",
   "conditionC53S.time24h")))## so this is what we need
res_c53s_24 = as.data.frame(res_c53s_24)
write.csv(res_c53s_24,
          "count_data/DESeq2_results_C53S_24h_vs_0h_vs_WT.csv",
          row.names = TRUE)

res_c53s_24 = na.omit(res_c53s_24)
deg_c53s_24 = res_c53s_24[abs(res_c53s_24$log2FoldChange) > 0.5 &
                            res_c53s_24$padj < 0.05,]
rownames(deg_c53s_24)

deg_c53s = read.csv("count_data/DESeq2_results_C53S_24h_vs_0h.csv")
head(deg_c53s)
deg_c53s = deg_c53s[abs(deg_c53s$log2FoldChange) > 0.5&
                   deg_c53s$padj < 0.05, ]
deg_c53s = na.omit(deg_c53s)
head(deg_c53s)
intersect(deg_c53s$X,
          rownames(deg_c53s_24))
nrow(deg_c53s)
nrow(deg_c53s_24)
setdiff(rownames(deg_c53s_24),
        deg_c53s$X)


deseq_c53s_24_2 = results(dds,
                          contrast = contraster(dds,
              group1 = list(c("condition", "C53S"),
                            c("time", "24h")),
              group2 = list(c("condition", "C53S"),
                            c("time", "0h"))))

deseq_c53s_24_2 = as.data.frame(deseq_c53s_24_2)
head(deseq_c53s_24_2)
identical(deseq_c53s_24_2, res_c53s_24)
deseq_c53s_24_2[2,]
res_c53s_24[2,]

res_c53s_24$gene = rownames(res_c53s_24)
deseq_c53s_24_2$gene = rownames(deseq_c53s_24_2)
deseq_c53s_24_2 = na.omit(deseq_c53s_24_2)
identical(res_c53s_24[do.call(order, res_c53s_24),],
          deseq_c53s_24_2[do.call(order, deseq_c53s_24_2),])

nrow(res_c53s_24)
nrow(deseq_c53s_24_2)
