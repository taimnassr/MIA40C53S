setwd("/mnt/d/Data_Analysis/NGS_KW_MIA40C53S/count_data")#where the data is stored

#now we do input the count_data
count_data = read.csv("valid_count_data.csv", header = T)

#we need to change the ensembl geneID to gene names using a file
#that was already provided

ensmbl_to_gene = read.delim("/mnt/c/NGS_files/genomes/gencode_human_release47/gencode.v47.primary_assembly.annotation.gene_id_to_name.tsv",
                            header = F)
#now match the ensembl id to the gene names according to the table provided
count_data_ann = merge(count_data, ensmbl_to_gene,
                       by.x = "Geneid", by.y = "V1", all.x = TRUE)

#change the arrangement 
count_data_ann = count_data_ann[,c(1,ncol(count_data_ann),
                                   2:ncol(count_data_ann)- 1)]
colnames(count_data_ann)[2] = "Gene_name"

write.csv(count_data_ann,"valid_count_data.csv", row.names = FALSE)# save
#annontated file 

#check for duplicated genes in the count_data
head(count_data_ann)
sum(duplicated(count_data_ann$Gene_name))# we have 15 duplicated genes
dupl_genes = count_data_ann$Gene_name[
  duplicated(count_data_ann$Gene_name)
]

#since thex correspond to different ensembl IDs, I will just then
# rename the genes

count_data2 =  as.data.frame(count_matrix)
gene_names = rownames(count_data2)
rownames(count_data2) = make.unique(gene_names)
###so  now after annotating the files, we want to do differential gene expression
#analysis 
# we have three different conditions accross different timepoints
#mutant C53S, wt, and mock
# for each condition we have at 0, 4, 12, 24 hours

#create_count matrix for the genes 
head(count_data_ann)
#the data should be as a matrix
count_matrix = as.matrix(count_data_ann[,4:ncol(count_data_ann)])
rownames(count_matrix) = count_data_ann$Gene_name

## the count matrix 2 will be used with the unique gene names
#using countdata2

count_matrix2 = as.matrix(count_data2)
#now do the col annotation
colnames = colnames(count_matrix)

condition = sub("\\..*", "", colnames)
library(stringr)
time = str_extract(colnames, "\\d+(?=h)")

#so now we got the time and condition
rownames = colnames


colData = data.frame(
  row.names = rownames,
  condition = factor(condition),
  time = factor(time, levels = c("0", "4","12","24"))
)
all(names(count_matrix) == rownames(colData)) # important checkup
#now we have colData, let us do the dds
library(DESeq2)
dds = DESeqDataSetFromMatrix(countData = count_matrix2,
                             colData = colData,
                             design = ~condition)
dds = DESeq(dds)
#we did the dds 

#now let us take the different comparisons
res_c53s_wt = results(dds, 
                      contrast = c("condition",
                                   "C53S",
                                   "WT"))
res_wt_mock = results(dds, 
                      contrast = c("condition",
                                   "WT",
                                   "Mock"))
res_c53_mock = results(dds, 
                       contrast = c("condition",
                                    "C53S",
                                    "Mock"))
######so now we did this let us see the results
#we change them into data.frames
df_res_c53_mock = as.data.frame(res_c53_mock)
df_res_wt_mock = as.data.frame(res_wt_mock)
df_res_c53s_wt = as.data.frame(res_c53s_wt)

### now let us save them 
write.csv(df_res_c53_mock, "c53s_vs_mock_log2fc_raw.csv",
  row.names = T )
write.csv(df_res_wt_mock, "wt_vs_mock_log2fc_raw.csv",
          row.names = T )
write.csv(df_res_c53s_wt, "c53s_vs_wt_log2fc_raw.csv",
          row.names = T )

#### now create the volcano plot for analysis
df_res_c53_mock # this is the dataframe

library(ggplot2)

volcano_plots = function(df, condition_name){
  
  #df should be containing DESEQ log2fc results as dataframe
  #there should be a column with the log2foldchange
  #called log2FoldChange  and a column called padj
  
  #conditon name is like mock vs wt 
  
  library(ggplot2)
  
  title = paste ("DESeq2 Results", condition_name )
  
          plot = ggplot(df, 
                       aes(x = log2FoldChange,
                           y = -log10(padj))) +
  geom_point() + labs(title = title,
                      x = "Log2FC",
                      y = "-log10(P. adjusted)")+
  geom_vline(xintercept = c(-1,1), color = "blue")+
  geom_hline(yintercept = 1.301, color = "blue") +
  ylim(c(0,200)) + 
  scale_x_continuous(breaks = seq(-5,5, by  = 1),
                     labels = seq(-5,5, by = 1),
                     limits = c(-5,5)) 
return(plot)
}

plot_c53_mock = volcano_plots(df_res_c53_mock,
                              "C53S vs Mock")
pdf("Volcano_plot_DESeq2_C53S_vs_Mock.pdf", width = 8,
    height = 6)
plot_c53_mock
dev.off()
plot_c53_wt = volcano_plots(df_res_c53s_wt,
                            "C53S vs WT")
pdf("Volcano_plot_DESeq2_C53S_vs_WT.pdf", width = 8,
    height = 6)
plot_c53_wt
dev.off()
plot_wt_mock = volcano_plots(df_res_wt_mock,
                             "WT vs Mock")
pdf("Volcano_plot_DESeq2_WT_vs_Mock.pdf", width = 8,
    height = 6)
plot_wt_mock
dev.off()


###now we extract the differentially expressed genes
##so the standard would be a log2fc larger than 1
#and p.adj value smaller than 0.5

sig_c53_wt = df_res_c53s_wt[abs(df_res_c53s_wt$log2FoldChange) >= 1 &
                             df_res_c53s_wt$padj < 0.05,  ]
sig_c53_mock = df_res_c53_mock[abs(df_res_c53_mock$log2FoldChange) >= 1 &
                              df_res_c53_mock$padj < 0.05,  ]
sig_wt_mock = df_res_wt_mock[abs(df_res_wt_mock$log2FoldChange) >= 1 &
                                 df_res_wt_mock$padj < 0.05,  ]

write.csv(sig_c53_mock, 
          "sig_DEG_C53S_vs_Mock.csv",
          row.names = TRUE)    

write.csv(sig_wt_mock, 
          "sig_DEG_WT_vs_Mock.csv",
          row.names = TRUE) 

write.csv(sig_c53_wt, 
          "sig_DEG_C53S_vs_WT.csv",
          row.names = TRUE) 

#######now first create venndiagrams based on the genes themselves

genes_c53s_wt = row.names(sig_c53_wt)
genes_c53s_mock = row.names(sig_c53_mock)
genes_wt_mock = row.names(sig_wt_mock)


library(VennDiagram)# for plotting

venn_digram = venn.diagram(x = list( Set1 = genes_c53s_wt,
                                     Set2 = genes_c53s_mock,
                                     Set3 = genes_wt_mock),
                           filename = NULL,
                           col = "black",
                           fill = c("red","blue", "green"),
                           main = "DEG among WT, MOCK, and C53S",
                           category.names = c("C53S vs WT",
                                              "C53S vs Mock",
                                              "WT vs Mock"),
                           main.cex = 2,
                           cat.cex = 1.2,
                           cex = 1.2)
grid.newpage()
grid.draw(venn_digram)

pdf("DEG_venndiagram_WT_Mock_C53S.pdf", width = 9,height = 9)
grid.draw(venn_digram)
dev.off()


####let us for now save the list of DEG
DEG_c53s_wt = sig_c53_wt[,c("log2FoldChange",
                            "padj")]
write.csv(DEG_c53s_wt, "sig_DEG_C53S_vs_WT.csv", row.names = T)
DEG_c53s_mock = sig_c53_mock[,c("log2FoldChange",
                              "padj")]
write.csv(DEG_c53s_mock, "sig_DEG_C53S_vs_Mock.csv", row.names = T)
DEG_wt_mock = sig_wt_mock[,c("log2FoldChange",
                            "padj")]
write.csv(DEG_wt_mock, "sig_DEG_WT_vs_Mock.csv", row.names = T)

### we also want to do GSEA analysis

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("KEGG.db") 
BiocManager::install("KEGGREST", force = TRUE)


library(clusterProfiler)
library(org.Hs.eg.db) # human
library(KEGGREST)

#now we start GSEA
#we need to map the genes to ensembl id
#but we already have the table
#we need all GENES


write.csv(as.data.frame(ensmbl_to_gene), 
          "esembl_to_gene_conversion.csv",
          row.names = TRUE)
##let us read the files then
df_res_c53s_wt = read.csv("c53s_vs_wt_log2fc_raw.csv",
                          header = TRUE)
df_res_c53s_mock = read.csv("c53s_vs_mock_log2fc_raw.csv",
                            header = TRUE)
df_res_wt_mock = read.csv("wt_vs_mock_log2fc_raw.csv",
                            header = TRUE)

#we need to convert gene names to ensembl IDS using the given
#file

ensmbl_to_gene = read.csv("ensembl_to_gene_conversion.csv",
                          header = TRUE)


###now we merge to create the dataframes needed for GSEA
#and we will use the stats in the results to generate
#the ranked gene list needed for GSEA

## so basically what is done down id for the ensembl ids
## which is not suitable for kegg.
##thereffore we need to redo it after for the entrez gene ids
GSEA_C53S_WT = merge(df_res_c53s_wt, ensmbl_to_gene,
                    by.x = "X",
                    by.y = "V2")


GSEA_C53S_WT_input = GSEA_C53S_WT[,"stat"]# this what we need
names(GSEA_C53S_WT_input) = GSEA_C53S_WT$V1 # row names as ensmebl id
#now we sort by stat
GSEA_C53S_WT_input = sort(GSEA_C53S_WT_input, decreasing = T)
#this is the input



GSEA_C53S_MOCK = merge(df_res_c53s_mock, ensmbl_to_gene,
                     by.x = "X",
                     by.y = "V2")

GSEA_C53S_MOCK_input = GSEA_C53S_MOCK[,"stat"]# this what we need
names(GSEA_C53S_MOCK_input) = GSEA_C53S_MOCK$V1 # row names as ensmebl id
#now we sort by stat
GSEA_C53S_MOCK_input = sort(GSEA_C53S_MOCK_input, decreasing = T)


GSEA_WT_MOCK = merge(df_res_wt_mock, ensmbl_to_gene,
                       by.x = "X",
                       by.y = "V2")

GSEA_WT_MOCK_input = GSEA_WT_MOCK[,"stat"]# this what we need
names(GSEA_WT_MOCK_input) = GSEA_WT_MOCK$V1 # row names as ensmebl id
#now we sort by stat
GSEA_WT_MOCK_input = sort(GSEA_WT_MOCK_input, decreasing = T)

## so now we have the sorted lists per conditions according
# to the stats columns and with ensembl ids as rows
##just quickly we need to modify the ensemblids
#from e4623919191.12 to ee462391919112



names(GSEA_C53S_WT_input) =  
  sub("\\.\\d+$", "", names(GSEA_C53S_WT_input))

names(GSEA_C53S_MOCK_input) =  
  sub("\\.\\d+$", "", names(GSEA_C53S_MOCK_input))

names(GSEA_WT_MOCK_input) =  
  sub("\\.\\d+$", "", names(GSEA_WT_MOCK_input))

### let us try mapping out ids
## since we need entrez IDS for the KEGG
#we will convert ensembl to entrez using biomart

BiocManager::install("biomaRt")

if (!dir.exists(Sys.getenv("R_LIBS_USER"))) dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))

BiocManager::install("biomaRt", lib=Sys.getenv("R_LIBS_USER"),
                     force = TRUE)
library(biomaRt)

##let us convert ensembl to entrez
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ids_c53s_wt = names(GSEA_C53S_WT_input)
mapp_c53s_wt = getBM(attributes = c("ensembl_gene_id",
                                    "entrezgene_id",
                                    "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = ids_c53s_wt,
                     mart = ensembl)
             
GSEA_C53S_WT_input = as.data.frame(GSEA_C53S_WT_input)
GSEA_C53S_WT_input$ensembl = row.names(GSEA_C53S_WT_input)
KEGG_C53S_WT = merge(GSEA_C53S_WT_input,
                           mapp_c53s_wt,
                           by.x = "ensembl",
                           by.y = "ensembl_gene_id")

#first remove na
KEGG_C53S_WT = na.omit(KEGG_C53S_WT)
KEGG_C53S_WT_input = KEGG_C53S_WT$GSEA_C53S_WT_input
names(KEGG_C53S_WT_input) = KEGG_C53S_WT$entrezgene_id
KEGG_C53S_WT_input = sort(KEGG_C53S_WT_input, decreasing = TRUE)


ids_c53s_mock = names(GSEA_C53S_MOCK_input)
mapp_c53s_mock = getBM(attributes = c("ensembl_gene_id",
                                    "entrezgene_id",
                                    "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = ids_c53s_mock,
                     mart = ensembl)

GSEA_C53S_MOCK_input = as.data.frame(GSEA_C53S_MOCK_input)
GSEA_C53S_MOCK_input$ensembl = row.names(GSEA_C53S_MOCK_input)
KEGG_C53S_MOCK = merge(GSEA_C53S_MOCK_input,
                     mapp_c53s_mock,
                     by.x = "ensembl",
                     by.y = "ensembl_gene_id")

#first remove na
KEGG_C53S_MOCK = na.omit(KEGG_C53S_MOCK)
KEGG_C53S_MOCK_input = KEGG_C53S_MOCK$GSEA_C53S_MOCK_input
names(KEGG_C53S_MOCK_input) = KEGG_C53S_MOCK$entrezgene_id
KEGG_C53S_MOCK_input = sort(KEGG_C53S_MOCK_input, decreasing = TRUE)


ids_wt_mock = names(GSEA_WT_MOCK_input)
mapp_wt_mock = getBM(attributes = c("ensembl_gene_id",
                                      "entrezgene_id",
                                      "hgnc_symbol"),
                       filters = "ensembl_gene_id",
                       values = ids_wt_mock,
                       mart = ensembl)

GSEA_WT_MOCK_input = as.data.frame(GSEA_WT_MOCK_input)
GSEA_WT_MOCK_input$ensembl = row.names(GSEA_WT_MOCK_input)
KEGG_WT_MOCK = merge(GSEA_WT_MOCK_input,
                       mapp_wt_mock,
                       by.x = "ensembl",
                       by.y = "ensembl_gene_id")

#first remove na
KEGG_WT_MOCK = na.omit(KEGG_WT_MOCK)
KEGG_WT_MOCK_input = KEGG_WT_MOCK$GSEA_WT_MOCK_input
names(KEGG_WT_MOCK_input) = KEGG_WT_MOCK$entrezgene_id
KEGG_WT_MOCK_input = sort(KEGG_WT_MOCK_input, decreasing = TRUE)

head(KEGG_C53S_MOCK_input)
head(KEGG_C53S_WT_input)
head(KEGG_WT_MOCK_input)

###now we have to remove duplicated entries, keeping
#only first occurence
KEGG_C53S_MOCK_input = KEGG_C53S_MOCK_input[!duplicated(names(KEGG_C53S_MOCK_input))]
KEGG_C53S_WT_input = KEGG_C53S_WT_input[!duplicated(names(KEGG_C53S_WT_input))]
KEGG_WT_MOCK_input = KEGG_WT_MOCK_input[!duplicated(names(KEGG_WT_MOCK_input))]
##we will do it accoGSEA_C53S_WT_input##we will do it according to KEGG pathways
KEGG_res_C53S_WT = gseKEGG(geneList = KEGG_C53S_WT_input,
                         organism= "hsa",
                         keyType = "kegg",
                         minGSSize = 5,
                         maxGSSize = 500,
                         pvalueCutoff = 0.05,
                         verbose = TRUE,
                         eps = 0) ##

KEGG_res_C53S_WT = as.data.frame(KEGG_res_C53S_WT)
write.csv(KEGG_res_C53S_WT, "KEGG_C53S_vs_WT.csv",
          row.names = TRUE)

KEGG_res_C53S_MOCK = gseKEGG(geneList = KEGG_C53S_MOCK_input,
                           organism= "hsa",
                           keyType = "kegg",
                           minGSSize = 5,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           verbose = TRUE,
                           eps = 0) ##

KEGG_res_C53S_MOCK = as.data.frame(KEGG_res_C53S_MOCK)
write.csv(KEGG_res_C53S_MOCK, "KEGG_C53S_vs_MOCK.csv",
          row.names = TRUE)


KEGG_res_WT_MOCK = gseKEGG(geneList = KEGG_WT_MOCK_input,
                             organism= "hsa",
                             keyType = "kegg",
                             minGSSize = 5,
                             maxGSSize = 500,
                             pvalueCutoff = 0.05,
                             verbose = TRUE,
                             eps = 0) ##

KEGG_res_WT_MOCK = as.data.frame(KEGG_res_WT_MOCK)
write.csv(KEGG_res_WT_MOCK, "KEGG_WT_vs_MOCK.csv",
          row.names = TRUE)


###let us try to graph

###I will do it using ggplot which is better
library(ggplot2)


dotplot_C53S_WT = GSEA_KEGG_df_to_dotplot(KEGG_res_C53S_WT,
                                    title =  "KEGG C53S vs WT")
pdf("Dotplot_KEGG_C53S_vs_WT.pdf", width = 10, height = 8)
dotplot_C53S_WT
dev.off()


dotplot_C53S_MOCK = GSEA_KEGG_df_to_dotplot(KEGG_res_C53S_MOCK,
                                          title =  "KEGG C53S vs MOCK")
pdf("Dotplot_KEGG_C53S_vs_MOCK.pdf", width = 10, height = 8)
dotplot_C53S_MOCK
dev.off()


dotplot_WT_MOCK = GSEA_KEGG_df_to_dotplot(KEGG_res_WT_MOCK,
                                            title =  "KEGG WT vs MOCK")
pdf("Dotplot_KEGG_WT_vs_MOCK.pdf", width = 10, height = 8)
dotplot_WT_MOCK
dev.off()

### now I want to change the ENTREZ IDS in the KEGG entries
#to gene symbols

###first let me find a mapping
library(biomaRt)

ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("hgnc_symbol",
                                "entrezgene_id"),
                 mart = ensembl)
mapping = na.omit(mapping) #now we have the mapping
#between gene symbols and entrezgene_id

#prepare the inputs
#basically the core-enrichment has the genes in each
#KEGG category. so we have to convert that into the
#gene symbols

#and then from each condition find out which genes
#in that KEGG was from the DEG
#we load the DEG per comparison
DEG_c53s_mock = read.csv("sig_DEG_C53S_vs_MOCK.csv",
                         header = TRUE)
DEG_wt_mock = read.csv("sig_DEG_WT_vs_MOCK.csv",
                     header = TRUE)
DEG_c53s_wt = read.csv("sig_DEG_C53S_vs_WT.csv",
                         header = TRUE)
task_entrez_to_DEG <- function(KEGG_df,
                               mapping,
                               DEG_df){
  #this is a function that takes a KEGG outputs (KEGG_df)
  #then converts the entrez ids per enriched KEGG to 
  #gene symbol according to (mapping)
  #then adds  zthe symbol and then tell us which 
  #one from these symbols were DEG according to the list
  #provided by DEG_df
  
  #first check the inputs
  #the KEGG_df should have the column #core_enrichment
  
  if(!("core_enrichment" %in% colnames(KEGG_df))){
    stop("KEGG_df missing the column \"core_enrichment\"")}
  
  if(!all(c("hgnc_symbol", "entrezgene_id") %in% colnames(mapping))){
    stop("\"mapping\" df missing the columns \"hgnc_symbol\" and/or \"entrezgene_id\"")
  }
  if(!all(c("X", "log2FoldChange") %in% colnames(DEG_df))){
    stop("\"DEG_df\" missing the columns \"X\" and/or \"log2FoldChange\"")
  }
  
 ###now we checked all columns are there
  ## let us start the mapping
  
  #first extract the entrez ids from the  KEGG_df
  entrez_ids = KEGG_df$core_enrichment
  ## we also need the positive and negative DEG
  positive_DEGs = DEG_df$X[DEG_df$log2FoldChange > 0]
  negative_DEGs = DEG_df$X[DEG_df$log2FoldChange < 0]
  ##now we need the following functions
  #for conversion of entrez_ to hgnc
  sep_entrez_to_symbol<- function(entrez_list, mapping){
    #this is a function that takes a list of entries  sepearated by /
    #and we convert each of these entries from entrez to hgnc_symbol
    #using mapping
    
    
    separated_entries = strsplit(entrez_list, "/")
    separated_entries = separated_entries[[1]]
    correspnding_symbols = mapping$hgnc_symbol[match(separated_entries, 
                                                     mapping$entrezgene_id)]
    correspnding_symbols = paste(correspnding_symbols,
                                 collapse = "/")
    return(correspnding_symbols)
  }
  hgnc_symbols = lapply(entrez_ids,
                        function(x) sep_entrez_to_symbol(x,
                                                  mapping))
  hgnc_symbols = unlist(hgnc_symbols)
  KEGG_df$hgnc_symbols = hgnc_symbols
  
  #this is to get the positive and negative DEG
  hgnc_to__DEG <- function(hgnc_entries, DEG){
    #this is a function that takes a vector of hgnc entries
    #separated by "/"
    #and keeps only those that are in the DEG list
    hgnc_entries = unlist(strsplit(hgnc_entries, "/"))
    relevant_entries = hgnc_entries[hgnc_entries %in% DEG]
    relevant_entries = paste(relevant_entries, collapse = "/")
    return(relevant_entries)
    
  }
  
  positive_hgnc = unlist(lapply(hgnc_symbols,
                                function(x) hgnc_to__DEG(x,
                                    positive_DEGs)))
  positive_hgnc[positive_hgnc == ""] = NA
  
  negative_hgnc = unlist(lapply(hgnc_symbols,
                                function(x) hgnc_to__DEG(x,
                                                         negative_DEGs)))
  negative_hgnc[negative_hgnc == ""] = NA
  KEGG_df$positive_hgnc = positive_hgnc
  KEGG_df$negative_hgnc = negative_hgnc
  return(KEGG_df)
  }

new_KEGG_c53s_vs_mock = task_entrez_to_DEG(KEGG_res_C53S_MOCK,
                                           mapping,
                                           DEG_c53s_mock)
write.csv(new_KEGG_c53s_vs_mock, "KEGG_C53S_vs_MOCK.csv",
          row.names = TRUE)

new_KEGG_c53s_vs_wt = task_entrez_to_DEG(KEGG_res_C53S_WT,
                                           mapping,
                                           DEG_c53s_wt)
write.csv(new_KEGG_c53s_vs_wt, "KEGG_C53S_vs_WT.csv",
          row.names = TRUE)

new_KEGG_wt_vs_mock = task_entrez_to_DEG(KEGG_res_WT_MOCK,
                                           mapping,
                                           DEG_wt_mock)
write.csv(new_KEGG_wt_vs_mock, "KEGG_WT_vs_MOCK.csv",
          row.names = TRUE)


