#first setting your work_directory
setwd("/mnt/d/Data_Analysis/NGS_KW_MIA40C53S/count_data")#where the data is stored

#load in the data
count_data = read.table("KW_MIA40C53S.counts.txt", header = T,
                        sep = "\t")

#keep the relevant info, which is gene_id, gene_length, and counts

count_data = count_data[,-c(2,3,4,5)]# remove the columns of the
#entries we don't want

#now we want to change the arrangment of the columns and the names
#so it is easier to work with

library(tidyverse)#load this library

#first change the names of some columns to only the sample names
col_names = names(count_data)
col_names[3:length(col_names)] = str_extract(col_names[3:length(col_names)],
                                             "S\\d+")# this will keep
names(count_data) = col_names#changed the names of columns

#S92- S127 from the big sample names

new_col_names= c(col_names[1], col_names[2], paste0("S",92:127)) # generate new order

#now we want to arrange the column names according to the new order

#to do that change the data into a dataframe
count_data_df= as.data.frame(count_data)
count_data_df2= count_data_df[,new_col_names]
#now we have count_data df2 with simplified and ordered sample names

###so count_data_df2 we work with

##now we try to match the samples to their names.
#i will open the sample_names file
sample_names =read.delim("Sample_Names.tab", header = TRUE)#
sample_names = sample_names$CCG.Sample.ID....Sample.Name
#change sample names to coloumn names
#now align sample names to the column names of the count_data_df2

#since we already arranged them it is very simple
names(count_data_df2)[3:ncol(count_data_df2)] = sample_names

write.csv(count_data_df2, "summarized_count_data.csv", row.names = FALSE)

#now we want to filter the data a bit, where we filter the genes based on
#at least 10 reads in 18/36 in the samples

gene_valid_10reads <- apply(count_data_df2[,3:ncol(count_data_df2)],
                            1, FUN = function(gene_data) {
  sum(gene_data >= 10) >= 18
})
#check how many genes pass the criteria
length(gene_valid_10reads[gene_valid_10reads == TRUE]) #18898 genes
#create dataframe with only these genes
count_data_valid = count_data_df2[gene_valid_10reads,]

#same the valid count_data which we would use
write.csv(count_data_valid, "valid_count_data.csv", row.names = FALSE)

#instead of doing rpkm, we need to do the deseq data analysis
#we use the vst, variance stable transformation to plot the PCA

library(DESeq2)
library(stringr)

count_matrix # we need a count matrix where rows are genes and columns
#are samples. this is already have done but if not should be done

colnames = colnames(count_matrix)#get the column names
colnames
sample_names = colnames
sample_type = sub("-.*", "", sample_names) # get sample type
sample_time = str_extract(sample_names, "\\d{1,2}(?=h)")

col_data = data.frame(
  sample_type = sample_type,
  time_point = sample_time,
  row.names = sample_names
)# important for deseq analysis

all(rownames(col_data) == colnames(count_matrix)) # check if true

dds_pca <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~1)

#now do the transformation needed for PCA
vsd = vst(dds_pca)
vsd_matrix = assay(vsd) # here you have the data that you use for PCA

pca = prcomp(t(vsd_matrix))
pca_input = as.data.frame(pca$x)
pca_input$PC1
rownames(pca_input)
library(ggplot2)
install.packages("ggrepel")#to seperate points from each other
library(ggrepel)

pca_plot = ggplot(pca_input, aes(x= PC1,y = PC3,
                                 label = rownames(pca_input))) +
  geom_point() + geom_text_repel(max.overlaps = Inf) +
  labs(x = "PC1", y = "PC5") +
  theme(axis.line= element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 20))+
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black")+
  theme_classic()
        
#save as pdf
pdf("PCA_plot_vst_PC1_PC5.pdf", width = 8, height = 6)
print(pca_plot)
dev.off()





