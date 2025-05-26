here is what is in each file

KW_MIA40C53S.counts.txt
the count data after using featurecounts (raw counts)

summarized_count_data.csv
count_data basically after arranging the samples according to their order and naming them properly (all genes are featured)

valid_count_data.csv
count data with sample names after keeping nly the genes that pass our criteria, which is the gene should have half of the samples at least having reads of more than 10. we are using this as an input for now

rpkm.csv
this is the calculated rpkm, which might be used for the PCA analysis

c53s_vs_mock_log2fc_raw.csv
here are the raw log2fc results from DESeq2 , where the comparison
is between C53S and mock across different time points

c53s_vs_wt_log2fc_raw.csv
here are the raw log2fc results from DESeq2 , where the comparison
is between C53S and wt across different time points

wt_vs_mock_log2fc_raw.csv
here are the raw log2fc results from DESeq2 , where the comparison
is between wt and mock across different time points

sig_DEG_C53S_vs_Mock.csv
the DEG between C53S and mock after DESeq2 analysis and using a filter
of p.adj < 0.05 and abs(log2fc) >= 1

sig_DEG_C53S_vs_WT.csv
the DEG between C53S and WT after DESeq2 analysis and using a filter
of p.adj < 0.05 and abs(log2fc) >= 1

sig_DEG_WT_vs_Mock.csv
the DEG between WT and mock after DESeq2 analysis and using a filter
of p.adj < 0.05 and abs(log2fc) >= 1

ensembl_to_gene_conversion
the one to one mapping of ensemble ids into gene names 

GSEA_results_C53S_vs_MOCK
results of GSEA analysis performed on the log2fc comparison between
C53S and MOCK using ranked gene lists based on the stats column from 
DESEQ results

GSEA_results_C53S_vs_WT
results of GSEA analysis performed on the log2fc comparison between
C53S and WT using ranked gene lists based on the stats column from 
DESEQ results

GSEA_results_WT_vs_MOCK
results of GSEA analysis performed on the log2fc comparison between
WT and MOCK using ranked gene lists based on the stats column from 
DESEQ results

KEGG_C53S_vs_MOCK
results of KEGG pathway enrichment performed on log2fc comparison between C53S and MOCK based on ranked gene list using stats from the deseq results

KEGG_C53S_vs_WT
results of KEGG pathway enrichment performed on log2fc comparison between C53S and WT based on ranked gene list using stats from the deseq results

KEGG_WT_vs_MOCK
results of KEGG pathway enrichment performed on log2fc comparison between WT and Mock based on ranked gene list using stats from the deseq results

DESeq2_results_C53S_24h_vs_0h.csv
this is the DESeq2 results outputting log2fc comparison between C53S at 24h timepoint and 0h timepoints using a model that takes into consideration the time and different conditions. 

DESeq2_results_WT_24h_vs_0h.csv
this is the DESeq2 results outputting log2fc comparison between WT at 24h timepoint and 0h timepoints using a model that takes into consideration the time and different conditions. 


DESeq2_results_MOCK_24h_vs_0h.csv
this is the DESeq2 results outputting log2fc comparison between MOCK at 24h timepoint and 0h timepoints using a model that takes into consideration the time and different conditions. 


DEG_DESeq2_C53S_24h_vs_0h.csv
these are the DEG from the comparison using a criteria of 
abs(log2fc) more than 0.5 and p.adj less than 0.05. we used 0.5 instead of 1 for abs(log2fc) because 1 only showed one DEG which is the MIA40 protein itself

KEGG_DESeq2_C53S_24h_vs_0h
This is the KEGG of the results of the DESeq2 comparison between C53S 24h vs 0h taking into consideration all timepoints and condition. Note that the DEG here are from the list above directly with a threshold of 0.5 for abs log2fc

DESeq2_results_C53S_24h_vs_0h_vs_WT
so here is the results of DESEQ analaysis using ONLY WT and C53S data. Here we did a contract between time_24h_vs_0h and ConditionC53S.time24h not only focusing on the time part according to the tutorial that was suggested by 
https://www.atakanekiz.com/technical/a-guide-to-designs-and-contrasts-in-DESeq2/#two-factors-with-interaction
