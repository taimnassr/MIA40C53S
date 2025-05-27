##here you have the functions used in the analysis
## to create volcano plots from DESEQ data results
volcano_plots = function(df, title,
                         log2fc_threshold = 1,
                         padj_threshold = 0.05){
  
  #df should be containing DESEQ log2fc results as dataframe
  #there should be a column with the log2foldchange
  #called log2FoldChange  and a column called padj
  ## also there should be a column called significance
  # to check whether a gene is considered DEG or not
  
  #also the rownames should contain the gene names
  col_names = c("log2FoldChange", "padj")
  if(!all( col_names %in% 
           colnames(df))){
    stop("Missing columns: log2FoldChange or padj")
  }
  
  #conditon name is like mock vs wt 
  
  #sig_threshold is what value should the abs log2fc  be higher
  #than. usually it is 1 unless you used another sig_threshold
  df$significance = "Non-Significant"
  df$significance[abs(df$log2FoldChange) > log2fc_threshold &
                    df$padj < padj_threshold] = "Significant"
  
  df$padj[df$padj == 0] = 1e-100
  sig_df = df[df$significance == "Significant",]#this is the 
  #significant df
  sig_df$gene = rownames(sig_df)
  library(ggplot2)
  library(ggrepel)
  title = paste ("DESeq2 Results", title )
  
  plot = ggplot(df, 
                aes(x = log2FoldChange,
                    y = -log10(padj))) +
    geom_point(aes(colour = significance)) +
    scale_color_manual(values = c("Significant" = "red",
                                  "Non-Significant" = "grey"))+
    labs(title = title,
         x = "Log2FC",
         y = "-log10(P. adjusted)")+
    geom_text_repel(data = sig_df,
                mapping = aes(x = log2FoldChange,
                        y = -log10(padj),
                        label = gene),
                    size = 3,
                    max.overlaps = 30)+
    geom_vline(xintercept = c(-log2fc_threshold,
                              log2fc_threshold), color = "blue")+
    geom_hline(yintercept = -log10(padj_threshold),
               ,color = "blue") +
    ylim(c(0,200)) + 
    scale_x_continuous(breaks = seq(-5,5, by  = 1),
                       labels = seq(-5,5, by = 1),
                       limits = c(-5,5)) 
  return(plot)
}

##to create dotplot for GSEA visualization
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
          y = "KEGG Pathways",
          title = title,
          size = "Number of genes per KEGG Entry",
          color = "Adjusted P.value")+ theme_minimal()+
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 14)) +
    scale_color_gradient(low = "blue",
                         high = "red",
                         limits = c(0,0.05)) 
  
  return(dotplot)
}

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
