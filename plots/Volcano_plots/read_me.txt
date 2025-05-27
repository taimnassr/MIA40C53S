here you have the volcano plots after the standard
DESeq2 analysis of the valid_count_data.csv data

you have the following pdfs
C53S_vs_Mock
C53S_vs_WT
WT_vs_Mock

you also have

Volcano_plot_C53S_24h_vs_0h_with_independent_filtering
for C53S comparing 24h vs 0h with the results function of dds using independent filtering for p adj value (independent filtering = TRUE) which is the default

Volcano_plot_C53S_24h_vs_0h_with_independent_filtering
you have the same plot but with (independent filtering = FALSE) since there are so many p.adj = NA. and we want to turn off the independent filtering to get their values

Volcano_plot_Mock_24h_vs_0h
Volcano_plot_C53S_24h_vs_0h
Volcano_plot_WT_24h_vs_0h
These are the volcano plots from the relative comparisons using DESeq2 between 24h and 0h points (independent_filtering = TRUE). taking into account all different conditions and timepoints and building a model only on that comparison
abs(log2fc) > 1 and padj< 0.05

Volcano_plot_C53S_vs_WT_24h_vs_0h
So this comparing C53S to WT (24h vs 0h) so this has two comparisons and it seeing whether C53S change between 24h and 0h would be different than WT (so taking two factors) and the filters here
are
abs(log2fc) > 0.5, and p.adj < 0.05