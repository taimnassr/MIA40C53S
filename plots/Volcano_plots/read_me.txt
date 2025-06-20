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

Volcano_plot_C53S_vs_WT_24h_vs_0h_p_adj_0.2
the same plot as directly above but here we changed the p.adj < 0.2, abs(log2fc) > 0.5

Volcano_plot_C53S_vs_WT_24h_vs_0h_no_log2fc_threshold
same plot as the "Volcano_plot_C53S_vs_WT_24h_vs_0h" but without any log2fc threshold, just p.adj < 0.05

Volcano_Plot_LIMMA_C53S_vs_WT_24h_vs_0h
this is plotting the LIMMA results from comparing
C53S to WT in the 24h vs 0h change 
abs(log2fc) > 0.5 and p.adj < 0.05

Volcano_Plot_DESeq2_C53S_vs_Mock_24h_vs_0h
This is the volcano plot for the log2fc between C53S and Mock comparing 24h to 0h (so what specific changes happens in C53S between 24h and 0h that doesn't happen in Mock). The thresholds put here is abs.log2fc more than 0.5 and p.adj < 0.05.

Volcano_Plot_DESeq2_Non_Mock_vs_Mock_24h_vs_0h
This is similar to the one directly before but here WT and C53S were clustered together into one condition (called non-mock) whose specific changes between 24h and 0h were compared to that in mock with the blue lines showing a threshold of abs log2fc more than 0.5 and p.adj < 0.05