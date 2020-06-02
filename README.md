# RNA-seq-DE-analysis
Rshiny application for differential expression analysis

This project allows to validate our first year of Master's degree in Bioinformatics. 

It consists in the creation of a Shiny application allowing the RNA sequencing differential expression analysis using Rstudio and the DESeq2 package.

The operation of the application is explained on the first page of the application.


Introduction
This is an R Shiny web interactive application developed as part of a course project. The purpose of this application is to perform a differential expression analysis from a counts table in order to help researchers getting interpretable results.
This application uses the package DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) from Bioconductor. It is a package to study differential gene expression analysis based on the negative binomial distribution. It allows a quick visualization of results in order to analyze the counts table data. The results will be given in different forms like graphics, heatmaps, MAplot or even Volcano plot.
1. Upload data
The input data files accepted for this App are 3 files in '.txt', '.csv' or '.tsv' format separated by comma, tabulation or semi-colon. This App necessarily requires a 'Count Data Table' and a 'Metadata Table'. An optional 'Annotation File' can be added
1.1 Count Data Table
The Count Data Table must contain the count for each sample of the experiment for each gene and the first column must be gene ID or gene name as below :
 GeneID Sample1
Sample2
486.00 0.00 523.00 258.00 81.00 0.00
Sample3
904.00 0.00 616.00 364.00 73.00 1.00
Sample4
445.00 0.00 371.00 237.00 66.00 0.00
Sample5
1170.00 0.00 582.00 318.00 118.00 2.00
Sample6
1097.00 0.00 781.00 447.00 94.00 0.00
       ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 ENSG00000000938
1.2 Metadata Table
723.00 0.00 467.00 347.00 96.00 0.00
                                   The Metadata table must contain the information
one corresponds to the samples in the same order as the columns of the Count Table. The second one is a condition column. You can add as many columns as you have factors in your experiment.
sample
sample1 sample2 sample3 sample4 sample5 sample6
condition ...
control treated control treated control treated
of the experiment with at least 2 columns. The first
                      https://gabin-coudray.shinyapps.io/m1project/ 1/3
02/06/2020 RNA-seq DE analysis
sample
sample7 sample8
condition ...
control treated
         1.2 Annotation File
The Annotation File contains informations about the genes. If you have one, it must contain a column named 'symbol' in which we can find the symbol of each gene.
GeneID entrez
symbol
chr start
X 100627109
X 100584802
20 50934867
1 169849631
1 169662007
1 27612064
end
100639991
100599885
50958555
169894267
169854080
27635277
strand bio
-1 prot
1 prot
-1 prot
-1 prot
1 prot
-1 prot
        ENSG00000000003
ENSG00000000005
ENSG00000000419
ENSG00000000457
ENSG00000000460
ENSG00000000938
7105 TSPAN6
64102 TNMD
8813 DPM1
57147 SCYL3
55732 C1orf112
2268 FGR
                                        2. Results
The results will be display after running DESeq2. You will obtain 9 differents results :
- Count distribution - Count by gene
- Depth of sample - Dispersion
- PCA
- MA plot
- Volcano plot
- Sample distance matrix
- Gene expression Heatmap
You can download all the results plots at the bottom of all these pages.
