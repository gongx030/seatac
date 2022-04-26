
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/gongx030/seatac/workflows/R-CMD-check/badge.svg)](https://github.com/gongx030/seatac/actions)
  [![R](https://github.com/gongx030/seatac/actions/workflows/r.yml/badge.svg)](https://github.com/gongx030/seatac/actions/workflows/r.yml)
  <!-- badges: end -->

## SeATAC

SeATAC: a tool for exploring the chromatin landscape and the role of pioneer factors

SeATAC is an R package to estimate the genomic regions with statistically differential chromatin accessibility from multiple ATAC-seq data. Using SeATAC, each genomic region is represented as a V-plot, a dot-plot showing how sequencing reads with different fragment sizes distribute surrounding one or a set of genomic region(s).  For a more detailed overview of the method, please see [the manuscript](https://www.biorxiv.org/content/10.1101/2022.04.25.489439v1).  

![Screen Shot 2022-04-20 at 3 00 15 PM](https://user-images.githubusercontent.com/16770031/164312915-c6239042-37d4-4fc9-967f-b548261ac32a.png)

SeATAC uses a conditional variational autoencoder (CVAE) model to learn the latent representation of the ATAC-seq V-plot.  With the probabilistic representation of the data, we developed a Bayesian method to evaluate the statistical difference between multiple V-plots.  

![Screen Shot 2022-04-20 at 3 00 44 PM](https://user-images.githubusercontent.com/16770031/164312993-0253633e-97cf-4f65-91b8-6595072f642b.png)

SeATAC has significantly better performance on four separate tasks compared to MACS2 and/or NucleoATAC on both synthetic and real ATAC-seq datasets, including (1) detection of differential V-plots; (2) definition of nucleosome positions; (3) detection of nucleosome changes and (4) designation of transcriptional factor binding sites (TFBS) with differential chromatin accessibility. 

![Screen Shot 2022-04-20 at 3 01 13 PM](https://user-images.githubusercontent.com/16770031/164313080-c8bf570c-91c2-4d2e-905c-d7ba63df6a7e.png)
