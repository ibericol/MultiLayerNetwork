# MultiLayerNetwork
MultiLayerNetwork package for R

This R package performs the Multilayer Network analysis from Yi Guan1, Chia Hsin Cheng, Luis Iberico, Sriman Narain, Sherman Bigornia, Mahdi O. Garelnabi, Tammy Scott, José M. Ordovás, Katherine L. Tucker, Rafeeque Bhadelia, and Bang-bon Koo.

Several functions are available in order to estimate the networks and analyse the key paths between nodes. In addition, you can compare different groups by doing the non parametric Mann-Whitney U Test on the bootstrapped weights of the resulting matrix. 

The package uses ppcor to estimate the partial correlations, WGCNA to perform the TOM, and cluster to carry out the k-medoids clustering. In addition, dplyr, purrr, tidyr, and ggplot2 are also required, so you may find tyverse useful.

