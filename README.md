# HCA-reproducibility
A case study on the detailed reproducibility of a human cell atlas project

## Data description
- **Data.RData**: data from the original paper, including:
	- **discovery.dataset** : discovery dataset with 1140 cells
	- **validation.dataset** : validation dataset with 1244 cells
	- **all.label** : refined cell type for 1077 out of the 1140 cells in the discovery datasets
	- **given.markers** : Maker genes for DC1--DC6

- **Reproduced.dc.markers.RData** : marker genes identified using data for DCs.
- **Reproduced.mono.dc.markers.RData** : marker genes identified using data for both DCs and monocytes.

## R packages
The most of reprodcution work were done using the Package `Seurat (2.3.4)`. There are also some other helpful packages. Run the following code to install all these packages.

```
install.packages(c("Seurat","plyr", 'dplyr', "reshape2", "ggpubr", "VennDiagram"))
```

## Reference
- A. C. Villani et al., "Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors," Science, vol. 356, no. 6335, Apr 21 2017.
- A. Butler, P. Hoffman, P. Smibert, E. Papalexi, and R. Satija, "Integrating single-cell transcriptomic data across different conditions, technologies, and species," Nat Biotechnol, vol. 36, no. 5, pp. 411-420, Jun 2018.

