# PeacoQC
Peak-based selection of high quality cytometry data

## Introduction
The PeacoQC package provides quality control functions that will check 
for monotonic increasing channels and that will remove outliers and unstable 
events introduced due to e.g. clogs, speed changes etc. during the measurement 
of your sample. It also provides the functionality of visualising the quality
control result of only one sample and the visualisation of the results of 
multiple samples in one experiment.

## Installation
You can install this package using the devtools library.

```{r}
BiocManager::install("ComplexHeatmap")
devtools::install_github("saeyslab/PeacoQC")
```
