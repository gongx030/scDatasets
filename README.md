# scDatasets

Processed single cell sequencing datasets

## 1. Installation

The recommended installation method for `scDatasets` is using `install_github` command from `devtools` library.  You will first have to have [devtools](https://github.com/hadley/devtools) package installed.

```r
library(devtools)
install_github('gongx030/scDatasets')
```

A number of needed packages are installed in this process.

## 2. Datasets

The scRNA-seq dataset could be loaded in this way:
```r
library(scDatasets)
library(SummarizedExperiment)
data(pollen)
pollen
```

Here summarizes the name, data type, expression unit and source of included datasets:

| Name | Data_Type | Expression_Unit | Source |
| :---: | :---: | :---: | :--- |
| guo | scPCR | exp(-Ct) | Guo et al. Resolution of cell fate decisions revealed by single-cell gene expression analysis from zygote to blastocyst. Dev Cell. 2010 Apr 20;18(4):675-85. ([20412781](https://www.ncbi.nlm.nih.gov/pubmed/20412781)). |
| guo2 | scRNA-seq | count | Guo et al. The Transcriptome and DNA Methylome Landscapes of Human Primordial Germ Cells. Cell. 2015 Jun 4;161(6):1437-52. ([26046443](https://www.ncbi.nlm.nih.gov/pubmed/26046443)).  |
| loh | scRNA-seq | count | Loh et al. Mapping the Pairwise Choices Leading from Pluripotency to Human Bone, Heart, and Other Mesoderm Cell Types. Cell. 2016 Jul 14;166(2):451-467. ([27419872](https://www.ncbi.nlm.nih.gov/pubmed/27419872)) |
| pollen | scRNA-seq | count | Pollen et al. Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity and activated signaling pathways in developing cerebral cortex. Nat Biotechnol. 2014 Oct;32(10):1053-8. ([25086649](https://www.ncbi.nlm.nih.gov/pubmed/25086649)).  |
| trapnell | scRNA-seq | FPKM | Trapnell et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat Biotechnol. 2014 Apr;32(4):381-386. ([24658644](https://www.ncbi.nlm.nih.gov/pubmed/24658644)). |
| usoskin | scRNA-seq | FPKM | Usoskin et al. Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nat Neurosci. 2015 Jan;18(1):145-53. ([25420068](https://www.ncbi.nlm.nih.gov/pubmed/25420068)). |
| blakeley | scRNA-seq | count | Blakeley et al. Defining the three cell lineages of the human blastocyst by single-cell RNA-seq. Development. 2015 Oct 15;142(20):3613.  ([26487783](https://www.ncbi.nlm.nih.gov/pubmed/26487783)). |
| iPSC | scRNA-seq | count | Single-Cell RNA-seq of iPSC-CM differentiations at iPSC stage (D0), differentiation day 6 (D6), D10, D30 and D60. |
| petropoulos | scRNA-seq | count | Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human
Preimplantation Embryos. Cell. 2016 May 5;165(4):1012-26. |
| deng | scRNA-seq | count | Deng et al. Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in
mammalian cells. Science. 2014 Jan 10;343(6167):193-6. |

