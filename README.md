
<h1 align="center">
  <br>
  <img src="https://github.com/MargaridaGoncalves/Shiny-App/blob/main/GeneSEA%20Explorer/official_banner.png" width="2000"></a>
  <br>

  
# GeneSEA Explorer
  
The GeneSEA Explorer is an innovative bioinformatics tool for conducting differential gene expression, using SEA. Its user-friendly interface enables users to effortlessly explore
and analyze diverse DEGs outputs and to conduct functional enrichment analysis. The innovative use of Shannon entropy to aggregate all DEGs outputs provides an informative selection of DEGs, enhancing researchers’ understanding of their RNA-Seq data results and minimizing challenges for researchers less confident in this domain.


More info on the usage of the application and some examples are detailed in the [GeneSEA Explorer User's Guide](https://github.com/MargaridaGoncalves/GeneSEA-Explorer/blob/main/GeneSEA%20Explorer%20Users%20Guide.pdf)

## Set Up

### R Session


> [!WARNING]  
> GeneSEA Explorer currently only supports version 4.3.3 of R. The current R version 4.4.1 is still not supported.

> [!IMPORTANT]  
> In RStudio open the Session in the GeneSEA Explorer folder.  

> For that go to: Session > Set Working Directory > GeneSEA Explorer



### Package Installation

Introduce in the Rstudio Console the following:

```r
# Package names
packages <- c("plotly", "DT", "shiny", "ggplot2", "dplyr", "readxl","R.utils","readr","bslib", "utils", "DOSE", "combinat",
              "tidyverse", "inops","densityClust","wesanderson","paletteer", "devtools", "shinyBS","gprofiler2",
              "ggupset", "shinyWidgets","shinythemes","GGally","ggbump", "genekitr", "ggVennDiagram","ggvenn","RColorBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



install.packages(c("plotly", "DT", "shiny", "ggplot2", "dplyr", "readxl","R.utils","readr","bslib", "utils", "DOSE","waiter",
                   "tidyverse", "inops","densityClust","wesanderson","paletteer", "devtools", "shinyBS","gprofiler2",
                   "ggupset", "shinyWidgets","shinythemes","GGally","ggbump", "genekitr", "ggVennDiagram","ggvenn","RColorBrewer","extrafont"))



# Package names
BiocManagerpackages <- c("clusterProfiler","enrichplot","DESeq2","edgeR","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db",
                         "org.Mmu.eg.db","AnnotationDbi","clustifyr","Biobase","GOSemSim","vidger","DEGreport","preprocessCore",
                         "sva","affydata","Glimma","limma","BiocGenerics","enrichplot","pathview","topGO","vidger","Biobase")

# Install packages not yet installed
installed_packages <- BiocManagerpackages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(BiocManagerpackages[!installed_packages])
}
```

```r
# Install PoissonSeq Package
if("PoissonSeq" %in% rownames(installed.packages()) == FALSE) {install.packages("PoissonSeq_1.1.2.tar.gz", repos = NULL, type ="source", dependencies = TRUE)}
```

### Packages loading

```r
invisible(lapply(packages, library, character.only = TRUE))

invisible(lapply(BiocManagerpackages, library, character.only = TRUE))
```


## Team 
- Ana M. Gonçalves [ORCiD](https://orcid.org/0009-0001-0800-0019)
- Pedro Macedo [ORCiD](https://orcid.org/0000-0002-4371-8069)
- Patrício Costa [ORCiD](https://orcid.org/0000-0002-1201-9177)
- Nuno S. Osório's [ORCiD](https://orcid.org/0000-0003-0949-5399) [GitHub](https://github.com/nunososorio)

## Author contributions 
Ana M. Gonçalves was the main developer, with support from the remaining authors. All authors revised the user's guide.







