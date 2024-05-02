# LatinCells pipeline 

**Usage:**
~~~
nextflow {Path_to_script}/main.nf --dataDir="{Path_to_Data}/RawData/*" --Project_Name="Chile_03" -c {Path_to_script}/nextflow.config -resume --NPlex=4
~~~
##### *Parameters*
###### -Mandatories
~~~
--dataDir = Sets the input sample table file.
-c = Sets the configuration file.
~~~
###### -Alternatives
~~~
--Project_Name = Sets a custom directory for the project
--NPlex = Sets the number of samples per run in multiplexed libraries.
--help = Prints this help message.
-resume = Continue an analysis after the last well executed part of a workflow, before a error or a change made in thee code of a process.[Nextflow parameter]
~~~
## Bioinformatics Tools
## Cellranger
Description: Cell Ranger is a set of analysis pipelines developed by 10x Genomics for single-cell RNA-seq data. It performs sample demultiplexing, barcode processing, gene counting, and V(D)J transcript sequence assembly.
GitHub Repository: [Cellranger](https://github.com/10XGenomics/cellranger)
## Freemuxlet
Description: Freemuxlet is a tool for demultiplexing and detecting sample multiplets in single-cell RNA-seq data.
GitHub Repository: [Popscle: Demuxlet/Freemuxlet](https://github.com/statgen/popscle)
## SoupX
Description: SoupX is an R package for estimating and removing cell-free mRNA contamination in droplet-based single-cell RNA-seq data.
GitHub Repository: [SoupX](https://github.com/immunogenomics/harmony)
## scDblFinder
Description: scDblFinder detects and handles doublets/multiplets in single-cell sequencing data.
GitHub Repository: [scDblFinder](https://github.com/plger/scDblFinder)
## Seurat
Description: Seurat is an R package for single-cell RNA-seq data analysis, including clustering, dimensionality reduction, and visualization.
GitHub Repository: [Seurat](https://github.com/satijalab/seurat)
## Harmony
Description: Harmony is an algorithm for batch correction and integration of single-cell RNA-seq data.
GitHub Repository: [Harmony](https://github.com/immunogenomics/harmony)
## Cellchat
Description: Cellchat is an R package for cell-cell communication analysis in single-cell RNA-seq data.
GitHub Repository: [Cellchat](https://github.com/jinworks/CellChat)
## Monocle 3
Description: Monocle 3 is an analysis toolkit for single-cell RNA-Seq experiments. To use this package, you will need the R statistical computing environment (version 3.0 or later) and several packages available through Bioconductor and CRAN.
GitHub Repository: [Monocle 3](https://github.com/cole-trapnell-lab/monocle3)

For more details on each tool, installation instructions, and usage examples, please refer to the respective GitHub repositories. Happy exploring!