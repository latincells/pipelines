# LatinCells pipeline 

**Usage:**
~~~
nextflow {Path_to_script}/main.nf --dataDir="{Path_to_Data}/RawData/*" --Project_Name="Chile_03" -c {Path_to_script}/nextflow.config --NPlex=4 -resume ...more parameters
# with params file:
nextflow {Path_to_script}/main.nf -params-file Tests.yml -profile NLHPC -resume
~~~
##### *Parameters*
###### -Mandatories
~~~
--dataDir = Sets the input sample table file.
--ref_data = Directory with the organism genome processed by cellranger tools.
-c = Sets the configuration file.
~~~
###### -Alternatives
~~~
--Project_Name = Sets a custom directory for the project. (def. "General")
--Sample_metadata = File with metadata of samples, must be a Tabular delimited file (TSV|TAB), and must have a column calles "Sample" matching the names of raw data folders.
--Batch = For Seurat analysis, column with the batches. (def. "orig.ident")
--Freemuxlet = Set use of Freemuxlet to true for demultiplexing without genotype. (def. false)
--ModelsCelltypist = Models for use of celltypist, . (def. "", no celltypist annotation)
--VCF_Files = Path to libraries VCFs, to use for genotyping guided demultiplexing. (def. "", no genotyping guided demultiplexing)
--NPlex = Sets the number of samples per run in multiplexed libraries.
--MainOutdir = Name of directory to save outputs. (def. "LatinCells_Results/${params.Project_Name}/")
-resume = Continue an analysis after the last well executed part of a workflow, before a error or a change made in thee code of a process.[Nextflow parameter].
--help = Prints this help message.
~~~

---

# LatinCells - Data Freeze and Integration

This repository contains the raw and processed data associated with the LatinCells project, focusing on single-cell transcriptomics and genotyping data across multiple Latin American populations.

## Directory Structure

```
02.Rawdata/
├── 01.scRNAseq
│   ├── 01.PBMC
│   │   ├── LC-BR-Lib01/
│   │   ├── LC-BR-Lib02/
│   │   ├── LC-BR-Lib03/
│   │   ├── LC-BR-Lib04/
│   │   ├── LC-BR-Lib05/
│   │   ├── LC-BRMX-Lib00/
│   │   ├── LC-CL-Lib03/
│   │   ├── LC-CL-Lib04/
│   │   ├── LC-CL-Lib05/
│   │   ├── LC-CL-Lib06/
│   │   ├── LC-MX-Lib01/
│   │   ├── LC-MX-Lib02/
│   │   ├── LC-MX-Lib03/
│   │   ├── LC-MX-Lib04/
│   │   ├── LC-MX-Lib05/
│   │   ├── LC-MX-Lib06/
│   │   └── PreviusIDs.txt
│   └── 02.Gallbladder/
├── 02.Genotypes
│   ├── 01.Array_GSA_V3/
│   ├── 02.Genomes_Fastq_vcf/
│   │   └── PLINKabril25v38.vcf
│   ├── LatinCells_PopGen/
│   │   ├── Analysis/
│   │   ├── Data/
│   │   ├── Figures/
│   │   ├── README.md
│   │   ├── Reports/
│   │   └── Scripts/
│   └── popgen/
│       ├── analysis/
│       ├── data/
│       ├── figures/
│       ├── test.genome
│       ├── test.hh
│       ├── test.log
│       ├── test.nosex
│       ├── test.prune.in
│       └── test.prune.out
└── 03.All_merged_VCF/
```

## Description

- `01.scRNAseq/`: Contains raw single-cell RNA sequencing data.
  - `01.PBMC/`: Peripheral Blood Mononuclear Cell libraries from Brazil, Mexico, Colombia, Ecuador, Peru, Uruguay, USA and Chile.
  - `02.Gallbladder/`: Placeholder for gallbladder-related datasets (in development).
- `02.Genotypes/`: Contains genotype data from arrays and sequencing.
  - `01.Array_GSA_V3/`: Genotyping array raw files (GSA v3).
  - `02.Genomes_Fastq_vcf/`: Whole genome FASTQ and merged VCF files.
  - `LatinCells_PopGen/`: Population genomics analyses, scripts, and figures.
  - `popgen/`: Legacy or preliminary analyses using PLINK.
- `03.All_merged_VCF/`: Final merged genotype VCF files used for downstream analyses.

## Notes

- All libraries are organized by sequencing pool.
- Population genomics scripts and figures are self-contained within their directory
- For reproducibility, please refer to the README files in subdirectories when available.

---


---

# LatinCells - Data Freeze v1 (2025-03-18)

This repository contains the first stable version of the data freeze for the LatinCells project. It includes processed single-cell RNA sequencing data and integration-ready objects following quality control and demultiplexing steps.

## Directory Structure

```
LatinCells_DataFreeze_v1_20250318/
├── 01.scRNAseq
│   ├── 01.Cellranger/
│   ├── 02.Demultiplexing/
│   ├── 03.Fastq_per_individual/
│   ├── 04.SoupX_per_library/
│   ├── 05.QC_doublets_per_library/
│   ├── 06.QC_per_Chip10x/
│   └── 07.Data_Freeze_Integration/
│       ├── 01.QC_Seurat/
│       ├── 02.Pre-integration/
│       ├── 03.Integration/
│       ├── 04.Post-integration/
│       └── 05.Seurat_Cell_Annotation/
└── README.txt
```

## Description of Contents

- `01.Cellranger/`: Outputs from Cell Ranger pipelines (count, multi, etc.).
- `02.Demultiplexing/`: Demux results and metadata to assign cells to individuals.
- `03.Fastq_per_individual/`: Demultiplexed FASTQ files per individual.
- `04.SoupX_per_library/`: Ambient RNA correction (SoupX) results per library.
- `05.QC_doublets_per_library/`: Doublet detection reports and filtered objects.
- `06.QC_per_Chip10x/`: Quality control metrics grouped by 10x chip.
- `07.Data_Freeze_Integration/`: Final processed datasets used for integration.
  - `01.QC_Seurat/`: Quality-controlled Seurat objects per sample.
  - `02.Pre-integration/`: Normalized and batch-corrected data before integration.
  - `03.Integration/`: Integrated Seurat objects across donors and sites.
  - `04.Post-integration/`: Post-integration cleanup and clustering.
  - `05.Seurat_Cell_Annotation/`: Final annotations of cell types and metadata.

## Notes

- The data here are ready for downstream analyses and represent a consistent freeze of processed single-cell data across all collected samples.
- For detailed metadata, see `README.txt` and subdirectory documentation.



---

### Help File Headers – Explained

These headers appear at the top of metadata or documentation files (e.g., `.txt` or `.tsv`) and are meant to provide clarity, traceability, and reproducibility for collaborators and analysts. Based on your example (`PreviusIDs.txt`), here’s what each line means:

---

- `# File:`  
  **The name of the file.**  
  Helps track which file is being documented.  
  Example: `PreviusIDs.txt`

- `# Description:`  
  **A brief explanation of what the file contains.**  
  Describes the purpose and context of the data.  
  Example: The relationship between internal IDs from the processing hub and the final project IDs.

- `# Created:`  
  **Date when the file was first created.**  
  Useful for version control and project timeline.  
  Example: `2025-03-26`

- `# Last Modified:`  
  **Date and time of the most recent modification.**  
  Indicates if the file is up to date.  
  Example: `2025-03-26 11:56`

- `# Author:`  
  **The person responsible for creating or maintaining the file.**  
  Helpful for identifying who to contact with questions.  
  Example: `Jose Antonio Corona-Gomez`

- `# Format:`  
  **The format of the data table.**  
  In this case: `TSV` (Tab-Separated Values), which is common for structured data.

- `# Columns:`  
  **A list of the column names in the table.**  
  Indicates what information each field contains and in what order.  
  Example: `NewID`, `OldID`, `Notes`

---

This type of standardized header is great for transparency, data sharing, and automating file parsing in pipelines or scripts.