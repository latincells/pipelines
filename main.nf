#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Mandatory Params
params.dataDir = "prueba/*" // Raw reads a analizar
params.outdir = "LatinCells_pipeline" // Directorio donde se guardaran archivos de interes
params.conda_envs_dir ="/media/storage2/software/anaconda3/envs/"
params.ref_data = "/media/storage2/software/refdata-gex-GRCh38-2020-A"
params.cellranger_path = "/media/storage2/software/cellranger-6.0.2/cellranger"
params.User = "$USER"
// Optional Params

// Channels definition
dataDir_ch = channel.fromFilePairs(params.dataDir, type: 'dir', size: -1) 
ref_data_ch = channel.fromPath(params.ref_data, type: 'dir') 
out_dir = file(params.outdir)
out_dir.mkdir()
params.Sample_metadata = ""
params.Mouse = false
process Cellranger { // 6.0.2
    publishDir "${params.outdir}/1-Counts", mode:'copy'  
    cpus 60
    maxForks 1
    tag "Cellranger on $sample_id"
    
    input: 
    tuple val(sample_id), path(reads)
    each path(ref_data)

    output:      
    tuple val(sample_id), path("Mapped_$sample_id")

    script:
    """
    cellranger count --id=Mapped_$sample_id --fastqs=$reads --sample=$sample_id --transcriptome=$ref_data --localcores $task.cpus
    """
}
process SoupX_scDblFinder { // Se necesita actualizar si es q la data no es de 10X, ya que el script de soupX esta simplificado y data de otro origen necesita la creacion de clusters 
    //conda "${params.conda_envs_dir}LatinCells"
    publishDir "${params.outdir}/2-SoupX_SCDblFinder", mode:'copy'
    cpus 30
    maxForks 4
    tag "SoupX on $sample_id"
    
    input: 
    tuple val(sample_id), path(reads)

    output:      
    path("$sample_id"), emit: clean_reads
    path("*png")

    script:
    """
    #!/usr/bin/env Rscript
    library(SoupX)
    library(DropletUtils)   
    sc = load10X("${reads}/outs/")
    png("soup_${sample_id}.png",width=480, height=480)
    sc = autoEstCont(sc)
    dev.off()
    out = adjustCounts(sc)
    DropletUtils:::write10xCounts("./${sample_id}_clean",out)
    library(Seurat)
    library(scDblFinder)
    library(dplyr) ###

    removeDoublets <- function(doubletFile, sobj) {
    to.remove <- read.table(doubletFile, sep = '\t', header = T) %>% 
        dplyr::filter(class == 'doublet') %>% 
        rownames()
    sobj <- subset(sobj, cells = setdiff(colnames(sobj), to.remove))
    return(sobj)
    }

    fnames <- "${sample_id}"
    fnames2 <- "${sample_id}_clean"

    data <- CreateSeuratObject(counts = Read10X(fnames2), min.cells = 3, min.features = 200)
    data <- RenameCells(data, add.cell.id = fnames)
    sce <- as.SingleCellExperiment(data)
    set.seed(123)
    results <- scDblFinder(sce, returnType = 'table') %>% 
        as.data.frame() %>% 
        filter(type == 'real')  
    write.table(results, file.path(paste0(fnames, '.txt')), sep = '\t', quote = F, col.names = T, row.names = T)

    datalist <- CreateSeuratObject(counts = Read10X(fnames2), project = fnames, min.cells = 3, min.features = 200)
    datalist <- RenameCells(datalist, add.cell.id = fnames)    
    datalist <- removeDoublets(file.path(paste0(fnames, '.txt')), datalist)
    write10xCounts(x = datalist@assays\$RNA@counts, path = fnames)
    """
}
process Scanpy {
    publishDir "${params.outdir}/Scanpy", mode:'copy'
    cpus 30
    maxForks 2
    conda "${params.conda_envs_dir}Scanpy"
    // conda create -n Scanpy -y
    // conda activate Scanpy
    // conda install -c conda-forge scanpy=1.8.1 -y
    // conda install -c bioconda scrublet=0.2.3 -y
    // conda install -c bioconda harmonypy=0.0.5 -y
    tag "Scanpy on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:      
    path("figures/*.png")

    script:
    """
    #!/usr/bin/env python3
    from itertools import groupby
    from pickle import TRUE
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import os
    import scrublet as scr
    import matplotlib.pyplot as plt
    from scipy import stats
    import scanpy.external as sce
    import harmonypy
    import seaborn as sns

    results_file = 'write/strainedCounts_${sample_id}.h5ad'  
    data_soup = sc.read_10x_mtx('$reads', var_names=('gene_symbols')) 
    data_soup.var['mt'] = data_soup.var_names.str.startswith('mt-') # mitochondrial genes
    data_soup.var['ribo'] = data_soup.var_names.str.startswith(("Rps","Rpl")) # ribosomal genes
    data_soup.var['hb'] = data_soup.var_names.str.contains(("^Hb[^(P)]")) # hemoglobin genes.
    sc.pp.calculate_qc_metrics(data_soup, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)

    lower_nmads_ribo=data_soup.obs.pct_counts_ribo.median() - (3 * stats.median_abs_deviation(data_soup.obs.pct_counts_ribo))
    upper_nmads_ribo=data_soup.obs.pct_counts_ribo.median() + (3 * stats.median_abs_deviation(data_soup.obs.pct_counts_ribo))
    lower_nmads_mt=data_soup.obs.pct_counts_mt.median() - (3 * stats.median_abs_deviation(data_soup.obs.pct_counts_mt))
    upper_nmads_mt=data_soup.obs.pct_counts_mt.median() + (3 * stats.median_abs_deviation(data_soup.obs.pct_counts_mt))
    lower_nmads_n_genes_by_counts=data_soup.obs.n_genes_by_counts.median() - (3 * stats.median_abs_deviation(data_soup.obs.n_genes_by_counts)) #“n_genes_by_counts”. The number of genes with at least 1 count in a cell. Calculated for all cells.
    upper_nmads_n_genes_by_counts=data_soup.obs.n_genes_by_counts.median() + (3 * stats.median_abs_deviation(data_soup.obs.n_genes_by_counts))
    lower_nmads_total_counts=data_soup.obs.total_counts.median() -(3 * stats.median_abs_deviation(data_soup.obs.total_counts))
    upper_nmads_total_counts=data_soup.obs.total_counts.median() +(3 * stats.median_abs_deviation(data_soup.obs.total_counts))
    filtered_data_soup = data_soup.obs.loc[(data_soup.obs.n_genes_by_counts > lower_nmads_n_genes_by_counts) & (data_soup.obs.n_genes_by_counts <  upper_nmads_n_genes_by_counts) & (data_soup.obs.pct_counts_mt < upper_nmads_mt) & (data_soup.obs.pct_counts_ribo< upper_nmads_ribo) & (data_soup.obs.total_counts > lower_nmads_total_counts)]    
    data_soup2 = data_soup[filtered_data_soup.index, :]
    sc.pl.violin(data_soup2, ['n_genes_by_counts'], jitter=0.4, multi_panel=True,save='_${sample_id}_soup.png') # Imagenes/WTs

    ## correr scrublet aqui, antes de concatenar, añadir al loop anterior? # se espera un 0.05
    #Predict doublets 
    # pip3 install scrublet --user 
    scrub = scr.Scrublet(data_soup2.X,expected_doublet_rate=0.06) # se trabaja con la matriz raw, no normalizada, trabajar en cada muestra por separado
    data_soup2.obs['doublet_scores'], data_soup2.obs['predicted_doublets'] = scrub.scrub_doublets()
    scrub.call_doublets(threshold=0.25)
    histo_doblet=scrub.plot_histogram()
    histo_doblet[0].savefig("figures/doublet_${sample_id}_soupscrub.png")  
    sum(data_soup2.obs['predicted_doublets'])
    # add in column with singlet/doublet instead of True/False
    data_soup2.obs['doublet_info'] = data_soup2.obs["predicted_doublets"].astype(str)
    sc.pl.violin(data_soup2, 'n_genes_by_counts',jitter=0.4, groupby = 'doublet_info', rotation=45,save='_${sample_id}_doublet.png')  
    #Now, lets run PCA and UMAP and plot doublet scores onto umap to check the doublet predictions
        
    data_soup3= data_soup2.copy()    
    sc.pp.log1p(data_soup3) #Logarithmize the data 
    sc.pp.highly_variable_genes(data_soup3, min_mean=0.0125, max_mean=3, min_disp=0.5)
    data_soup3 = data_soup3[:, data_soup3.var.highly_variable]
    sc.pp.regress_out(data_soup3, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(data_soup3, max_value=10)
    sc.tl.pca(data_soup3, svd_solver='arpack')
    sc.pp.neighbors(data_soup3, n_neighbors=10, n_pcs=40)
    sc.tl.umap(data_soup3)
    sc.pl.umap(data_soup3, color=['doublet_scores','doublet_info'],save='_${sample_id}_doublet2.png')
    """
}
process Monocle {
    cpus 30
    maxForks 2
    conda "${params.conda_envs_dir}Monocle"
    // conda create -n Monocle -y
    // conda activate Monocle
    // conda install -c bioconda r-monocle3=1.0.0 -y
    tag "Monocle on $sample_id"
    
    input:
    path(sample_id)

    output:    
    //path "*.html"
    //path "*zip", emit: multiqc_input    
    path "strainedCounts_*", emit: consensus_fasta  

    script:
    """
    Rscript --vanilla /media/storage/Adolfo/singleCell/SoupX.nf.R $sample_id
    """
}
process Seurat_Object_creation {
    //conda "${params.conda_envs_dir}LatinCells"
    publishDir "${params.outdir}/3-Seurat-Object", mode:'copy'

    input:
    path(reads)

    output:       
    path("datalist.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat) ###
    library(dplyr) ###
    fnames <- basename(list.dirs('./', recursive = F))
    files <- file.path(fnames)

    datalist <- list()
    if ("$params.Sample_metadata" != ""){
        meta <- read.delim("$params.Sample_metadata",sep=",")
        for (i in 1:length(files)) {
            datalist[[i]] <- CreateSeuratObject(counts = Read10X(files[i]), project = fnames[i], min.cells = 3, min.features = 200)  
            for (k in colnames(meta)){
                datalist[[i]] <- AddMetaData(datalist[[i]], metadata = meta[meta\$Sample == fnames[i],][[k]], col.name = k)
            }
            print(paste(i, 'out', length(files)))}
    } else {
        for (i in 1:length(files)) {
            datalist[[i]] <- CreateSeuratObject(counts = Read10X(files[i]), project = fnames[i], min.cells = 3, min.features = 200) 
            print(paste(i, 'out', length(files)))}
    }

    datalist <- Reduce(merge, datalist)
    saveRDS(datalist, 'datalist.rds')
    """
}
process Seurat_QC_integration {
    //conda "${params.conda_envs_dir}SingleCell"
    publishDir "${params.outdir}/4-Seurat_QC-Integration", mode:'copy'
    
    input:
    path(Datalist)

    output:    
    path("datalist.postQC-int*.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat) ###
    library(preprocessCore)
    library(dplyr) ###
    library(ggplot2) 
    #colors <- clusterExperiment::bigPalette
    colors <- c("#E31A1C","#1F78B4","#33A02C","#FF7F00","#6A3D9A","#B15928","#A6CEE3","#bd18ea","cyan","#B2DF8A","#FB9A99","deeppink4","#00B3FFFF","#CAB2D6","#FFFF99","#05188a","#CCFF00FF","cornflowerblue","#f4cc03","black","blueviolet","#4d0776","maroon3","blue","#E5D8BD","cadetblue4","#e5a25a","lightblue1","#F781BF","#FC8D62","#8DA0CB","#E78AC3","green3","#E7298A","burlywood3","#A6D854","firebrick","#FFFFCC","mediumpurple","#1B9E77","#FFD92F","deepskyblue4","yellow3","#00FFB2FF","#FDBF6F","#FDCDAC","gold3","#F4CAE4","#E6F5C9","#FF00E6FF","#7570B3","goldenrod","#85848f","lightpink3","olivedrab","cadetblue3")
    hk <- c("C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29")
    ribo <- c('RPLP0', 'RPLP1', 'RPLP2', 'RPL3', 'RPL3L', 'RPL4', 'RPL5', 'RPL6', 'RPL7', 'RPL7A', 'RPL7L1', 'RPL8', 'RPL9', 'RPL10', 'RPL10A', 'RPL11', 'RPL12', 'RPL13', 'RPL13A', 'RPL14', 'RPL15', 'RPL17', 'RPL18', 'RPL18A', 'RPL19', 'RPL21', 'RPL22', 'RPL22L1', 'RPL23', 'RPL23A', 'RPL24', 'RPL26', 'RPL26L1', 'RPL27', 'RPL27A', 'RPL28', 'RPL29', 'RPL30', 'RPL31', 'RPL32', 'RPL34', 'RPL35', 'RPL35A', 'RPL36', 'RPL36A', 'RPL36AL', 'RPL37', 'RPL37A', 'RPL38', 'RPL39', 'RPL39L', 'RPL41', 'RPSA', 'RPS2', 'RPS3', 'RPS3A', 'RPS4X', 'RPS4Y1', 'RPS5', 'RPS6', 'RPS7', 'RPS8', 'RPS9', 'RPS10', 'RPS11', 'RPS12', 'RPS13', 'RPS14', 'RPS15', 'RPS15A', 'RPS16', 'RPS17', 'RPS18', 'RPS19', 'RPS20', 'RPS21', 'RPS23', 'RPS24', 'RPS25', 'RPS26', 'RPS27', 'RPS27A', 'RPS28', 'RPS29')
    optimizePCA <- function(sobj, csum){
    dp <- Stdev(sobj)^2
    for (z in 1:length(dp)) {
        soma <- sum(dp[1:z])/sum(dp)
        if (soma >= csum) {
        best_pc <- z
        break()
        }
    }
    return(best_pc)
    }
    seurat_theme <- function(){
    theme_bw() +
        theme(panel.background = element_rect(colour = "black", size=0.1), 
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    datalist <- readRDS("${Datalist}")
    mt <- PercentageFeatureSet(datalist, pattern = "^MT-")
    
    tmp1 <- datalist %>% 
        AddModuleScore(features = list(hk), name = "housekeeping") %>% 
        AddModuleScore(features = list(ribo), name = "ribo.percent") %>% 
        AddMetaData(metadata = mt, col.name = 'percent.mt') 
    
    VlnPlot(tmp1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3) #,group.by = "Batch"
    ggsave(paste0("PreQC.png"),heigh=9,width=16) # .byBatch
    
    tmp1 <- tmp1 %>% 
        subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) %>% 
        NormalizeData(verbose = FALSE) %>% 
        FindVariableFeatures(verbose = FALSE) %>% 
        ScaleData(vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "housekeeping1", "ribo.percent1")) %>% 
        CellCycleScoring(s.features = cc.genes\$s.genes, g2m.features = cc.genes\$g2m.genes, set.ident = TRUE, search = TRUE) %>% 
        NormalizeData(verbose = FALSE) %>% 
        FindVariableFeatures(verbose = FALSE) %>% 
        ScaleData(vars.to.regress = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'housekeeping1', 'ribo.percent1', 'S.Score', 'G2M.Score')) %>% 
        RunPCA(verbose = FALSE)
    
    genes_remove <- grep("^HBB|^HBA|^MTRNR|^MT-|^MALAT1|^NEAT1|^RPS|^RPL", x = rownames(GetAssayData(tmp1, 'scale.data')), value = TRUE)
    genes_keep <- setdiff(rownames(GetAssayData(tmp1, 'scale.data')), genes_remove)
    tmp1@assays\$RNA@scale.data <- tmp1@assays\$RNA@scale.data[genes_keep, ]
    
    mat <- as.matrix(GetAssayData(tmp1, 'scale.data'))
    test <- normalize.quantiles(mat)
    rownames(test) <- rownames(mat)
    colnames(test) <- colnames(mat)
    rm(mat)
    tmp1@assays\$RNA@scale.data <- test
    rm(test)
    
    gc(); gc(); gc()
    
    tmp1 <- tmp1 %>% 
        RunPCA() %>% 
        FindNeighbors(reduction = 'pca', dims = 1:optimizePCA(tmp1, 0.8), k.param = 30) %>% 
        FindClusters(resolution = seq(0.6, 2, 0.2)) %>% 
        RunUMAP(reduction = 'pca', dims = 1:optimizePCA(tmp1, 0.8))
    
    if (sum(grepl('Seurat_clusters', colnames(tmp1@meta.data))) > 0) {
        tmp1@meta.data <- tmp1@meta.data[,-grep('Seurat_clusters', colnames(tmp1@meta.data))]
    }
    
    colnames(tmp1@meta.data) <- gsub('RNA_snn', 'Seurat_clusters', colnames(tmp1@meta.data))
    datalist <- tmp1 
    rm(tmp1)

    #datalist <- ppSeurat(datalist, res = seq(0.6, 2, 0.1))
    saveRDS(datalist, 'datalist.postQC-int.rds')
    """
}
process Seurat_Cell_Annotation {
    //conda "${params.conda_envs_dir}SingleCell"
    publishDir "${params.outdir}/5-Seurat-Cell_Annotation", mode:'copy'
    
    input:
    path(Datalist)

    output:    
    path("datalist.Cell_Annotation*.rds"), emit: rds 
    path("*.png")

    script:
    """
    #!/usr/bin/env Rscript
    options(timeout=10000000000000000000000000)
    library(rols)
    cl <- rols::Ontology("cl")
    library(celldex)
    library(Rtsne)
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(cowplot)
    library(clustree)
    library(BiocParallel)
    library(ggrepel)
    library(SingleR)
    ## Run just one time
    ref1 <- celldex::BlueprintEncodeData()
    ref2 <- celldex::DatabaseImmuneCellExpressionData()
    ref3 <- celldex::HumanPrimaryCellAtlasData()
    ref4 <- celldex::MonacoImmuneData()
    shared <- Reduce(intersect, list(rownames(ref1), rownames(ref2), rownames(ref3), rownames(ref4)))
    ref1 <- ref1[shared,]; ref2 <- ref2[shared,]; ref3 <- ref3[shared,]; ref4 <- ref4[shared,]
    save(ref1, ref2, ref3, ref4, cl, file = 'SingleR_refs.RData')
    #load("Data/SingleR_refs.RData")
    set.seed(123)
    colors <- c("#E31A1C","#1F78B4","#33A02C","#FF7F00","#6A3D9A","#B15928","#A6CEE3","#bd18ea","cyan","#B2DF8A","#FB9A99","deeppink4","#00B3FFFF","#CAB2D6","#FFFF99","#05188a","#CCFF00FF","cornflowerblue","#f4cc03","black","blueviolet","#4d0776","maroon3","blue","#E5D8BD","cadetblue4","#e5a25a","lightblue1","#F781BF","#FC8D62","#8DA0CB","#E78AC3","green3","#E7298A","burlywood3","#A6D854","firebrick","#FFFFCC","mediumpurple","#1B9E77","#FFD92F","deepskyblue4","yellow3","#00FFB2FF","#FDBF6F","#FDCDAC","gold3","#F4CAE4","#E6F5C9","#FF00E6FF","#7570B3","goldenrod","#85848f","lightpink3","olivedrab","cadetblue3")
    seurat_theme <- function(){
        theme_bw() +
            theme(panel.background = element_rect(colour = "black", size=0.1), 
                plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
                axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    ppSingleR <- function(sobj, cluster){
    sce <- as.SingleCellExperiment(sobj)
    
    pred <- SingleR(test = sce, 
                    ref = list(ref1, ref2, ref3, ref4), 
                    labels = list(ref1\$label.ont, ref2\$label.ont, ref3\$label.ont, ref4\$label.ont), 
                    clusters = sce[[cluster]], BPPARAM = MulticoreParam(1))
    
    sobj <- SetIdent(sobj, value = cluster)
    new.cluster.ids <- pred\$labels
    names(new.cluster.ids) <- levels(Idents(sobj))
    sobj <- RenameIdents(sobj, new.cluster.ids)
    sobj[[paste0('SingleR_ont', gsub('Seurat_clusters', '', cluster))]] <- Idents(sobj)
    
    onts <- pred\$labels
    new.cluster.ids <- c()
    for (z in 1:length(onts)) {
        new.cluster.ids <- c(new.cluster.ids, termLabel(term(cl, onts[z])))
    }
    sobj <- SetIdent(sobj, value = cluster)
    names(new.cluster.ids) <- levels(Idents(sobj))
    sobj <- RenameIdents(sobj, new.cluster.ids)
    sobj[[paste0('SingleR', gsub('Seurat_clusters', '', cluster))]] <- Idents(sobj)
    
    return(sobj)
    }
    datalist <- readRDS("${Datalist}")
    ###################################################################################################################################################
    #                                                        SingleR                                                                                 #
    ###################################################################################################################################################
    datalist
    clusters <- grep('Seurat_clusters', names(datalist@meta.data), value = T)
    for (cluster in clusters) {
    datalist <- ppSingleR(datalist, cluster)
    print(cluster)
    }
    clustree(datalist,prefix="SingleR_res.",node_text_angle = 15)
    ggsave("clustree_singleR.png", width = 16,height = 9)
    ###################################################################################################################################################
    DimPlot(datalist,reduction = "umap", group.by = "SingleR_res.1.4", cols = colors,raster=FALSE) + seurat_theme()
    ggsave(paste0("SingleR_res.1.4.png"),heigh=9,width=16)
    DimPlot(datalist,reduction = "umap", group.by = "SingleR_res.0.8", cols = colors,raster=FALSE) + seurat_theme()
    ggsave(paste0("SingleR_res.0.8.png"),heigh=9,width=16)
    ###################################################################################################################################################
    clustree(datalist,prefix="Seurat_clusters_res.",node_text_angle = 15)
    ggsave("clustree_Seurat_clusters.png", width = 16,height = 9)
    DimPlot(datalist,reduction = "umap", group.by = "Seurat_clusters_res.1.4", cols = colors,raster=FALSE) + seurat_theme()
    ggsave("Seurat_clusters_res.1.4.png",heigh=9,width=16)
    ###################################################################################################################################################
    saveRDS(datalist, file = "datalist.Cell_Annotation.rds")
    """
}
process Seurat_create_object_mouse {
    //conda "${params.conda_envs_dir}SingleCell"
    publishDir "${params.outdir}/Seurat", mode:'copy'
    
    input:
    path(reads)
    path(opt)

    output:     
    path("datalist.rds"), emit: rds 
    path("datalist.h5ad"), emit: h5ad 

    script:
    def Sample_metadata_file = opt.name != "" ? "$opt" : ''
    """
    #!/usr/bin/env Rscript
    library(Seurat) ###
    library(preprocessCore)
    library(SeuratDisk)
    library(scDblFinder)###
    library(dplyr) ###
    library(gprofiler2)
    colors <- clusterExperiment::bigPalette

    hk <- c("C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4", "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29")
    hk <-  gorth(hk, source_organism = "hsapiens", target_organism = "mmusculus")\$ortholog_name
    ribo <- c('RPLP0', 'RPLP1', 'RPLP2', 'RPL3', 'RPL3L', 'RPL4', 'RPL5', 'RPL6', 'RPL7', 'RPL7A', 'RPL7L1', 'RPL8', 'RPL9', 'RPL10', 'RPL10A', 'RPL11', 'RPL12', 'RPL13', 'RPL13A', 'RPL14', 'RPL15', 'RPL17', 'RPL18', 'RPL18A', 'RPL19', 'RPL21', 'RPL22', 'RPL22L1', 'RPL23', 'RPL23A', 'RPL24', 'RPL26', 'RPL26L1', 'RPL27', 'RPL27A', 'RPL28', 'RPL29', 'RPL30', 'RPL31', 'RPL32', 'RPL34', 'RPL35', 'RPL35A', 'RPL36', 'RPL36A', 'RPL36AL', 'RPL37', 'RPL37A', 'RPL38', 'RPL39', 'RPL39L', 'RPL41', 'RPSA', 'RPS2', 'RPS3', 'RPS3A', 'RPS4X', 'RPS4Y1', 'RPS5', 'RPS6', 'RPS7', 'RPS8', 'RPS9', 'RPS10', 'RPS11', 'RPS12', 'RPS13', 'RPS14', 'RPS15', 'RPS15A', 'RPS16', 'RPS17', 'RPS18', 'RPS19', 'RPS20', 'RPS21', 'RPS23', 'RPS24', 'RPS25', 'RPS26', 'RPS27', 'RPS27A', 'RPS28', 'RPS29')
    ribo <-  gorth(ribo, source_organism = "hsapiens", target_organism = "mmusculus")\$ortholog_name
    cc.genes\$s.genes <-  gorth(cc.genes\$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")\$ortholog_name
    cc.genes\$g2m.genes <-  gorth(cc.genes\$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")\$ortholog_name
    optimizePCA <- function(sobj, csum){
    dp <- Stdev(sobj)^2
    for (z in 1:length(dp)) {
        soma <- sum(dp[1:z])/sum(dp)
        if (soma >= csum) {
        best_pc <- z
        break()
        }
    }
    return(best_pc)
    }

    removeDoublets <- function(doubletFile, sobj) {
    to.remove <- read.table(doubletFile, sep = '\t', header = T) %>% 
        dplyr::filter(class == 'doublet') %>% 
        rownames()
    sobj <- subset(sobj, cells = setdiff(colnames(sobj), to.remove))
    return(sobj)
    }

    seurat_theme <- function(){
    theme_bw() +
        theme(panel.background = element_rect(colour = "black", size=0.1), 
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    ppSeurat <- function(sobj, res){
    mt <- PercentageFeatureSet(sobj, pattern = "^mt-")
    
    tmp1 <- sobj %>% 
        AddModuleScore(features = list(hk), name = "housekeeping") %>% 
        AddModuleScore(features = list(ribo), name = "ribo.percent") %>% 
        AddMetaData(metadata = mt, col.name = 'percent.mt') %>% 
        subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10) %>% 
        NormalizeData(verbose = FALSE) %>% 
        FindVariableFeatures(verbose = FALSE) %>% 
        ScaleData(vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt", "housekeeping1", "ribo.percent1")) %>% 
        CellCycleScoring(s.features = cc.genes\$s.genes, g2m.features = cc.genes\$g2m.genes, set.ident = TRUE, search = TRUE) %>% 
        NormalizeData(verbose = FALSE) %>% 
        FindVariableFeatures(verbose = FALSE) %>% 
        ScaleData(vars.to.regress = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'housekeeping1', 'ribo.percent1', 'S.Score', 'G2M.Score')) %>% 
        RunPCA(verbose = FALSE)
    
    genes_remove <- grep("^Hbb|^Hba|^Mtrnr|^mt-|^Malat1|^Neat1|^Rps|^Rpl", x = rownames(GetAssayData(tmp1, 'scale.data')), value = TRUE)
    genes_keep <- setdiff(rownames(GetAssayData(tmp1, 'scale.data')), genes_remove)
    tmp1@assays\$RNA@scale.data <- tmp1@assays\$RNA@scale.data[genes_keep, ]
    
    mat <- as.matrix(GetAssayData(tmp1, 'scale.data'))
    test <- normalize.quantiles(mat)
    rownames(test) <- rownames(mat)
    colnames(test) <- colnames(mat)
    rm(mat)
    tmp1@assays\$RNA@scale.data <- test
    rm(test)
    
    gc(); gc(); gc()
    
    tmp1 <- tmp1 %>% 
        RunPCA() %>% 
        FindNeighbors(reduction = 'pca', dims = 1:optimizePCA(tmp1, 0.8), k.param = 30) %>% 
        FindClusters(resolution = res) %>% 
        RunUMAP(reduction = 'pca', dims = 1:optimizePCA(tmp1, 0.8))
    
    if (sum(grepl('Seurat_clusters', colnames(tmp1@meta.data))) > 0) {
        tmp1@meta.data <- tmp1@meta.data[,-grep('Seurat_clusters', colnames(tmp1@meta.data))]
    }
    
    colnames(tmp1@meta.data) <- gsub('RNA_snn', 'Seurat_clusters', colnames(tmp1@meta.data))
    
    return(tmp1)
    }

    fnames <- basename(list.dirs('./', recursive = F))
    files <- file.path(fnames)

    for (i in 1:length(files)) {
    data <- CreateSeuratObject(counts = Read10X(files[i]), min.cells = 3, min.features = 200)
    data <- RenameCells(data, add.cell.id = fnames[i])
    
    sce <- as.SingleCellExperiment(data)
    set.seed(123)
    results <- scDblFinder(sce, returnType = 'table') %>% 
        as.data.frame() %>% 
        filter(type == 'real')  
    write.table(results, file.path(paste0(fnames[i], '.txt')), sep = '\t', quote = F, col.names = T, row.names = T)
    }

    fnames <- basename(list.dirs('./', recursive = F))
    files <- file.path(fnames)

    datalist <- list()
    if ("$Sample_metadata_file" != ""){
        meta <- read.delim("$Sample_metadata_file",sep=",")
        for (i in 1:length(files)) {
            datalist[[i]] <- CreateSeuratObject(counts = Read10X(files[i]), project = fnames[i], min.cells = 3, min.features = 200)
            datalist[[i]] <- RenameCells(datalist[[i]], add.cell.id = fnames[i])    
            datalist[[i]] <- removeDoublets(file.path(paste0(fnames[i], '.txt')), datalist[[i]])
            for (k in colnames(meta)){
                datalist[[i]] <- AddMetaData(datalist[[i]], metadata = meta[meta\$Sample == fnames[i],][[k]], col.name = k)
            }
            print(paste(i, 'out', length(files)))}
    } else {
        for (i in 1:length(files)) {
            datalist[[i]] <- CreateSeuratObject(counts = Read10X(files[i]), project = fnames[i], min.cells = 3, min.features = 200)
            datalist[[i]] <- RenameCells(datalist[[i]], add.cell.id = fnames[i])    
            datalist[[i]] <- removeDoublets(file.path(paste0(fnames[i], '.txt')), datalist[[i]])
            print(paste(i, 'out', length(files)))}
    }

    datalist <- Reduce(merge, datalist)
    datalist <- ppSeurat(datalist, res = seq(0.3, 1, 0.1))
    saveRDS(datalist, 'datalist.rds')
    SaveH5Seurat(datalist, filename = "datalist.h5Seurat")
    Convert("datalist.h5Seurat", dest = "h5ad")
    file.remove("datalist.h5Seurat")
    """
}
workflow {    
    Cellranger(dataDir_ch,ref_data_ch)
    SoupX_scDblFinder(Cellranger.out)
    if (params.Mouse != false){ 
    Seurat_create_object_mouse(SoupX_scDblFinder.out.clean_reads.collect())
    } else {
    Seurat_Object_creation(SoupX_scDblFinder.out.clean_reads.collect())   
    Seurat_QC_integration(Seurat_Object_creation.out.rds)    
    Seurat_Cell_Annotation(Seurat_QC_integration.out.rds)
    }
    //Scanpy(SoupX.out.clean_reads)
}
workflow prueba {    
    print(params.User)
}
//nextflow /media/storage2/Adolfo2/LatinCells/LatinCells/main.nf '--dataDir=/media/storage2/Adolfo2/LatinCells/workflow_test/RawData_test/*' -resume