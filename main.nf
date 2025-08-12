#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Mandatory Params
params.dataDir = "Prueba/*" // Raw reads a analizar
params.ref_data = "/media/storage2/software/refdata-gex-GRCh38-2020-A"
params.Project_Name= "General"
params.MainOutdir = "LatinCells_Results"  // Directorio donde se guardaran archivos de interes
params.outdir = "${params.MainOutdir}/${params.Project_Name}/"  // Directorio donde se guardaran archivos de interes
// Optional Params
params.Batch = "orig.ident"
params.NPlex= ""
params.Sample_metadata = ""
params.Mouse = false
params.Freemuxlet = false
params.VCF_Files = "" 
params.ModelsCelltypist = ""
params.help = false
process Cellranger {
    storeDir "${params.outdir}/1-Counts"
    maxForks 2
    tag "on $sample_id"
    errorStrategy 'ignore'
    
    input: 
    tuple val(sample_id), path(reads)
    each path(ref_data)

    output:      
    tuple val(sample_id), path("Mapped_$sample_id")

    script:
    """
    cellranger count --id=Mapped_${sample_id} --fastqs=${reads} --sample=${sample_id} --transcriptome=${ref_data} --create-bam=true --localcores=${task.cpus}
    mv Mapped_${sample_id}/outs/filtered_feature_bc_matrix.h5 Mapped_${sample_id}/outs/${sample_id}.filtered_feature_bc_matrix.h5
    mv Mapped_${sample_id}/outs/web_summary.html Mapped_${sample_id}/outs/${sample_id}.web_summary.html
    mv Mapped_${sample_id}/outs/raw_feature_bc_matrix.h5 Mapped_${sample_id}/outs/${sample_id}.raw_feature_bc_matrix.h5
    """
}
process Freemuxlet {
    label 'Demultiplex'
    publishDir "${params.outdir}/Metadata", mode:'copy'  
    maxForks 2
    tag "on $sample_id"
    errorStrategy 'ignore'
    
    input: 
    tuple val(sample_id), path(reads)

    output:      
    tuple val(sample_id), path("Freemuxlet_outs/${sample_id}.Freemuxlet.samples")

    script:
    """
    bedtools merge -i /opt/GRCh38_1000G_MAF0.05_ExonFiltered_ChrEncoding.sorted.vcf | samtools view -@ 8 --write-index -L - -D CB:<(zcat ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz) -o filtered_bam.bam ${reads}/outs/possorted_genome_bam.bam
    mkdir Pileup_files
    popscle dsc-pileup --sam filtered_bam.bam --vcf /opt/GRCh38_1000G_MAF0.05_ExonFiltered_ChrEncoding.sorted.vcf --group-list ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --out Pileup_files/${sample_id}
    mkdir Freemuxlet_outs
    popscle freemuxlet --plp Pileup_files/${sample_id} --out Freemuxlet_outs/${sample_id} --group-list ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --nsample ${params.NPlex}
    gunzip Freemuxlet_outs/${sample_id}.clust1.samples.gz
    mv Freemuxlet_outs/${sample_id}.clust1.samples Freemuxlet_outs/${sample_id}.Freemuxlet.samples
    """
}
process CombinedFremuxletDemuxlet {
    label 'Demultiplex'
    publishDir "${params.outdir}/Metadata", mode:'copy'  
    maxForks 2
    tag "on $sample_id"
    errorStrategy 'ignore'
    
    input: 
    tuple val(sample_id), path(reads), path(vcf)

    output:      
    tuple val(sample_id), path("${sample_id}.CombinedDemultiplex.samples")

    script:
    """
    NPlex=\$(gunzip -c ${vcf} | grep -m 1 "^#CHROM" | awk '{print NF-9}')
    bedtools merge -i ${vcf} | samtools view -@ 8 --write-index -L - -D CB:<(zcat ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz) -o filtered_Demuxlet_bam.bam ${reads}/outs/possorted_genome_bam.bam &
    bedtools merge -i /opt/GRCh38_1000G_MAF0.05_ExonFiltered_ChrEncoding.sorted.vcf | samtools view -@ 8 --write-index -L - -D CB:<(zcat ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz) -o filtered_Freemuxlet_bam.bam ${reads}/outs/possorted_genome_bam.bam
    mkdir Pileup_files
    popscle dsc-pileup --sam filtered_Freemuxlet_bam.bam --vcf /opt/GRCh38_1000G_MAF0.05_ExonFiltered_ChrEncoding.sorted.vcf --group-list ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --out Pileup_files/${sample_id}
    rm filtered_Freemuxlet_bam.bam &
    popscle freemuxlet --plp Pileup_files/${sample_id} --out ${sample_id} --group-list ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --nsample \$NPlex
    gunzip ${sample_id}.clust1.samples.gz
    popscle demuxlet --sam filtered_Demuxlet_bam.bam --group-list ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --vcf ${vcf} --field GT --out ${sample_id} --alpha 0 --alpha 0.5
    (head -n 1 ${sample_id}.clust1.samples && tail -n +2 ${sample_id}.clust1.samples && tail -n +2 ${sample_id}.best) > ${sample_id}.CombinedDemultiplex.samples
    rm filtered_Demuxlet_bam.bam 
    """
}
process Demuxlet {
    label 'Demultiplex'
    publishDir "${params.outdir}/Metadata", mode:'copy'  
    maxForks 2
    tag "on $sample_id"
    
    input: 
    tuple val(sample_id), path(reads), path(vcf)

    output:      
    tuple val(sample_id), path("${sample_id}.Demuxlet.best")

    script:
    """
    mkdir Demuxlet_outs
    popscle demuxlet --sam ${reads}/outs/possorted_genome_bam.bam --group-list ${reads}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --vcf ${vcf} --field GT --out ${sample_id} --alpha 0 --alpha 0.5
    mv ${sample_id}.best ${sample_id}.Demuxlet.best
    """
}
process MappingStats {
    publishDir "${params.outdir}/Metadata", mode:'copy'  
    
    input: 
    path(CellrangerOuts)

    output:      
    path("CellrangerRunsSummary.tsv")

    script:
    """
    #!/usr/bin/env Rscript
    DataFrame <- data.frame()
    for (dir in list.dirs(path = "./", full.names = TRUE, recursive = F)) {
        print(dir)
        metrics_summary <- read.csv(paste0(dir, "/outs/metrics_summary.csv"))
        print(metrics_summary)
        RunID <- data.frame("RunID" = basename(dir))
        DataFrame <- rbind(DataFrame,cbind(RunID, metrics_summary))
    }

    write.table(DataFrame, "CellrangerRunsSummary.tsv", row.names = FALSE, quote = FALSE, sep = "\\t")
    """
}
process SoupX_scDblFinder {
    publishDir "${params.outdir}/2-SoupX_SCDblFinder", mode:'copy', pattern: '*png'
    maxForks 4
    tag "on $sample_id"
    errorStrategy 'ignore'
    
    input: 
    tuple val(sample_id), path(reads)
    file(demultiplexData)

    output:      
    path("$sample_id"), emit: clean_reads
    path("${sample_id}_StatsTable.csv"), emit: Stats
    path("*png")

    script:
    """
    #!/usr/bin/env Rscript
    library(SoupX)
    library(DropletUtils)
    set.seed(123)
    DemultiplexingCombinations <- sort(c("AMB-AMB", "AMB-DBL", "AMB-SNG", 
                                         "DBL-AMB", "DBL-DBL", "DBL-SNG", 
                                         "SNG-AMB", "SNG-DBL", "SNG-SNG"))
    colors <- c("#E31A1C","#33A02C","#FF7F00","#6A3D9A","#B15928","olivedrab","#bd18ea","cyan","green3")
    names(colors) <- DemultiplexingCombinations
    seurat_theme <- function(){
    theme_bw() +
        theme(panel.background = element_rect(colour = "black", size=0.1), 
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(colour = NA, fill = NA))
    }   
    sc = load10X("${reads}/outs/")
    StatsTable <- data.frame("RunID"="${sample_id}", "CellRanger" = ncol(sc\$toc))
    png("SoupX_${sample_id}.png",width=1080, height=1080)
    sc = autoEstCont(sc)
    dev.off()
    rhoEst <- sc\$fit\$rhoEst
    StatsTable\$rhoEst <- rhoEst
    out = adjustCounts(sc, roundToInt = TRUE)
    DropletUtils:::write10xCounts("./${sample_id}_clean",out)
    library(Seurat)
    library(scDblFinder)
    library(dplyr) ###
    library(ggalluvial)

    removeDoublets <- function(doubletFile, sobj) {
    to.remove <- read.table(doubletFile, sep = '\\t', header = T) %>% 
        dplyr::filter(class == 'doublet') %>% 
        rownames()
    sobj <- subset(sobj, cells = setdiff(colnames(sobj), to.remove))
    return(sobj)
    }

    fnames <- "${sample_id}"
    fnames2 <- "${sample_id}_clean"

    data <- CreateSeuratObject(counts = Read10X(fnames2), min.cells = 3, min.features = 200)
    StatsTable\$SeuratInitial <- ncol(data)
    if (file.exists(paste0("${sample_id}",".Freemuxlet.samples")) | file.exists(paste0("${sample_id}",".Demuxlet.best"))){
        if (file.exists(paste0("${sample_id}",".Freemuxlet.samples"))){
            MultiPlexData <- read.delim(paste0("${sample_id}",".Freemuxlet.samples"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]
        } else {
            MultiPlexData <- read.delim(paste0("${sample_id}",".Demuxlet.best"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]}
        rownames(MultiPlexData) <- MultiPlexData\$BARCODE
        MultiPlexData\$BARCODE <- NULL
        data <- AddMetaData(data, metadata = MultiPlexData)
        data <- subset(data, DROPLET.TYPE == "SNG")
        StatsTable\$nSNG_Plex <- ncol(data)
    }
    if (file.exists(paste0("${sample_id}",".CombinedDemultiplex.samples"))){
        Combined <- read.delim(paste0("${sample_id}",".CombinedDemultiplex.samples"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]
        FreemuxletData <- Combined[!is.na(as.numeric(Combined\$SNG.BEST.GUESS)),]
        names(FreemuxletData) <- c("BARCODE","Freemuxlet.DROPLET.TYPE","Freemuxlet.BEST.GUESS")
        DemuxletData <- Combined[is.na(as.numeric(Combined\$SNG.BEST.GUESS)),]
        names(DemuxletData) <- c("BARCODE","Demuxlet.DROPLET.TYPE","Demuxlet.BEST.GUESS")
        Joined_table <- merge(DemuxletData, FreemuxletData, by="BARCODE", all=TRUE)
        Joined_table\$DROPLET.TYPE <- paste(Joined_table\$Demuxlet.DROPLET.TYPE, Joined_table\$Freemuxlet.DROPLET.TYPE, sep = "-")
        PlotData <- as.data.frame(table(Joined_table\$Demuxlet.BEST.GUESS,Joined_table\$Freemuxlet.BEST.GUESS,Joined_table\$DROPLET.TYPE))
        names(PlotData) <- c("Demuxlet.BEST.GUESS","Freemuxlet.BEST.GUESS","DROPLET.TYPE","Freq")
        PlotData <- PlotData[PlotData\$Freq > 0,]
        Frequencies <- Joined_table %>% filter(Freemuxlet.DROPLET.TYPE == "SNG") %>% group_by(Freemuxlet.BEST.GUESS, Demuxlet.BEST.GUESS) %>% 
                        summarise(Frequency = n(), .groups = 'drop')
        Assigned_identity <- Frequencies %>% 
                                group_by(Freemuxlet.BEST.GUESS) %>% 
                                slice_max(Frequency, n = 1) %>% 
                                ungroup() %>% select(Freemuxlet.BEST.GUESS, Demuxlet.BEST.GUESS)
        for (i in Assigned_identity\$Freemuxlet.BEST.GUESS){
            FreemuxletData\$Freemuxlet.BEST.GUESS[FreemuxletData\$Freemuxlet.BEST.GUESS == i] <- Assigned_identity\$Demuxlet.BEST.GUESS[Assigned_identity\$Freemuxlet.BEST.GUESS == i]
        }
        names(FreemuxletData) <- c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")
        ggplot(as.data.frame(PlotData),
            aes(y = Freq, axis1 = Demuxlet.BEST.GUESS, axis2 = Freemuxlet.BEST.GUESS)) +
            geom_alluvium(aes(fill = DROPLET.TYPE), width = 1/12) +
            geom_stratum(width = 1/12, fill = "grey", color = "black") +
            geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limits = c("Demuxlet.BEST.GUESS", "Freemuxlet.BEST.GUESS"), expand = c(.05, .05)) +
            scale_fill_manual(values = colors) +
            ggtitle("Demuxlet vs Freemuxlet") + seurat_theme()
        ggsave("DemultiplexingComparison_${sample_id}.png", width = 21, height = 20, dpi = 300)
        MultiPlexData <- FreemuxletData
        rownames(MultiPlexData) <- MultiPlexData\$BARCODE
        MultiPlexData\$BARCODE <- NULL
        data <- AddMetaData(data, metadata = MultiPlexData)
        data <- subset(data, DROPLET.TYPE == "SNG")
        StatsTable\$nSNG_Plex <- ncol(data)
    }
    data <- RenameCells(data, add.cell.id = fnames)
    sce <- as.SingleCellExperiment(data)
    results <- scDblFinder(sce, returnType = 'table') %>% as.data.frame() %>% filter(type == 'real')  
    write.table(results, file.path(paste0(fnames, '.txt')), sep = '\\t', quote = F, col.names = T, row.names = T)

    datalist <- CreateSeuratObject(counts = Read10X(fnames2), project = fnames, min.cells = 3, min.features = 200)
    if (file.exists(paste0("${sample_id}",".Freemuxlet.samples")) | file.exists(paste0("${sample_id}",".Demuxlet.best"))){
        if (file.exists(paste0("${sample_id}",".Freemuxlet.samples"))){
            MultiPlexData <- read.delim(paste0("${sample_id}",".Freemuxlet.samples"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]
        } else {
            MultiPlexData <- read.delim(paste0("${sample_id}",".Demuxlet.best"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]}
        rownames(MultiPlexData) <- MultiPlexData\$BARCODE
        MultiPlexData\$BARCODE <- NULL
        datalist <- AddMetaData(datalist, metadata = MultiPlexData)
        datalist <- subset(datalist, DROPLET.TYPE == "SNG")
    }
    if (file.exists(paste0("${sample_id}",".CombinedDemultiplex.samples"))){
        Combined <- read.delim(paste0("${sample_id}",".CombinedDemultiplex.samples"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]
        FreemuxletData <- Combined[!is.na(as.numeric(Combined\$SNG.BEST.GUESS)),]
        FreemuxletDataSNG <- FreemuxletData[FreemuxletData\$DROPLET.TYPE == "SNG",]
        DemuxletData <- Combined[is.na(as.numeric(Combined\$SNG.BEST.GUESS)),][c("BARCODE","SNG.BEST.GUESS")]
        names(DemuxletData) <- c("BARCODE","Identity")
        Joined_table <- inner_join(FreemuxletDataSNG, DemuxletData, by = "BARCODE")
        Frequencies <- Joined_table %>% group_by(SNG.BEST.GUESS, Identity) %>% summarise(Frequency = n(), .groups = 'drop')
        Assigned_identity <- Frequencies %>% group_by(SNG.BEST.GUESS) %>% slice_max(Frequency, n = 1) %>% ungroup() %>% select(SNG.BEST.GUESS, Identity)
        for (i in Assigned_identity\$SNG.BEST.GUESS){
            FreemuxletData\$SNG.BEST.GUESS[FreemuxletData\$SNG.BEST.GUESS == i] <- Assigned_identity\$Identity[Assigned_identity\$SNG.BEST.GUESS == i]}
        MultiPlexData <- FreemuxletData
        rownames(MultiPlexData) <- MultiPlexData\$BARCODE
        MultiPlexData\$BARCODE <- NULL
        datalist <- AddMetaData(datalist, metadata = MultiPlexData)
        datalist <- subset(datalist, DROPLET.TYPE == "SNG")
    }
    datalist <- RenameCells(datalist, add.cell.id = fnames)    
    datalist <- removeDoublets(file.path(paste0(fnames, '.txt')), datalist)
    StatsTable\$nSNG_scDblFinder <- ncol(datalist)
    write.csv(StatsTable, "${sample_id}_StatsTable.csv", row.names = F)
    sce2 <- as.SingleCellExperiment(datalist)
    DropletUtils:::write10xCounts(x = sce2@assays@data\$counts, path = fnames)
    """
}
process Seurat_Object_creation {
    publishDir "${params.outdir}/3-Seurat-Object", mode:'copy'

    input:
    path(reads)
    file(demultiplexData)
    file(metadata)

    output:       
    path("datalist.rds"), emit: rds 

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(dplyr)
    set.seed(123)   
    fnames <- basename(list.dirs('./', recursive = F))
    files <- file.path(fnames)

    meta <- data.frame()
    datalist <- list()
    if ("${params.Sample_metadata}" != ""){
        if(grepl(".tsv","${params.Sample_metadata}")){
        meta <- read.delim(basename("${params.Sample_metadata}"),sep="\\t")}
        if(grepl(".tab","${params.Sample_metadata}")){
        meta <- read.delim(basename("${params.Sample_metadata}"),sep="\\t")}
        if(grepl(".csv","${params.Sample_metadata}")){
        meta <- read.delim(basename("${params.Sample_metadata}"),sep=",")}
        }

    Multiplex <- "No"

    for (i in 1:length(files)) {
        datalist[[i]] <- CreateSeuratObject(counts = Read10X(files[i]), project = fnames[i], min.cells = 3, min.features = 200)  
        if (file.exists(paste0(files[i],".Freemuxlet.samples")) | file.exists(paste0(files[i],".Demuxlet.best"))){
            if (file.exists(paste0(files[i],".Freemuxlet.samples"))){
                MultiPlexData <- read.delim(paste0(files[i],".Freemuxlet.samples"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]
            } else {
                MultiPlexData <- read.delim(paste0(files[i],".Demuxlet.best"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]}
            MultiPlexData\$BARCODE <- paste0(files[i],"_",MultiPlexData\$BARCODE)
            if ("${params.Sample_metadata}" != ""){
                MultiPlexData <- merge(MultiPlexData, meta, by.x = "SNG.BEST.GUESS", by.y = "Sample")
            }
            rownames(MultiPlexData) <- MultiPlexData\$BARCODE
            MultiPlexData\$BARCODE <- NULL
            datalist[[i]] <- AddMetaData(datalist[[i]], metadata = MultiPlexData)
            Multiplex <- "Si"
        }
        if (file.exists(paste0(files[i],".CombinedDemultiplex.samples"))){
                Combined <- read.delim(paste0(files[i],".CombinedDemultiplex.samples"),sep="\\t")[c("BARCODE","DROPLET.TYPE","SNG.BEST.GUESS")]
                FreemuxletData <- Combined[!is.na(as.numeric(Combined\$SNG.BEST.GUESS)),]
                FreemuxletDataSNG <- FreemuxletData[FreemuxletData\$DROPLET.TYPE == "SNG",]
                DemuxletData <- Combined[is.na(as.numeric(Combined\$SNG.BEST.GUESS)),][c("BARCODE","SNG.BEST.GUESS")]
                names(DemuxletData) <- c("BARCODE","Identity")
                Joined_table <- inner_join(FreemuxletDataSNG, DemuxletData, by = "BARCODE")
                Frequencies <- Joined_table %>% group_by(SNG.BEST.GUESS, Identity) %>% summarise(Frequency = n(), .groups = 'drop')
                Assigned_identity <- Frequencies %>% group_by(SNG.BEST.GUESS) %>% slice_max(Frequency, n = 1) %>% ungroup() %>% select(SNG.BEST.GUESS, Identity)
                for (k in Assigned_identity\$SNG.BEST.GUESS){
                    FreemuxletData\$SNG.BEST.GUESS[FreemuxletData\$SNG.BEST.GUESS == k] <- Assigned_identity\$Identity[Assigned_identity\$SNG.BEST.GUESS == k]}
                MultiPlexData <- FreemuxletData
                MultiPlexData\$BARCODE <- paste0(files[i],"_",MultiPlexData\$BARCODE)
                if ("${params.Sample_metadata}" != ""){
                    MultiPlexData <- merge(MultiPlexData, meta, by.x = "SNG.BEST.GUESS", by.y = "Sample")
                }
                rownames(MultiPlexData) <- MultiPlexData\$BARCODE
                MultiPlexData\$BARCODE <- NULL
                datalist[[i]] <- AddMetaData(datalist[[i]], metadata = MultiPlexData)
                Multiplex <- "Si"
        }
        if (Multiplex == "No"){
            if ("${params.Sample_metadata}" != ""){
                for (k in colnames(meta)){
                    datalist[[i]] <- AddMetaData(datalist[[i]], metadata = meta[meta\$Sample == fnames[i],][[k]], col.name = k)
                }
            }
        }
        print(paste(i, 'out', length(files)))}
    datalist <- Reduce(merge, datalist)
    if (file.exists(paste0(files[i],".Freemuxlet.samples"))){
        best_guess_factor <- factor(datalist@meta.data\$SNG.BEST.GUESS)
        best_guess_numbers <- as.numeric(best_guess_factor)
        datalist@meta.data\$SNG.BEST.GUESS <- LETTERS[best_guess_numbers]
        datalist@meta.data\$SNG.BEST.GUESS <- paste0(datalist@meta.data\$orig.ident,"-",datalist@meta.data\$SNG.BEST.GUESS)
    }
    saveRDS(datalist, 'datalist.rds')
    """
}
process Seurat_QC_integration {
    publishDir "${params.outdir}/4-Seurat_QC-Integration", mode:'copy'
    
    input:
    path(Datalist)
    path(Stats)

    output:    
    path("datalist.postQC-int*.rds"), emit: rds
    path("Analysis/") 
    path("StatsTable.csv"), emit: Stats
    path("DataForCelltypist/"), emit: DataForCelltypist, optional: true

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
    library(Seurat)
    library(preprocessCore)
    library(dplyr)
    library(ggplot2) 
    library(harmony)
    })
    colors <- clusterExperiment::bigPalette
    colors <- colors[colors!="black"]
    set.seed(123)
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
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(colour = NA, fill = NA))
    }
    Stats <- list()
    for (CSV in list.files("./",pattern="*.csv")) {
    print(CSV)
    Stats[[CSV]] <- read.csv(CSV)
    }
    Stats <- Reduce(rbind,Stats)

    system("mkdir -p Analysis/")
    system("mkdir -p Analysis/Images/")
    system("mkdir -p Analysis/Images/1-QC/")
    system("mkdir -p Analysis/Images/2-PreIntegration/")
    datalist <- readRDS("${Datalist}") %>% JoinLayers()
    datalist[["percent.mt"]] <- PercentageFeatureSet(datalist, pattern="^MT-")
    VlnPlot(datalist, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,group.by = "${params.Batch}",raster=FALSE)
    ggsave(paste0("Analysis/Images/1-QC/Pre-QC.byBatch.png"), width = 16, height = 9)
    datalist <- datalist %>% NormalizeData(verbose = FALSE) %>% subset(subset = nFeature_RNA > 200 & percent.mt < 15)
    datalist <- datalist %>% AddModuleScore(features = list(hk), name = "housekeeping") %>% AddModuleScore(features = list(ribo), name = "ribo.percent")
    print("Splitting for QC per sample")

    PostQC <- as.data.frame(table(datalist@meta.data\$orig.ident))
    colnames(PostQC) <- c("RunID","PostQC")
    Stats <- merge(Stats,PostQC,by="RunID")
    Stats <- Stats[Stats\$PostQC > 0,]

    if (any(colnames(datalist@meta.data) == "SNG.BEST.GUESS")){
        PostQC <- as.data.frame(table(datalist@meta.data\$orig.ident,datalist@meta.data\$SNG.BEST.GUESS))
        colnames(PostQC) <- c("RunID","SNG.BEST.GUESS","PostQC2")
        Stats <- merge(Stats,PostQC,by="RunID")
        Stats <- Stats[Stats\$PostQC2 > 0,]
    } 
    write.csv(Stats, "StatsTable.csv", row.names = F)
    gc()
    datalist <- NormalizeData(datalist, verbose = F)
    VlnPlot(datalist, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,group.by = "${params.Batch}",raster=FALSE)
    ggsave(paste0("Analysis/Images/1-QC/QC.byBatch.png"), width = 16, height = 9)
    ###################################################################################################################################################
    datalist <- FindVariableFeatures(datalist, selection.method = "vst", nfeatures = 2000, verbose = F)
    LabelPoints(plot=VariableFeaturePlot(datalist), points = head(VariableFeatures(datalist), 15), repel=T, xnudge=0, ynudge=0) + theme(legend.position="none") + seurat_theme()
    ggsave(paste0("Analysis/Images/1-QC/VariableFeatures.png"), width = 9, height = 9)
    #all.genes <- rownames(datalist)
    datalist <- ScaleData(datalist, verbose = F) # ,features = all.genes
    datalist <- CellCycleScoring(datalist,s.features = cc.genes\$s.genes, g2m.features = cc.genes\$g2m.genes, set.ident = TRUE, search = TRUE)
    gc()
    datalist <- ScaleData(datalist,vars.to.regress = c('percent.mt')) #, 'housekeeping1', 'ribo.percent1', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'
    datalist <- RunPCA(datalist, verbose = F)

    datalist <- SetIdent(datalist,value = "orig.ident")
    DimPlot(datalist, reduction = "pca",raster=FALSE)
    ggsave(paste0("Analysis/Images/1-QC/PCA.png"), width = 9, height = 9)
    Best_PC <- optimizePCA(datalist, 0.9) # 0.9 is % of variance explained by these PCs
    ElbowPlot(datalist, ndims = length(datalist@reductions\$pca), reduction = "pca") + geom_vline(xintercept=Best_PC, linetype = "dashed", color = "red")
    ggsave(paste0("Analysis/Images/1-QC/ElbowPlot.png"), width = 9, height = 9, bg="white")
    datalist <- RunUMAP(datalist, reduction = "pca", dims = 1:Best_PC, verbose = F, seed.use = 123)
    ###################################################################################################################################################
    if ("${params.Batch}" != "orig.ident"){
        DimPlot(datalist,reduction = "umap", group.by = "${params.Batch}", cols = colors,raster=FALSE) + seurat_theme()
        ggsave(paste0("Analysis/Images/2-PreIntegration/Batch.png"), width = 9, height = 9)
    } else {
        DimPlot(datalist,reduction = "umap", group.by = "orig.ident", cols = colors,raster=FALSE) + seurat_theme()
        ggsave(paste0("Analysis/Images/2-PreIntegration/Orig.ident.png"), width = 9, height = 9)
    }
    DimPlot(datalist, group.by = 'Phase', reduction = 'umap',raster=FALSE) + seurat_theme()
    ggsave("Analysis/Images/2-PreIntegration/Phase.png", width = 9, height = 9)
    datalist %>% FeaturePlot(reduction = "umap",features="nCount_RNA",raster=FALSE) + seurat_theme()
    ggsave(paste0("Analysis/Images/2-PreIntegration/nCount_RNA.png"), width = 9, height = 9)
    datalist %>% FeaturePlot(reduction = "umap",features="ribo.percent1",raster=FALSE) + seurat_theme()
    ggsave(paste0("Analysis/Images/2-PreIntegration/Ribo.percent1.png"), width = 9, height = 9)
    datalist %>% FeaturePlot(reduction = "umap",features="percent.mt",raster=FALSE) + seurat_theme()
    ggsave(paste0("Analysis/Images/2-PreIntegration/Percent.mt.png"), width = 9, height = 9)
    if (any(colnames(datalist@meta.data) == "SNG.BEST.GUESS")){
        if (length(unique(datalist@meta.data\$SNG.BEST.GUESS)) <= length(colors)){
            DimPlot(datalist, group.by = 'SNG.BEST.GUESS', reduction = 'umap', cols = colors,raster=FALSE) + seurat_theme()
            ggsave(paste0("Analysis/Images/2-PreIntegration/Demultiplex.png"), width = 9, height = 9)
        }
    }
    ###################################################################################################################################################
    if (length(unique(datalist@meta.data\$orig.ident)) == 1){
        if ("SNG.BEST.GUESS" %in% colnames(datalist@meta.data)){
            system("mkdir -p Analysis/Images/3-PosIntegration/")
            datalist <- datalist %>% RunHarmony("SNG.BEST.GUESS", plot_convergence = T, max_iter = 50)
            datalist <- datalist %>% 
            RunUMAP(reduction = "harmony", verbose = F, dims = 1:Best_PC, seed.use = 123) %>% 
            FindNeighbors(reduction = "harmony", dims = 1:Best_PC) %>% 
            FindClusters(resolution = seq(0.4,1.4,0.2)) %>% 
            identity()
            FeaturePlot(datalist,reduction = "umap",features="nCount_RNA",raster=FALSE) + seurat_theme()
            ggsave(paste0("Analysis/Images/3-PosIntegration/nCount_RNA.png"), width = 9, height = 9)
            FeaturePlot(datalist, reduction = "umap",features="ribo.percent1",raster=FALSE) + seurat_theme()
            ggsave(paste0("Analysis/Images/3-PosIntegration/Ribo.percent1.png"), width = 9, height = 9)
            DimPlot(datalist, group.by = 'Phase', reduction = 'umap',raster=FALSE) + seurat_theme()
            ggsave("Analysis/Images/3-PosIntegration/Phase.png", width = 9, height = 9)
            datalist %>% FeaturePlot(reduction = "umap",features="percent.mt",raster=FALSE) + seurat_theme()
            ggsave(paste0("Analysis/Images/3-PosIntegration/Percent.mt.png"), width = 9, height = 9)
            if (length(unique(datalist@meta.data\$SNG.BEST.GUESS)) <= length(colors)){
                DimPlot(datalist, group.by = 'SNG.BEST.GUESS', reduction = 'umap', cols = colors,raster=FALSE) + seurat_theme()
                ggsave(paste0("Analysis/Images/3-PosIntegration/Demultiplex.png"), width = 9, height = 9)
            }
            colnames(datalist@meta.data) <- gsub('RNA_snn', 'Seurat_clusters', colnames(datalist@meta.data))
        } else {
            datalist <- datalist %>% 
            FindNeighbors(reduction = "pca", dims = 1:Best_PC) %>% 
            FindClusters(resolution = seq(0.4,1.4,0.2)) %>% 
            identity()
            colnames(datalist@meta.data) <- gsub('RNA_snn', 'Seurat_clusters', colnames(datalist@meta.data))
        }
    } else {
        system("mkdir -p Analysis/Images/3-PosIntegration/")
        if ("${params.Batch}" != "orig.ident"){
            datalist <- datalist %>% RunHarmony(c("${params.Batch}","orig.ident"), plot_convergence = T, max_iter = 50)
        } else {
            datalist <- datalist %>% RunHarmony("orig.ident", plot_convergence = T, max_iter = 50)
        }
        datalist <- datalist %>% 
        RunUMAP(reduction = "harmony", verbose = F, dims = 1:Best_PC, seed.use = 123) %>% 
        FindNeighbors(reduction = "harmony", dims = 1:Best_PC) %>% 
        FindClusters(resolution = seq(0.4,1.4,0.2)) %>% 
        identity()
        colnames(datalist@meta.data) <- gsub('RNA_snn', 'Seurat_clusters', colnames(datalist@meta.data))
        ###################################################################################################################################################
        if ("${params.Batch}" != "orig.ident"){
            DimPlot(datalist,reduction = "umap", group.by = "${params.Batch}", cols = colors,raster=FALSE) + seurat_theme()
            ggsave(paste0("Analysis/Images/3-PosIntegration/Batch.png"), width = 9, height = 9)
        } else {
            DimPlot(datalist,reduction = "umap", group.by = "orig.ident", cols = colors,raster=FALSE) + seurat_theme()
            ggsave(paste0("Analysis/Images/3-PosIntegration/Orig.ident.png"), width = 9, height = 9)
        }
        FeaturePlot(datalist,reduction = "umap",features="nCount_RNA",raster=FALSE) + seurat_theme()
        ggsave(paste0("Analysis/Images/3-PosIntegration/nCount_RNA.png"), width = 9, height = 9)
        FeaturePlot(datalist, reduction = "umap",features="ribo.percent1",raster=FALSE) + seurat_theme()
        ggsave(paste0("Analysis/Images/3-PosIntegration/Ribo.percent1.png"), width = 9, height = 9)
        DimPlot(datalist, group.by = 'Phase', reduction = 'umap',raster=FALSE) + seurat_theme()
        ggsave("Analysis/Images/3-PosIntegration/Phase.png", width = 9, height = 9)
        datalist %>% FeaturePlot(reduction = "umap",features="percent.mt",raster=FALSE) + seurat_theme()
        ggsave(paste0("Analysis/Images/3-PosIntegration/Percent.mt.png"), width = 9, height = 9)
        if (any(colnames(datalist@meta.data) == "SNG.BEST.GUESS")){
            if (length(unique(datalist@meta.data\$SNG.BEST.GUESS)) <= length(colors)){
                DimPlot(datalist, group.by = 'SNG.BEST.GUESS', reduction = 'umap', cols = colors,raster=FALSE) + seurat_theme()
                ggsave(paste0("Analysis/Images/3-PosIntegration/Demultiplex.png"), width = 9, height = 9)
            }
        }
    }
    saveRDS(datalist, 'datalist.postQC-int.rds')
    if ("${params.ModelsCelltypist}" != ""){
        datalist_sce <- as.SingleCellExperiment(datalist)
        DropletUtils:::write10xCounts(x = datalist_sce@assays@data\$counts, path = paste0("DataForCelltypist"),overwrite = TRUE)
        write.table(datalist@meta.data,file = paste0("DataForCelltypist/Metadata.csv"), row.names = T, quote = F, col.names = T, sep = "\\t")
        write.table(datalist@reductions\$umap@cell.embeddings,file = paste0("DataForCelltypist/Cell.embeddings.UMAP.csv"), row.names = T, quote = F, col.names = T, sep = "\\t")
        write.table(datalist@reductions\$harmony@cell.embeddings,file = paste0("DataForCelltypist/Cell.embeddings.PCA.csv"), row.names = T, quote = F, col.names = T, sep = "\\t")
    }
    """
}
process CellTypist {
    publishDir "${params.outdir}/5-Seurat-Cell_Annotation", mode:'copy'
    
    input:
    path(Datalist)

    output:    
    path("StatsCelltypist/")
    path("*_Celltypist_annotation.csv"), emit: Celltypist_annotations

    script:
    """
    #!/usr/bin/env python
    import celltypist
    from celltypist import models
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import os

    def Celltypist_annotation(Dir):
        os.system("mkdir -p StatsCelltypist")
        models.download_models(model = "${params.ModelsCelltypist}".replace(" ","").split(","))
        adata = sc.read_10x_mtx(Dir)
        metadata = pd.read_csv( Dir + 'Metadata.csv',sep="\\t",low_memory=False)
        adata.obs = adata.obs.join(metadata)
        umap = pd.read_csv( Dir + 'Cell.embeddings.UMAP.csv',sep="\\t",low_memory=False)
        pca = pd.read_csv( Dir + 'Cell.embeddings.PCA.csv',sep="\\t",low_memory=False)
        adata.obsm['X_umap'] = umap.to_numpy()
        adata.obsm['X_pca'] = pca.to_numpy()
        clusters = adata.obs.columns[adata.obs.columns.str.contains("^Seurat_clusters_res")].tolist()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        for ModelFile in "${params.ModelsCelltypist}".replace(" ","").split(","):
            ModelName = ModelFile.replace(".pkl","")
            print(ModelName)
            Model = pd.DataFrame()
            for cluster in clusters:
                adata.obs[cluster] = adata.obs[cluster].astype('category')
                sc.pl.umap(adata, color = cluster)  
                predictions = celltypist.annotate(adata, model = ModelFile, majority_voting = True, over_clustering = cluster)
                Model[ModelName + "_"+cluster.replace("Seurat_clusters_res","res")] = predictions.predicted_labels["majority_voting"]
                predictions.predicted_labels.to_csv("StatsCelltypist/" + ModelName + '_predicted_labels-'+cluster+'.csv')
                predictions.probability_matrix.to_csv("StatsCelltypist/" + ModelName + '_probability_matrix-'+cluster+'.csv')
            Model.to_csv(ModelName +'_Celltypist_annotation.csv')

    Celltypist_annotation("DataForCelltypist/")
    """
}
process RunAzimuth {
    publishDir "${params.outdir}/5-Seurat-Cell_Annotation", mode:'copy'
    errorStrategy "ignore"

    input:
    path(Datalist)

    output:    
    path("Azimuth/")
    path("*.png")

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(Azimuth)
    library(ggplot2)
    library(patchwork)
    colors <- clusterExperiment::bigPalette
    options(timeout=100000000)
    options(future.globals.maxSize = 4 * 1024^3)  # 4 GB
    seurat_theme <- function(){
        theme_bw() +
            theme(panel.background = element_rect(colour = "black", size=0.1), 
                plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
                axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(colour = NA, fill = NA))
        }
    datalist <- RunAzimuth("${Datalist}", reference = "pbmcref")
    DimPlot(datalist, reduction = "ref.umap", group.by = c("seurat_clusters","predicted.celltype.l2"),label=T,repel=T,raster=FALSE, cols = colors) + seurat_theme()
    ggsave("Azimuth.png",width = 24, height = 9)

    AzimuthAnnotation <- as.data.frame(datalist@meta.data[(grepl("predicted.celltype.l",colnames(datalist@meta.data))) & !(grepl("score",colnames(datalist@meta.data)))])
    rm(datalist)
    gc()
    Sobj <- readRDS("${Datalist}")
    Sobj <- AddMetaData(Sobj, AzimuthAnnotation)
    system("mkdir -p Azimuth")
    for (Annotation in names(Sobj@meta.data)[grepl("predicted.celltype.",names(Sobj@meta.data))]){
        print(Annotation)
        DimPlot(Sobj,reduction = "umap", group.by = Annotation, cols = colors,raster=FALSE) + seurat_theme()
        ggsave(paste0("Azimuth/Azimuth.",Annotation,".png"),width = 9, height = 9)
    }
    """
}
process Seurat_Cell_Annotation {
    publishDir "${params.outdir}/5-Seurat-Cell_Annotation", mode:'copy'
    
    input:
    path(Datalist)
    file(ModelsCelltypist)

    output:    
    path("datalist.Cell_Annotation*.rds"), emit: rds 
    path("*.png")
    path("Seurat_clusters/")
    path("SingleR/")
    path("Celltypist/"), optional: true

    script:
    """
    #!/usr/bin/env Rscript
    options(timeout=10000000000000000000000000)
    suppressPackageStartupMessages({
    library(rols)})
    cl <- rols::Ontology("cl")
    suppressPackageStartupMessages({
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
    })
    ## Run just one time
    ref1 <- celldex::BlueprintEncodeData()
    ref2 <- celldex::DatabaseImmuneCellExpressionData()
    ref3 <- celldex::HumanPrimaryCellAtlasData()
    ref4 <- celldex::MonacoImmuneData()
    shared <- Reduce(intersect, list(rownames(ref1), rownames(ref2), rownames(ref3), rownames(ref4)))
    ref1 <- ref1[shared,]; ref2 <- ref2[shared,]; ref3 <- ref3[shared,]; ref4 <- ref4[shared,]
    #save(ref1, ref2, ref3, ref4, cl, file = 'SingleR_refs.RData')
    #load("Data/SingleR_refs.RData")
    set.seed(123)
    colors <- c("#E31A1C","#1F78B4","#33A02C","#FF7F00","#6A3D9A","#B15928","#A6CEE3","#bd18ea","cyan","#B2DF8A","#FB9A99","deeppink4","#00B3FFFF","#CAB2D6","#FFFF99","#05188a","#CCFF00FF","cornflowerblue","#f4cc03","black","blueviolet","#4d0776","maroon3","blue","#E5D8BD","cadetblue4","#e5a25a","lightblue1","#F781BF","#FC8D62","#8DA0CB","#E78AC3","green3","#E7298A","burlywood3","#A6D854","firebrick","#FFFFCC","mediumpurple","#1B9E77","#FFD92F","deepskyblue4","yellow3","#00FFB2FF","#FDBF6F","#FDCDAC","gold3","#F4CAE4","#E6F5C9","#FF00E6FF","#7570B3","goldenrod","#85848f","lightpink3","olivedrab","cadetblue3")
    seurat_theme <- function(){
        theme_bw() +
            theme(panel.background = element_rect(colour = "black", size=0.1), 
                plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
                axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(colour = NA, fill = NA))
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
    CellAnnotations <- function(RDS){
        print(paste0("Working in ",RDS))
        Sobj <- readRDS(RDS)
        system(paste0("mkdir -p Seurat_clusters"))
        system(paste0("mkdir -p SingleR"))
        clusters <- grep('Seurat_clusters', names(Sobj@meta.data), value = T)
        for (cluster in clusters) {
            Sobj <- ppSingleR(Sobj, cluster)
            print(cluster)
            DimPlot(Sobj,reduction = "umap", group.by = cluster, cols = colors,raster=FALSE) + seurat_theme()
            ggsave(paste0("Seurat_clusters/",cluster,".png"),heigh=9,width=9)
        }
        for (Cluster in names(Sobj@meta.data)[grepl("SingleR_res.",names(Sobj@meta.data))]){
            print(Cluster)
            DimPlot(Sobj,reduction = "umap", group.by = Cluster, cols = colors,raster=FALSE) + seurat_theme()
            ggsave(paste0("SingleR/SingleR.",Cluster,".png"),width = 9, height = 9)
        }
        ###################################################################################################################################################
        if ("${params.ModelsCelltypist}" != ""){
            system(paste0("mkdir -p Celltypist"))
            for (CelltypistAnnot in list.files("./",pattern="Celltypist_annotation.csv")){
                CelltypistModel <- read.delim(CelltypistAnnot,sep = ",",row.names=1)
                Sobj <- AddMetaData(Sobj, CelltypistModel)
                ModelName <- gsub("_Celltypist_annotation.csv","",CelltypistAnnot)
                for (Cluster in names(Sobj@meta.data)[grepl(paste0(ModelName,"_res"),names(Sobj@meta.data))]){
                    print(Cluster)
                    DimPlot(Sobj,reduction = "umap", group.by = Cluster, cols = colors,raster=FALSE) + seurat_theme()
                    ggsave(paste0("Celltypist/Celltypist.",Cluster,".png"),width = 9, height = 9)
                }
                clustree(Sobj,prefix=paste0(ModelName,"_res."),node_text_angle = 15)
                ggsave(paste0("Clustree_Celltypist_",ModelName,".png"), width = 21,height = 9)
            }
        }
        ###################################################################################################################################################
        clustree(Sobj,prefix="SingleR_res.",node_text_angle = 15)
        ggsave(paste0("Clustree_SingleR.png"), width = 21,height = 9)
        clustree(Sobj,prefix="Seurat_clusters_res.",node_text_angle = 15)
        ggsave(paste0("Clustree_Seurat_clusters.png"), width = 21,height = 9)
        ###################################################################################################################################################
        return(Sobj)
    }
    datalist <- CellAnnotations("${Datalist}")
    if ("${params.ModelsCelltypist}" != ""){
        Idents(datalist) <- "Immune_All_Low_res.0.8"
    } else {
        Idents(datalist) <- "SingleR_res.0.8"
    }
    saveRDS(datalist, file = "datalist.Cell_Annotation.rds")
    """
}
process CellChat {
    publishDir "${params.outdir}/6-Cell_Communication", mode:'copy'

    input:
    path(Datalist)

    output:    
    path("*.pdf")

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(tibble)
    library(data.table)
    library(circlize)
    library(pheatmap)
    library(reshape2)
    library(ComplexHeatmap)
    library(CellChat)
    library(reticulate)})
    set.seed(123)
    data <- readRDS('${Datalist}')

    cellChat <- createCellChat(object = data, group.by = "ident", assay = "RNA")
    ##############################################################################################
    #Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)
    # use a subset of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    cellChat@DB <- CellChatDB.use
    ##############################################################################################
    cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
    cellChat <- identifyOverExpressedGenes(cellChat)
    cellChat <- identifyOverExpressedInteractions(cellChat)
    cellChat <- smoothData(cellChat, adj = PPI.human)
    cellChat <- computeCommunProb(cellChat, type = "triMean", raw.use = F) # LR pairs
    ##############################################################################################
    # Infer the cell-cell communication at a signaling pathway level
    cellChat <- computeCommunProbPathway(object = cellChat) # Pathways
    # Calculate the aggregated cell-cell communication network
    #CellChat calculates the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. Users can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
    cellChat <- aggregateNet(cellChat)
    ##############################################################################################
    groupSize <- as.numeric(table(cellChat@idents))
    pdf(paste0("InteractionNumber.pdf"))
    netVisual_circle(cellChat@net\$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste0("Number of interactions"))
    dev.off()
    pdf(paste0("InteractionStrength.pdf"))
    netVisual_circle(cellChat@net\$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste0("Interaction weights/strength"))
    dev.off()

    groupSize <- as.numeric(table(cellChat@idents))
    mat <- cellChat@net\$weight
    par(mfrow = c(3,4), xpd=TRUE)
    pdf(paste0("InteractionStrength.perCellType.pdf"))
    for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]

    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()

    par(mfrow=c(1,1))
    pdf(paste0("Heatmap.pdf"))
    netVisual_heatmap(cellChat, color.heatmap = "Reds")
    dev.off()
    """
}
process Seurat_create_object_mouse {
    publishDir "${params.outdir}/Seurat", mode:'copy'
    
    input:
    path(reads)
    file(demultiplexData)
    file(metadata)

    output:     
    path("datalist.rds"), emit: rds 
    path("datalist.h5ad"), emit: h5ad 

    script:
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
    to.remove <- read.table(doubletFile, sep = '\\t', header = T) %>% 
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
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(colour = NA, fill = NA))
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
    write.table(results, file.path(paste0(fnames[i], '.txt')), sep = '\\t', quote = F, col.names = T, row.names = T)
    }

    fnames <- basename(list.dirs('./', recursive = F))
    files <- file.path(fnames)

    datalist <- list()
    if ("${params.Sample_metadata}" != ""){
        if(grepl(".tsv","${params.Sample_metadata}")){
        meta <- read.delim(basename("${params.Sample_metadata}"),sep="\\t")}
        if(grepl(".tab","${params.Sample_metadata}")){
        meta <- read.delim(basename("${params.Sample_metadata}"),sep="\\t")}
        if(grepl(".csv","${params.Sample_metadata}")){
        meta <- read.delim(basename("${params.Sample_metadata}"),sep=",")}
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
workflow Demultiplex {
    take:
        data
    main:
        if (params.Freemuxlet != false & params.VCF_Files == ""){
            print "Demultiplex activado ${params.NPlex}Plex"
            Freemuxlet(data)
            Demultiplex_ch = Freemuxlet.out.map{it[1]}.collect()
        } else {
            if (params.Freemuxlet != false & params.VCF_Files != ""){
                VCFs = channel.fromFilePairs(params.VCF_Files, size: -1)
                CombinedFremuxletDemuxlet(data.combine(VCFs, by: 0))
                Demultiplex_ch = CombinedFremuxletDemuxlet.out.map{it[1]}.collect()
            } else {
                if (params.VCF_Files != ""){
                    VCFs = channel.fromFilePairs(params.VCF_Files, size: -1)
                    Demuxlet(data.combine(VCFs, by: 0))
                    Demultiplex_ch = Demuxlet.out.map{it[1]}.collect()
                } else {
                    Demultiplex_ch = channel.of("None").combine(data.map{it[0]}).map{it[0]}
                }
            }
        }
    emit:
        Demultiplex_ch
}
workflow ObjectCreation {
    take:
        Clean_outs
        DemultiplexData
    main:
        if (params.Sample_metadata != ""){
            MetaDataFiles_ch = channel.fromPath(params.Sample_metadata)
            MetaDataFiles_ch.view()
        } else {
            MetaDataFiles_ch = channel.of("None")
        }
        if (params.Mouse != false){ 
        Seurat_create_object_mouse(Clean_outs.collect(),DemultiplexData,MetaDataFiles_ch)
        scObject = Seurat_create_object_mouse.out.rds
        } else {
        Seurat_Object_creation(Clean_outs.collect(),DemultiplexData,MetaDataFiles_ch)
        scObject = Seurat_Object_creation.out.rds}
        emit:
            scObject
}
workflow CellAnnotation {
    take:
        data
        data2
    main:
        if (params.ModelsCelltypist != ""){
            CellTypist(data2)
            CellTypistData = CellTypist.out.Celltypist_annotations
        } else {
            CellTypistData = channel.of("None")
        }
        Seurat_Cell_Annotation(data,CellTypistData)
        RunAzimuth(data)
    emit:
        Seurat_Cell_Annotation.out.rds
}
workflow Mapping { 
    take:
        data
        data2
    main: 
        Cellranger(data,data2)
    emit: 
        Cellranger.out
}
workflow {
    if (params.help == true){
        print '\nBasic Usage:\n\tnextflow {Path_to_script}/main.nf --dataDir="{Path_to_Data}/RawData/*" --Project_Name="Chile_03" -c {Path_to_script}/nextflow.config --NPlex=4 -resume'
        print 'Parameters:\n\tMandatories:\n\t\t--dataDir = Sets the input sample table file.\n\t\t--ref_data = Directory with the organism genome processed by cellranger tools.\n\t\t-c = Sets the configuration file.'
        print '\n\tAlternatives:\n\t\t--Project_Name = Sets a custom directory for the project.\n\t\t--Sample_metadata = File with metadata of samples, must be a Tabular delimited file (TSV|TAB),\n\t\t\tand must have a column calles "Sample" matching the names of raw data folders.\n\t\t--Batch = For Seurat analysis, column with the batches. (def. "orig.ident")\n\t\t--Freemuxlet = Set use of Freemuxlet to true for demultiplexing without genotype. (def. false)\n\t\t--ModelsCelltypist = Models for use of celltypist, . (def. "", no celltypist annotation)\n\t\t--VCF_Files = Path to libraries VCFs, to use for genotyping guided demultiplexing. (def. "", no genotyping guided demultiplexing)\n\t\t--MainOutdir = Name of directory to save outputs. (def. "LatinCells_Results/${params.Project_Name}/")\n\t\t--NPlex = Sets the number of samples per run in multiplexed libraries.\n\t\t-resume = Continue an analysis after the last well executed part of a workflow,\n\t\t\tbefore a error or a change made in thee code of a process.[Nextflow parameter].\n\t\t--help = Prints this help message.'
        print 'Output directory:\n\tBy default this workflow creates a directory called "LatinCells_Results/General/" in the working directory.'
    } else{   
        dataDir_ch = channel.fromFilePairs(params.dataDir, type: 'dir', size: -1) 
        ref_data_ch = channel.fromPath(params.ref_data, type: 'dir') 
        out_dir = file(params.outdir)
        out_dir.mkdir() 
        print "Proyecto: ${params.Project_Name}"
        print "Directorio de Salida: ${params.outdir}"
        print "Muestras en: ${params.dataDir}"
        Mapping(dataDir_ch,ref_data_ch)
        MappingStats(Mapping.out.map{it[1]}.collect())
        Demultiplex(Mapping.out)
        SoupX_scDblFinder(Mapping.out,Demultiplex.out)
        ObjectCreation(SoupX_scDblFinder.out.clean_reads,Demultiplex.out)
        Seurat_QC_integration(ObjectCreation.out, SoupX_scDblFinder.out.Stats.collect())
        CellAnnotation(Seurat_QC_integration.out.rds,Seurat_QC_integration.out.DataForCelltypist)
        CellChat(CellAnnotation.out)
        }
}
workflow NoIntegration { 
    dataDir_ch = channel.fromFilePairs(params.dataDir, type: 'dir', size: -1) 
    ref_data_ch = channel.fromPath(params.ref_data, type: 'dir') 
    out_dir = file(params.outdir)
    out_dir.mkdir() 
    print "Proyecto: ${params.Project_Name}"
    print "Directorio de Salida: ${params.outdir}"
    print "Muestras en: ${params.dataDir}"
    Mapping(dataDir_ch,ref_data_ch)
    MappingStats(Mapping.out.map{it[1]}.collect())
    Demultiplex(Mapping.out)
    SoupX_scDblFinder(Mapping.out,Demultiplex.out)
    ObjectCreation(SoupX_scDblFinder.out.clean_reads,Demultiplex.out)
}
workflow JustMapping {
    dataDir_ch = channel.fromFilePairs(params.dataDir, type: 'dir', size: -1) 
    ref_data_ch = channel.fromPath(params.ref_data, type: 'dir') 
    out_dir = file(params.outdir)
    out_dir.mkdir() 
    print "Proyecto: ${params.Project_Name}"
    print "Directorio de Salida: ${params.outdir}"
    print "Muestras en: ${params.dataDir}"
    Mapping(dataDir_ch,ref_data_ch)
    MappingStats(Mapping.out.map{it[1]}.collect())
}
workflow prueba {    
}