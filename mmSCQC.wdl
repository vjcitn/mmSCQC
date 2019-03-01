
workflow task1 {
  call dommSCQC { }
}

task dommSCQC {
  command {
## ------------------------------------------------------------------------
      R -e "ii = rownames(installed.packages()); \
      condinst = function(x,ins) if (!(x %in% ins)) BiocManager::install(x); \
      condinst('BiocFileCache',ii); \
      condinst('SingleCellExperiment',ii); \
      condinst('org.Mm.eg.db',ii); \
      condinst('scater',ii); \
      condinst('TxDb.Mmusculus.UCSC.mm10.ensGene',ii); \
      library('BiocFileCache'); \
      library('SingleCellExperiment'); \
      library('org.Mm.eg.db'); \
      library('scater'); \
      library('TxDb.Mmusculus.UCSC.mm10.ensGene'); \
      bfc <- BiocFileCache('raw_data', ask = FALSE); \
      lun.zip <- bfcrpath(bfc,  \
          file.path('https://www.ebi.ac.uk/arrayexpress/files', \
              'E-MTAB-5522/E-MTAB-5522.processed.1.zip')); \
      lun.sdrf <- bfcrpath(bfc,  \
          file.path('https://www.ebi.ac.uk/arrayexpress/files', \
              'E-MTAB-5522/E-MTAB-5522.sdrf.txt')); \
      unzip(lun.zip, exdir=tempdir()); \
      plate1 <- read.delim(file.path(tempdir(), 'counts_Calero_20160113.tsv'),  \
          header=TRUE, row.names=1, check.names=FALSE); \
      plate2 <- read.delim(file.path(tempdir(), 'counts_Calero_20160325.tsv'),  \
          header=TRUE, row.names=1, check.names=FALSE); \
      gene.lengths <- plate1[['Length']] ; \
      plate1 <- as.matrix(plate1[,-1]) ; \
      plate2 <- as.matrix(plate2[,-1]);  \
      rbind(Plate1=dim(plate1), Plate2=dim(plate2)); \
      stopifnot(identical(rownames(plate1), rownames(plate2))); \
      all.counts <- cbind(plate1, plate2); \
      sce <- SingleCellExperiment(list(counts=all.counts)); \
      rowData(sce)[['GeneLength']] <- gene.lengths; \
      isSpike(sce, 'ERCC') <- grepl('^ERCC', rownames(sce)); \
      summary(isSpike(sce, 'ERCC')); \
      is.sirv <- grepl('^SIRV', rownames(sce)); \
      sce <- sce[!is.sirv,];  \
      summary(is.sirv); \
      metadata <- read.delim(lun.sdrf, check.names=FALSE, header=TRUE); \
      m <- match(colnames(sce), metadata[['Source Name']]);
      stopifnot(all(!is.na(m))); 
      metadata <- metadata[m,]; \
      head(colnames(metadata)); \
      colData(sce)[['Plate']] <- factor(metadata[['Factor Value[block]']]); \
      pheno <- metadata[['Factor Value[phenotype]']]; \
      levels(pheno) <- c('induced', 'control'); \
      colData(sce)[['Oncogene']] <- pheno; \
      table(colData(sce)[['Oncogene']], colData(sce)[['Plate']]); \
      symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype='ENSEMBL', column='SYMBOL'); \
      rowData(sce)[['ENSEMBL']] <- rownames(sce); \
      rowData(sce)[['SYMBOL']] <- symb; \
      head(rowData(sce)); \
      rownames(sce) <- uniquifyFeatureNames(rowData(sce)[['ENSEMBL']], rowData(sce)[['SYMBOL']]); \
      head(rownames(sce)); \
      location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)[['ENSEMBL']],  \
          column='CDSCHROM', keytype='GENEID'); \
      rowData(sce)[['CHR']] <- location; \
      summary(location=='chrM'); \
      mito <- which(rowData(sce)[['CHR']]=='chrM'); \
      sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito)); \
      head(colnames(colData(sce)), 10); \
      sce[['PlateOnco']] <- paste0(sce[['Oncogene']], '.', sce[['Plate']]); \
      multiplot( \
          plotColData(sce, y='total_counts', x='PlateOnco'), \
          plotColData(sce, y='total_features_by_counts', x='PlateOnco'), \
          plotColData(sce, y='pct_counts_ERCC', x='PlateOnco'), \
          plotColData(sce, y='pct_counts_Mt', x='PlateOnco'), \
          cols=2); \
       \
      par(mfrow=c(1,3)); \
      plot(sce[['total_features_by_counts']], sce[['total_counts']]/1e6, xlab='Number of expressed genes', \
          ylab='Library size (millions)'); \
      plot(sce[['total_features_by_counts']], sce[['pct_counts_ERCC']], xlab='Number of expressed genes', \
          ylab='ERCC proportion (%)'); \
      plot(sce[['total_features_by_counts']], sce[['pct_counts_Mt']], xlab='Number of expressed genes', \
          ylab='Mitochondrial proportion (%)');"  
   }  
}  
