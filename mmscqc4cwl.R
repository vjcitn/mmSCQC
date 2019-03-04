## Code is from the Bioconductor simpleSingleCell workflow by Aaron Lun
## ------------------------------------------------------------------------
library(BiocFileCache)
bfc <- BiocFileCache('raw_data', ask = FALSE)
lun.zip <- bfcrpath(bfc, 
    file.path('https://www.ebi.ac.uk/arrayexpress/files',
        'E-MTAB-5522/E-MTAB-5522.processed.1.zip'))
lun.sdrf <- bfcrpath(bfc, 
    file.path('https://www.ebi.ac.uk/arrayexpress/files',
        'E-MTAB-5522/E-MTAB-5522.sdrf.txt'))
unzip(lun.zip, exdir=tempdir())

## ------------------------------------------------------------------------
plate1 <- read.delim(file.path(tempdir(), 'counts_Calero_20160113.tsv'), 
    header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim(file.path(tempdir(), 'counts_Calero_20160325.tsv'), 
    header=TRUE, row.names=1, check.names=FALSE)

gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2))

## ------------------------------------------------------------------------
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)

## ------------------------------------------------------------------------
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths
sce

## ------------------------------------------------------------------------
isSpike(sce, 'ERCC') <- grepl('^ERCC', rownames(sce))
summary(isSpike(sce, 'ERCC'))

## ------------------------------------------------------------------------
is.sirv <- grepl('^SIRV', rownames(sce))
sce <- sce[!is.sirv,] 
summary(is.sirv)

## ------------------------------------------------------------------------
metadata <- read.delim(lun.sdrf, check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[['Source Name']]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))

## ------------------------------------------------------------------------
colData(sce)$Plate <- factor(metadata[['Factor Value[block]']])
pheno <- metadata[['Factor Value[phenotype]']]
levels(pheno) <- c('induced', 'control')
colData(sce)$Oncogene <- pheno
table(colData(sce)$Oncogene, colData(sce)$Plate)

## ------------------------------------------------------------------------
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype='ENSEMBL', column='SYMBOL')
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

## ------------------------------------------------------------------------
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce))

## ------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
    column='CDSCHROM', keytype='GENEID')
rowData(sce)$CHR <- location
summary(location=='chrM')

## ------------------------------------------------------------------------
mito <- which(rowData(sce)$CHR=='chrM')
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)

## ----qcplot416b, fig.wide=TRUE, fig.cap='Distributions of various QC metrics for all cells in the 416B dataset. This includes the library sizes, number of expressed genes, and proportion of reads mapped to spike-in transcripts or mitochondrial genes.'----
sce$PlateOnco <- paste0(sce$Oncogene, '.', sce$Plate)
multiplot(
    plotColData(sce, y='total_counts', x='PlateOnco'),
    plotColData(sce, y='total_features_by_counts', x='PlateOnco'),
    plotColData(sce, y='pct_counts_ERCC', x='PlateOnco'),
    plotColData(sce, y='pct_counts_Mt', x='PlateOnco'),
    cols=2)

## ----qcbiplot416b, fig.width=10, fig.asp=0.5, fig.cap='Behaviour of each QC metric compared to the total number of expressed features. Each point represents a cell in the 416B dataset.'----
par(mfrow=c(1,3))
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab='Number of expressed genes',
    ylab='Library size (millions)')
plot(sce$total_features_by_counts, sce$pct_counts_ERCC, xlab='Number of expressed genes',
    ylab='ERCC proportion (%)')
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab='Number of expressed genes',
    ylab='Mitochondrial proportion (%)')

