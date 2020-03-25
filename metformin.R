library( tidyverse )

## Preliminaries
synapser::synLogin()
syn <- synExtra::synDownloader("/data/metformin", ifcollision="overwrite.local")
read_csv2 <- partial( read_csv, col_types=cols() )
read_tsv2 <- partial( read_tsv, col_types=cols() )

## Download counts data and matching metadata
X <- syn("syn21763099") %>% read_tsv2
Y <- syn("syn21781943") %>% read_csv2

## Identify protein-coding genes
GM <- read_csv2( "GRCh38.94.genemap.csv" ) %>%
    filter( gene_biotype == "protein_coding" )

## Standardizes Metformin / Control annotations
metcontrol <- function( x ) {
    case_when(grepl("Control", x) ~ "Control",
              grepl("Met 40", x) ~ "Metformin",
              TRUE ~ as.character(NA)) %>%
        factor( levels=c("Control", "Metformin") )
}

## Identify Metformin (40 uM/72hr) samples and matching controls
Y1 <- Y %>% filter( grepl("72hr", Treatment) ) %>%
    mutate_at( "Treatment", metcontrol ) %>%
    na.omit() %>% select( Sample, Treatment ) %>%
    as.data.frame %>% column_to_rownames("Sample")

## Retrieve the corresponding expression slice
X1 <- X %>% select( symbol, one_of(rownames(Y1)) ) %>%
    filter( symbol %in% GM$gene_name ) %>%
    group_by( symbol ) %>% summarize_all( sum ) %>%
    as.data.frame %>% column_to_rownames("symbol")

## Create the design matrix and estimate dispersion
dl <- edgeR::DGEList( counts = X1, samples = Y1 )
dl <- edgeR::calcNormFactors( dl )
mmx <- model.matrix( ~Treatment, data = dl$samples )
dl <- edgeR::estimateDisp( dl, mmx )

## Compute differential expression and
##   identify Metformin's drug-associated gene list (DGL)
gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
dgl <- edgeR::topTags( gf, nrow(X1) ) %>% as.data.frame %>%
    rownames_to_column( "Gene" ) %>% as_tibble %>%
    filter( FDR < 0.05 ) %>% pull( Gene )

## Write to file
cat( dgl, file="metformin-dgl.txt", sep="\n" )
