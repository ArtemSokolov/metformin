library( tidyverse )

## Preliminaries
synapser::synLogin()
syn <- synExtra::synDownloader("/data/metformin", ifcollision="overwrite.local")
read_csv2 <- partial( read_csv, col_types=cols() )
read_tsv2 <- partial( read_tsv, col_types=cols() )

## Identify protein-coding genes
GM <- read_csv2( "data/GRCh38.94.genemap.csv" ) %>%
    filter( gene_biotype == "protein_coding" )

## Download counts data and matching metadata
X <- syn("syn21763099") %>% read_tsv2 %>% select(-id) %>%
    filter(symbol %in% GM$gene_name) %>%
    group_by(symbol) %>% summarize_all(sum)
Y <- syn("syn21781943") %>% read_csv2 %>%
    select( Sample, Well, Treatment )

## Treatment label patterns
ptr <- c("Metformin"    = "Met 40-[1-3] 72hr",
         "Glyburide"    = "Glyb 40-[1-3] 72hr",
         "Control"      = "Control-[1-3] 72hr",
         "VM"           = "VM \\(lipo\\)",
         "VMpic"        = "VM \\+ pic",
         "TBK1#8"       = "TBK1 #8 lipo",
         "TBK1#10"      = "TBK1 #10 lipo",
         "TBK1#8dsRNA"  = "TBK1 #8 \\+ dsRNA",
         "TBK1#10dsRNA" = "TBK1 #10 \\+ dsRNA")

## Extract sampleIDs associated with each treatment pattern
sid <- map( ptr, ~filter(Y, grepl(.x, Treatment))$Sample )

## Dichotomies of interest:
## 1. Metformin (40uM, 72hr) vs Control (72hr)
## 2. Glyburide (40uM, 72hr) vs Control (72hr)
## 3. Metformin (40uM, 72hr) vs Glyburide (40uM, 72hr)
## 4. VM + pic vs VM (lipo)
## 5. TBK1 #8 (lipo) vs VM (lipo)
## 6. TBK1 #10 (lipo) vs VM (lipo)
## 7. TBK1 #8 + dsRNA vs TBK1 #8 (lipo)
## 8. TBK1 #10 + dsRNA vs TBK1 #10 (lipo)
dich <- list(c("Metformin",   "Control"),
             c("Glyburide",   "Control"),
             c("Metformin",    "Glyburide"),
             c("VMpic",        "VM"),
             c("TBK1#8",       "VM"),
             c("TBK1#10",      "VM"),
             c("TBK1#8dsRNA",  "TBK1#8"),
             c("TBK1#10dsRNA", "TBK1#10"))

## Composes a meta-data matrix based on a dichotomy definition
## smap - mapping of sample IDs to wells
## trt  - label corresponding to the "treatment"
## ctrl - label corresponding to the "control"
makeMeta <- function(smap, trt, ctrl) {
    smap[c(ctrl,trt)] %>% stack %>%
        rename(Treatment = ind) %>%
        column_to_rownames("values")
}

## Generate Metadata matrices for all dichotomies
M <- map(dich, ~makeMeta(sid, .x[1], .x[2]))

## Extract the corresponding slices of expression data
sliceX <- function(v) {
    select(X, symbol, one_of(v)) %>%
        as.data.frame %>% column_to_rownames("symbol")
}
Xs <- map(M, ~sliceX(rownames(.x)))

## Wrapper for edgeR
edger <- function(X1, Y1) {
    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = X1, samples = Y1 )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~Treatment, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential expression and
    ##   identify Metformin's drug-associated gene list (DGL)
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    edgeR::topTags( gf, nrow(X1) ) %>% as.data.frame %>%
        rownames_to_column( "Gene" ) %>% as_tibble
}

## Apply differential gene expression to each dichotomy
DFX <- map2( Xs, M, edger )
saveRDS( DFX, file="dfx.rds" )

## Compose weight vectors for GSEA
W <- map( DFX, with, set_names(logFC, Gene) )

## Load pathways
c2cp <- DRIAD::read_gmt( "data/c2.cp.v7.0.symbols.gmt" )

## Run GSEA
GSEA <- map(W, ~fgsea::fgsea(c2cp, .x, 1e5)) %>%
    map( as_tibble )
saveRDS( GSEA, file="gsea.rds" )

## Generate dfx sheets for an xlsx file
snm <- map_chr( dich, str_flatten, "-" )
S1 <- map( DFX, filter, FDR < 0.05 ) %>%
    set_names( str_c(snm, "-dfx") )

## Generate gsea sheets
S2 <- map(GSEA, select, pathway, pval, FDR=padj, NES, size) %>%
    map( filter, FDR < 0.05 ) %>% map( arrange, NES ) %>%
    set_names( str_c(snm, "-gsea") )

## Compose a summary sheet
SM <- enframe(sid, "Name", "Sample") %>% unnest(Sample) %>%
    inner_join(Y, by="Sample")

c( list(Summary=SM), S1, S2 ) %>%
    openxlsx::write.xlsx( file="metformin.xlsx" )
