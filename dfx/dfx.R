library( tidyverse )

## Wrapper for edgeR
edger <- function(X1, Y1, cf=2) {
    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = X1, samples = Y1 )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~Treatment, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential expression and
    ##   identify Metformin's drug-associated gene list (DGL)
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = cf )
    edgeR::topTags( gf, nrow(X1) ) %>% as.data.frame %>%
        rownames_to_column( "Gene" ) %>% as_tibble
}

## Wrapper for GSEA
gsea <- function(w, pw) {
    fgsea::fgsea(pw, w, 1e5) %>% as_tibble() %>%
        arrange( desc(abs(NES)) )
}

## NSE wrapper for column-to-rownames
c2rn <- function(.df, .x) {
    as.data.frame(.df) %>% column_to_rownames( str_c(ensym(.x)) )
}

## Preliminaries
synapser::synLogin()
syn <- synExtra::synDownloader("/data/metformin", ifcollision="overwrite.local")
read_csv2 <- partial( read_csv, col_types=cols() )
read_tsv2 <- partial( read_tsv, col_types=cols() )

## Load pathways
c2cp <- DRIAD::read_gmt( "../data/c2.cp.v7.0.symbols.gmt" )

## Identify protein-coding genes
GM <- read_csv2( "../data/GRCh38.94.genemap.csv" ) %>%
    filter( gene_biotype == "protein_coding" )

## Download counts data and matching metadata
Xraw <- syn("syn21763099") %>% read_tsv2 %>% select(-id) %>%
    filter(symbol %in% GM$gene_name) %>%
    group_by(symbol) %>% summarize_all(sum)
Yraw <- syn("syn21781943") %>% read_csv2 %>%
    select( Sample, Well, Treatment )

## Name cleanup
ncu <- c("Glyb 10" = "Glyburide (10uM)",
         "Glyb 40" = "Glyburide (40uM)",
         "Met 10"  = "Metformin (10uM)",
         "Met 40"  = "Metformin (40uM)",
         "control" = "Control")

## Isolate the slice of interest
Y <- Yraw %>% select(-Well) %>%
    mutate_at( "Treatment", str_replace, "24 hr", "24hr" ) %>%
    filter( grepl("24hr", Treatment) | grepl("72hr", Treatment) ) %>%
    mutate(Time      = str_sub(Treatment, -4),
           Treatment = str_sub(Treatment, 1, -8)) %>%
    mutate_at( "Treatment", recode, !!!ncu )
X <- Xraw %>% select( symbol, one_of(Y$Sample) ) %>% c2rn( symbol )

## Drop genes that have no variance across the wells
X <- X[(pmap_dbl(X,max) - pmap_dbl(X,min)) > 0,]

## Compose the differential expression tasks
ftask <- function( Z, trt, tm ) {
    Z %>% filter(grepl(trt, Treatment) |
                 grepl("Control", Treatment),
                 grepl(tm, Time) ) %>%
        mutate_at( "Treatment", as.factor ) %>%
        select( -Time ) %>% c2rn( Sample ) %>%
        { set_names( list(.), str_c(trt, tm) ) }
}

## Isolate data slices
Ydfx <- c(ftask( Y, "Met", "24" ),
          ftask( Y, "Met", "72" ),
          ftask( Y, "Gly", "24" ),
          ftask( Y, "Gly", "72" ))
Xdfx <- map( Ydfx, ~X[,rownames(.x)] )

## Compute differential gene expression
DFX <- map2( Xdfx, Ydfx, edger, cf=2:3 ) %>%
    map( rename, logFC10=2, logFC40=3 )

save( Xdfx, Ydfx, DFX, file="all-dfx.RData" )

## Gene set enrichment analysis
W <- map( DFX, with, set_names(LR, Gene) )
EGS <- map( W, gsea, c2cp )

map( EGS, select, -nMoreExtreme, -leadingEdge ) %>%
    map( arrange, pval )
save( EGS, file="all-gsea.RData" )

## Filters to dose-dependent perturbations
fdosedep <- function( Z ) {
    Z %>% filter(sign(logFC40) == sign(logFC10),
                 abs(logFC40) > abs(logFC10))
}

DDX <- map( DFX, fdosedep ) %>%
    map( arrange, FDR ) %>%
    map( slice, 1:80 )

map( EGS, select, -leadingEdge, -nMoreExtreme ) %>%
    map( arrange, pval )

## IFN-gamma highlight
ifg <- c2cp[["REACTOME_INTERFERON_GAMMA_SIGNALING"]]
vifg <- DFX$Met24 %>% filter( Gene %in% ifg, FDR < 0.01 ) %>% pull(Gene)
Xdfx$Met24[vifg,]
