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
         "Met 10" = "Metformin (10uM)",
         "Met 40" = "Metformin (40uM)")

## Isolate the slice of interest
Y <- Yraw %>% select(-Well) %>%
    filter( grepl("72hr", Treatment) ) %>%
    mutate(Treatment = str_sub(Treatment, 1, -8)) %>%
    mutate_at( "Treatment", recode, !!!ncu )
X <- Xraw %>% select( symbol, one_of(Y$Sample) ) %>%
    as.data.frame %>% column_to_rownames( "symbol" )

## Drop genes that have no variance across the wells
X <- X[(pmap_dbl(X,max) - pmap_dbl(X,min)) > 0,]

## Compose the differential expression tasks
fsplit <- function( Z, name ) {
    list( -(7:9), -(4:6), -(1:3) ) %>%
        set_names( map_chr(c("12", "13", "23"), ~str_c(name,.x)) ) %>%
        map( ~Z[.x,] ) %>% map( mutate_at, "Treatment", as.factor ) %>%
        map( as.data.frame ) %>% map( column_to_rownames, "Sample" )
}

Ydfx <- c(fsplit(filter( Y, !grepl("Gly", Treatment) ), "Met"),
          fsplit(filter( Y, !grepl("Met", Treatment) ), "Gly"))
Xdfx <- map( Ydfx, ~X[,rownames(.x)] )

## Compute differential gene expression
DFX <- map2( Xdfx, Ydfx, edger )

R12 <- DFX$Met12 %>% arrange( Gene )
R23 <- DFX$Met23 %>% arrange( Gene )
R13 <- DFX$Met13 %>% arrange( Gene )

RR <- inner_join(select(R12, -LR, -PValue),
                 select(R13, -LR, -PValue),
                 by="Gene", suffix=c(".12", ".13") )

Z <- RR %>% filter(FDR.12 < 0.05 | FDR.13 < 0.05,
                   sign(logFC.13) == sign(logFC.12),
                   abs(logFC.13) > abs(logFC.12))

w <- with( R13, set_names( LR, Gene ) )
res <- fgsea::fgsea(c2cp, w, 1e5) %>% as_tibble() %>%
    arrange( desc(abs(NES)) )
