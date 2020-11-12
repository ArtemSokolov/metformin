library( tidyverse )
library( grid )
library( cowplot )

load( "all-dfx.RData" )

## Filter to top 80 dose-dependent perturbations
DDX <- DFX %>%
    map( filter, sign(logFC40) == sign(logFC10),
        abs(logFC40) > abs(logFC10) ) %>%
    map( arrange, FDR ) %>%
    map( slice, 1:40 )

## Ceiling to the closest 0.1
ceil.1 <- function( v ) {ceiling( v * 10 ) / 10}

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
lmt <- 2

## Combine everything into a single data frame
RX <- bind_rows( DDX, .id = "Treatment" ) %>%
    mutate(logFCD = logFC40 - logFC10,
           GeneID = str_c(Treatment, Gene)) %>%
    arrange( logFCD )

## Per Mark's request, distill the figure to top4 and bottom4
RX4 <- RX %>%
    group_by( Treatment ) %>%
    slice( c(1:4, 37:40) ) %>%
    ungroup() %>%
    mutate( GeneID = factor(GeneID, GeneID) ) %>%
    select( Treatment, Gene, GeneID, logFC10, logFC40 ) %>%
    gather( Comparison, logFCraw, logFC10, logFC40 ) %>%
    mutate_at( "Comparison", recode,
              logFC10 = "10uM",
              logFC40 = "40uM" ) %>%
    mutate_at( "Treatment", recode,
              Gly24 = "  24h",
              Gly72 = "  72h",
              Met24 = "   24h",
              Met72 = "   72h" ) %>%
    mutate(Timepoint = str_sub(Treatment, -3),
           Treatment = str_sub(Treatment, 1, -5)) %>%
    mutate( logFC = case_when(logFCraw > lmt ~ lmt,
                              logFCraw < -lmt ~ -lmt,
                              TRUE ~ logFCraw),
           lbl = str_c(round(logFCraw, 2)) )

flbl <- function( v ) {with(RX, set_names(Gene, GeneID))[v]}

## Plot as a faceted heatmap
gg <- ggplot( RX4, aes(x=Comparison, y=GeneID, fill=logFC) ) +
    theme_minimal() + geom_tile() +
    geom_text( aes(label=lbl), data=filter(RX4, abs(logFCraw) > 2), color="gold" ) +
    facet_wrap( ~Treatment+Timepoint, scales="free_y", ncol=4 ) +
    scale_y_discrete( labels = flbl, name="" ) + xlab("") +
    ggtitle("Genes with greatest dose-dependent increase or decrease relative to DMSO") +
    scale_fill_gradientn( colors=pal, name="Log Fold\nChange", limits=c(-lmt, lmt) ) +
    theme(axis.text = element_text(face="bold"),
          strip.text  = element_text(size=12),
          panel.grid  = element_blank(),
          plot.title = element_text(hjust = 0.5))

ffacet <- function( text )
    {
        ggplot() + theme_void() +
            geom_line( aes(x=x, y=y), data=data.frame(x=c(0,1), y=c(0.75,0.75)) ) +
            geom_label( aes(x=x, y=y, label=txt),
                       data=data.frame(x=0.5, y=0.75, txt=text),
                       label.size = NA, size=5 ) +
            scale_y_continuous( limits=c(0,1) )
    }

ggg <- ggplotGrob(gg) %>%
    gtable::gtable_add_grob( ggplotGrob(ffacet("Glyburide")), t=7, l=5, r=9 ) %>%
    gtable::gtable_add_grob( ggplotGrob(ffacet("Metformin")), t=7, l=13, r=17 )

ggsave( "Fig-dfx4.pdf", ggg, width=9, height=3 )

## Save all differential gene expression as a supplementary table
openxlsx::write.xlsx( DFX, file="Suppl-Table-dfx.xlsx" )

## Save top 10 gene sets enriched in each comparison to a suppl. table
load( "all-gsea.RData" )
EGS %>% map( select, pathway, pval, NES ) %>%
    map( arrange, pval ) %>% map( slice, 1:10 ) %>%
    openxlsx::write.xlsx( file="Suppl-Table-gsea.xlsx" )

## Highlight ifn-g signaling
c2cp <- DRIAD::read_gmt( "../data/c2.cp.v7.0.symbols.gmt" )
ifng <- c2cp$REACTOME_INTERFERON_GAMMA_SIGNALING
v <- DFX$Met24 %>% filter( Gene %in% ifng, FDR < 0.01 ) %>% pull(Gene)
X <- Xdfx$Met24[v,] %>% rownames_to_column("Gene") %>% as_tibble
Y <- Ydfx$Met24 %>% rownames_to_column("Sample") %>% as_tibble

## Normalize relative to the median control
RR <- X %>% mutate( Med = pmap_dbl(list(S31,S32,S33), lift_vd(median)) ) %>%
    gather( Sample, Counts, -Gene, -Med ) %>%
    inner_join( Y, by="Sample" ) %>%
    mutate(logFC  = log(Counts / Med),
           Sample = str_c(Sample, " ", Treatment))

## Hierarchical clustering on rows
library( seriation )
DM <- RR %>% select( Gene, Sample, logFC ) %>% spread( Sample, logFC ) %>%
    as.data.frame %>% column_to_rownames("Gene") %>% dist
lvl <- hclust(DM) %>% reorder(DM) %>% dendextend::order.hclust() %>% labels(DM)[.]
RR <- RR %>% mutate_at( "Gene", factor, levels=lvl )

## Plot log fold change relative to median control
ggifg <- ggplot( RR, aes(x=Sample, y=Gene, fill=logFC) ) +
    theme_minimal() + geom_tile() +
    scale_fill_gradientn( colors=pal, name="logFC", limits=c(-1.8, 1.8) ) +
    theme(axis.text.x = element_text( angle=90, vjust=0.5, hjust=1 ),
          panel.grid = element_blank() ) +
    ggsave( "Suppl-Fig-IFNg.pdf", width=5, height=5 )
