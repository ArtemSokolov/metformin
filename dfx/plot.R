library( tidyverse )

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
    arrange( logFCD ) %>%
    mutate( GeneID = factor(GeneID, GeneID) ) %>%
    select( Treatment, Gene, GeneID, logFC10, logFC40 ) %>%
    gather( Comparison, logFCraw, logFC10, logFC40 ) %>%
    mutate_at( "Comparison", recode,
              logFC10 = "10uM vs DMSO",
              logFC40 = "40uM vs DMSO" ) %>%
    mutate_at( "Treatment", recode,
              Gly24 = "Glyburide 24h",
              Gly72 = "Glyburide 72h",
              Met24 = "Metformin 24h",
              Met72 = "Metformin 72h" ) %>%
    mutate( logFC = case_when(logFCraw > lmt ~ lmt,
                              logFCraw < -lmt ~ -lmt,
                              TRUE ~ logFCraw),
           lbl = str_c(round(logFCraw, 2)) )

flbl <- function( v ) {with(RX, set_names(Gene, GeneID))[v]}

## Plot as a faceted heatmap
gg <- ggplot( RX, aes(x=Comparison, y=GeneID, fill=logFC) ) +
    theme_minimal() + geom_tile() +
    geom_text( aes(label=lbl), data=filter(RX, abs(logFCraw) > 2), color="gold" ) +
    facet_wrap( ~Treatment, ncol=4, scales="free_y" ) +
    scale_y_discrete( labels = flbl, name="Gene" ) +
    scale_fill_gradientn( colors=pal, name="logFC", limits=c(-lmt, lmt) ) +
    theme(axis.text.x = element_text( angle=90, vjust=0.5, hjust=1 ),
          strip.text  = element_text(size=12),
          panel.grid  = element_blank() )

ggsave( "FigA.pdf", gg, width=9, height=7 )

## Save all differential gene expression as a supplementary table
openxlsx::write.xlsx( file="Suppl-Table-A.xlsx" )
