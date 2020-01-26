######################################
#########  MANHATTAN PLOTS    ########
#####         H Steiner         ######
### steiner@pharmacy.arizona.edu #####
######################################


## load packages
library(readr)
library(ggrepel)
library(tidyverse)
library(data.table)

## load data desktop
dat <- read_table2("HIPA2_AgeSexPCA12_imputed1K_whites_8_7_17.assoc.logistic")

## load gene refernce data
gene_result <- read_delim("Genes_grch37.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


## remove letters from chromosome column
gene_result$chrom = str_replace(gene_result$chrom, "chr", "")


## define significance 
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

## define genes for annotating
## just top snp
topsnp = dat[which.min(dat$P),]

## top SNP per chromosome
topsnps = dat %>%  group_by(CHR) %>% slice(which.min(P)) %>% filter(P < sugg)

## this gives just TOP Gene
topgene = gene_result %>% filter(topsnp$CHR == gene_result$chrom & 
                                    topsnp$BP >  gene_result$txStart &
                                    gene_result$txEnd > topsnp$BP)

## top genes per chromosome 
topgene_out_fn = function(topsnps){
  gene_result %>% 
    filter(topsnps$CHR == gene_result$chrom & 
                                    topsnps$BP >  gene_result$txStart - 500000 &
                                    gene_result$txEnd + 500000 > topsnps$BP)
  
}
topgene_in_fn = function(topsnps){
  gene_result %>% 
    filter(topsnps$CHR == gene_result$chrom & 
             topsnps$BP >  gene_result$txStart &
             gene_result$txEnd  > topsnps$BP)
  
}

topgenes_out = by(topsnps, topsnps$CHR, function(topsnps) topgene_out_fn(topsnps))
topgenes_in = by(topsnps, topsnps$CHR, function(topsnps) topgene_in_fn(topsnps))

topgenes_out = do.call(rbind.data.frame, topgenes_out)
topgenes_in = do.call(rbind.data.frame, topgenes_in)


###### HELP ######################################
## this doesn't work when SNP is OUTSIDE  gene  
topgenes_out = topgenes_out %>% group_by(chrom) %>%
  arrange(txStart) %>% 
  slice(which.min(abs(topsnps$BP - txEnd) | abs( txStart - topsnps$BP)))


topgenes_out$chrom = as.numeric(topgenes_out$chrom)
topgenes_in$chrom = as.numeric(topgenes_in$chrom)

topgenes_out = topgenes_out[order(topgenes_out$chrom), ]
topgenes_in = topgenes_in[order(topgenes_in$chrom), ]
topgenes_in = topgenes_in %>% distinct(name2, .keep_all = T)


## define snps of interest
snpsOfInterest = dat %>% filter(topsnp$BP -50000 <  dat$BP & dat$BP < topsnp$BP +50000 & 
                                  dat$CHR == topsnp$CHR)

topsnps_in = topsnps %>% filter(CHR %in% topgenes_in$chrom)


don <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  filter(CHR <= 22 & CHR > 0) %>%
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>% 
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest$SNP, "yes", "no")) %>% 
  mutate( is_annotate=ifelse(SNP %in% topsnps_in$SNP, "yes", "no")) 

print("CONGRATULATIONS! Data was cleaned")

### prepare x axis 
## chromosome name 
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


## make the plot 
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#81D3EB", "#0C234B"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 2, 4, 6, 8)) +     # remove space between plot area and x axis
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -log10(sig), color = "gray") +
  geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "gray") +
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="#EF4056", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label= topgenes_in$name2), size=3,
                    segment.color = "transparent", nudge_x = 4) +

  expand_limits(y = 10) + 
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "", y = "") 


ggsave("manhattan.png", 
plot = plot, 
device = "png", 
width = 12, height = 6, units = "in")
