######################################
#########  MANHATTAN PLOTS    ########
#####         H Steiner         ######
### steiner@pharmacy.arizona.edu #####
######################################

# change the name of the output file at the bottom of this page 

## load packages
library(readr)
library(ggrepel)
library(tidyverse)
library(data.table)


## load data desktop
dat <- read_table2("HIPA2_AgeSexPCA12_nonimp_whites_12_5_19.assoc.logistic")

## load gene refernce data
gene_result <- read_delim("Genes_grch37.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

gene_result = gene_result %>% 
  mutate(chrom = str_replace(gene_result$chrom, "chr", "")) %>% 
  select(-name)


## define significance 
sig = 5e-8 # significant threshold line
sugg = 1e-5 # suggestive threshold line




don <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP, -desc(P)) %>%
  mutate( BPcum=BP+tot) %>% 
  distinct(SNP, .keep_all = T)


## REMOVE P FILTER BEFORE PUBLISHING
don2 = don %>% filter(CHR <= 22 & CHR > 0,
                      P < .1)

print("CONGRATULATIONS! Data was cleaned")

## for testing only comment out for plots 

## define genes for annotating
## just top snp
topsnp = don2[which.min(don2$P),]

## top SNP per chromosome
topsnps = don2 %>%  group_by(CHR) %>% slice(which.min(P)) %>% filter(P < sugg)

## this gives just TOP Gene
topgene = gene_result %>% filter(topsnp$CHR == gene_result$chrom & 
                                    topsnp$BP >  gene_result$txStart &
                                    gene_result$txEnd > topsnp$BP)


## top genes per chromosome
topgene_out_fn = function(topsnps){
  gene_result %>%     filter(topsnps$CHR == gene_result$chrom &
                               topsnps$BP >  gene_result$txStart - 500000 &                                     gene_result$txEnd + 500000 > topsnps$BP)
  
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

## remove dups? 
is_duplicate_out <- lapply(X = topgenes_out, FUN = duplicated, incomparables = FALSE)


drop_idx_out <- which(is_duplicate_out$name)
topgenes_out = topgenes_out[-drop_idx_out, ]


is_duplicate_in <- lapply(X = topgenes_in, FUN = duplicated, incomparables = FALSE)


drop_idx_in<- which(is_duplicate_in$name)
topgenes_in = topgenes_in[-drop_idx_in, ]





topgenes = list()

for (i in 1:nrow(topsnps)) {
  if(topsnps$CHR[i] %in% topgenes_in$chrom){
    topgenes[[i]] = topgenes_in[topgenes_in$chrom == topsnps$CHR[i],]}
  else{
   targetgenes =  topgenes_out[topgenes_out$chrom == topsnps$CHR[i],]
   
   ##Get the distance to each target gene
   distances<-list()
   for(j in 1:nrow(targetgenes)){
     if(targetgenes$exonStarts[j] > topsnps$BP[i]){
       distances[[j]] <-targetgenes$exonStarts[j]-topsnps$BP[i]
     }else{
       distances[[j]] <- topsnps$BP[i]-targetgenes$exonEnds[j]
     }
   }
   
   topgenes[[i]] =targetgenes[which.min(unlist(distances)),]
  }
  print(i) 
}

topgenes = do.call(rbind,topgenes)


## define snps of interest
snpsOfInterest = don2 %>% filter(topsnp$BP -50000 <  don2$BP & don2$BP < topsnp$BP +50000 & 
                                  don2$CHR == topsnp$CHR)

don2 = don2 %>% 
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(BP %in% snpsOfInterest$BP, "yes", "no")) %>% 
  mutate( is_annotate=ifelse(SNP %in% topsnps$SNP, "yes", "no")) 


### prepare x axis 
## chromosome name 
axisdf = don2 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


## make the plot 
ggplot(don2, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#81D3EB", "#0C234B"), each = 1, len = 22)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 2, 4, 6, 8)) +     # remove space between plot area and x axis
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -log10(sig), color = "gray") +
  geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "gray") +
  
  # Add highlighted points
  geom_point(data=subset(don2, is_highlight=="yes"), color="#EF4056", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don2, is_annotate=="yes"), aes(label= topgenes$name2), size=3,
                  segment.color = "transparent", nudge_x = -2) +

  expand_limits(y = 10) + 
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(x = "", y = "") 



ggsave("manhattan.png", 
plot = last_plot(), 
device = "png", 
width = 12, height = 6, units = "in")
