

library(phyloseq)
library(ggplot2)
library(vegan)
library(directlabels)
library(gridExtra)
library(dplyr)

############################ rarefaction curves ########################
############################ rarefaction curves ########################
############################ rarefaction curves ########################

setwd("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\PHYLOSEQ OBJECTS\\Rarefaction curves")

merged_8species<-readRDS( "C:\\Users\\risel\\Dropbox\\Sommer postdoc\\PHYLOSEQ OBJECTS\\Rarefaction curves\\merged_8species_unrarefied")
  
#######################################################


p <- ranacapa::ggrare(merged_8species, step = 500,   se = FALSE)+ xlim(c(0,50000))

p + xlim(c(0,50000))


p1<- p + facet_wrap(~Species, scales="free_y", ncol=4)+
  geom_vline(xintercept=10000)+
  theme_bw()+
  theme(legend.position = "none")


rarefaction_fig<- p1+geom_line(aes(col=Species))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(text=element_text(size=14))+
  ylab("Number of ASVs")+
  xlab("Number of reads")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Rarefaction curves with sequencing depth")


rarefaction_fig




###############################  species accumulation curve ######
###############################  species accumulation curve ######
###############################  species accumulation curve ######
###############################  species accumulation curve ######
###############################  species accumulation curve ######

SAClist<-list()

uniq <- unique(sample_data(merged_8species)$Species)



for (i in 1:length(uniq)){
  
  
  data_1<-subset_samples(merged_8species, Species == uniq[i])
  
  data_1<-prune_taxa(taxa_sums(data_1)>0, data_1)
  data_1_matrix<-data.frame(t(data.frame(otu_table(data_1))))
  data_1_specaccum<-specaccum(data_1_matrix, method="random", permutations = 500)
  
  
  sac_df<- data_1_specaccum$sites
  sac_df<-data.frame(sac_df)
  names(sac_df)[1]<-"Site"
  sac_df$Richness <-  data_1_specaccum$richness
  sac_df$SD <-  data_1_specaccum$sd
  
  
  sac_total_estimated<-specpool(data_1_matrix)
  sac_df$Total <- sac_total_estimated$boot
  
  sac_df$Species <- as.character(sample_data(data_1)$Species[1])
  head(sac_df)
  SAClist[[i]]<-sac_df
  
}


sac_df_all<-do.call(rbind, SAClist)

head(sac_df_all)


######### species accumulation figure
######### species accumulation figure
######### species accumulation figure

### add extra lines


species_totals<-sac_df_all %>% distinct(Total, .keep_all = T)

species_totals[,1]<-400
species_totals[,2]<-species_totals$Total
species_totals[,3]<-NA
species_totals[,4]<-NA

species_totals

sac_df_fig<-rbind(sac_df_all, species_totals)




sub<-ggplot(sac_df_fig, aes(x = Site, y = Richness, group = Species))+geom_line(alpha=0.7, linetype = "dashed")+
  geom_point( aes(col=Species), size = 1)+theme_bw()+
  xlab("Number of individuals sampled")+
  ylab("Number of ASVs")+
  theme(text=element_text(size=14))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(legend.position = "none")+ scale_color_brewer(palette = "Dark2")+
  ylim(c(0,3000))+xlim(c(0,50))



sac_fig<-ggplot(sac_df_fig, aes(x = Site, y = Richness, group = Species))+
  geom_line(alpha=0.7, linetype = "dashed")+
  geom_point( aes(col=Species), size = 2, alpha = 0.5)+theme_bw()+
  xlab("Number of individuals sampled")+
  ylab("Number of ASVs")+
  theme(text=element_text(size=14))+
  theme(legend.position = "none")+ 
  scale_color_brewer(palette = "Dark2")+
  annotation_custom(ggplotGrob(sub), xmin = 150, xmax=400, ymin = 3500, ymax = 6000)+
  geom_dl(aes(label = Species), method = list("last.points", cex = 0.8, vjust = -0.2, hjust = 0.5))+
  ggtitle("ASV accumulaiton curves with sample size")



grid.arrange(rarefaction_fig, sac_fig, ncol=2)
