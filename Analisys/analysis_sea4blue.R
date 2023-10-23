

#load packages
library(MiscMetabar)
library(formattable)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(metagenomeSeq)
library(ggrepel)
library(ranacapa)
library(ggside)
library(ggpubr)
library(grid)

#import data for the phyloseq object ( ASVs table, Tax table, metadata table)
ps <- read_pq(path = "sea4blue_phyloseq") #this function import all the files from a specified folder 

ps

# the metadata 
sample_data(ps)


##################################################
####### LET'S VISUAALIZE RAREFACTION CURVES ###### 
##################################################

set.seed(42)


# define a new ps object called psdata that will be used in the following function
ps


#### define  the function calculate_rarefaction_curves

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'none', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}


#### use the function calculate_rarefaction_curves on ps
rarefaction_curve_data <- calculate_rarefaction_curves(ps, c('Observed'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
#print the summary
summary(rarefaction_curve_data)


# do some calculations  and merge the data for plotting 
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, 
                                                data.frame(sample_data(ps)), 
                                                by.x = 'Sample', by.y = 'row.names')



### plot with ggplot

ggplot(data=rarefaction_curve_data_summary_verbose, aes(x = Depth, 
                                                        y = Alpha_diversity_mean, 
                                                        group = Sample, color = Zone)) +
  geom_line(aes(group = Sample), size = 1) +
  geom_point(size = 2.5) +
  geom_label_repel(data = rarefaction_curve_data_summary_verbose %>% 
      group_by(Sample) %>% 
      filter(Alpha_diversity_mean == max(Alpha_diversity_mean)), 
    aes(Depth, Alpha_diversity_mean,  label = ( Sample)),
    color = "black"  )+ggtitle("Rarefaction Curves")



## there is something strange with sample 16 and 17 could you guess why


ggplot(data=rarefaction_curve_data_summary_verbose, aes(x = Depth, 
                                                        y = Alpha_diversity_mean, 
                                                        group = Sample, color = Zone)) +
  geom_line(aes(group = Sample), size = 1) +
  geom_point(size = 2.5) +
  geom_label_repel(data = rarefaction_curve_data_summary_verbose %>% 
                     group_by(Sample) %>% 
                     filter(Alpha_diversity_mean == max(Alpha_diversity_mean)), 
                   aes(Depth, Alpha_diversity_mean,  label = ( SampleName )), #changing the parameter of the labels
                   color = "black"  )+ggtitle("Rarefaction Curves")


#### do you see some patterns???????


#######################################
#######       NORMALIZATION     ###### 
######################################


phyloseq_obj<-ps
#search for Mitochondria
grep(pattern = "Mitochondria", tax_table(phyloseq_obj)) 
#search for Chloroplast
grep(pattern = "Chloroplast", tax_table(phyloseq_obj)) 

#remove Chloroplast and Mitochondria from phyloseq_obj
phyloseq_obj <- phyloseq_obj %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Order) ) 

phyloseq_obj


#filtering singletons
doubleton <- genefilter_sample(phyloseq_obj, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj) 
doubleton

## rimuovo i pos e neg

sample_variables(doubleton)


#rimuovo i positivi e negativi
doubleton = subset_samples(doubleton,!( SampleName=="Positive" | SampleName=="Negative"))
#check it
sample_data(doubleton)$SampleName




# transforming 
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)

p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
head(data.CSS)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
phyloseq_obj_css <- phyloseq_obj
otu_table(phyloseq_obj_css) <- otu_table(data.CSS, taxa_are_rows = T)

phyloseq_obj_css

#change sample names to the correct ones
sample_names(phyloseq_obj_css)<-sample_data(phyloseq_obj_css)$SampleName 

sample_names(phyloseq_obj_css)


############ this is the phyloseq object normalized
phyloseq_obj_css




#######################################
#######     ALPHA DIVERSITY      ###### 
#######################################


# round normalized counts for alpha diversity
phyloseq_obj_css_round <- phyloseq_obj_css
otu_table(phyloseq_obj_css_round) <- round(otu_table(phyloseq_obj_css), digits = 0)

# check it out
head(otu_table(phyloseq_obj_css))
head(otu_table(phyloseq_obj_css_round))




#let's focus on richeness (number of species) and let's check its variation along the longitudinal transect

plot_richness(physeq_normalized,x="Longitude",color = "Zone",measures=c("Observed"))
p1<-plot_richness(physeq_normalized,x="Longitude",color = "Zone",measures=c("Observed"))

#change the order of 
newSTorder = c( "ANW","ANC","ANE")
p1$data$Zone<- as.character(p1$data$Zone)
p1$data$Zone <- factor(p1$data$Zone, levels=newSTorder) 



#plot
ggplot(p1$data,aes(Longitude,value))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(colour = Zone),alpha =0.45)+
  scale_color_manual(values=c("#4280fc", "#ffb452","#f7170a","#000000"))+
  geom_ysideboxplot(aes(x=Zone,y=value,colour = Zone), orientation = "x")+
  theme(        ggside.panel.scale.y = .4)+scale_ysidex_discrete()+
  geom_smooth(aes(colour = variable ),linewidth=1.5)+ ylab('Richness') 




#let's add the richness to our metadata, maybe will be useful later...
sample_data(phyloseq_obj_css)$Richness<-p1$data$value


#######################################
#######     BETA DIVERSITY      #######
#######################################




### calculate bray-curtis dissimilarity among our samples

otu.ord <- ordinate(physeq = phyloseq_obj_css, "PCoA", distance = "bray")


#plot ordination

plot_ordination(phyloseq_obj_css, otu.ord,  color="Zone",axes =c(1,2))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ geom_point(size=3)+ 
  geom_text_repel(aes(label=SampleName),max.overlaps = Inf, show.legend = FALSE)+ 
  theme_void()+
  theme_bw() + 
  geom_xsidedensity(aes(y=stat(density),fill=Zone), alpha = 0.5, show.legend = FALSE) +
  geom_ysidedensity(aes(x=stat(density),fill=Zone), alpha = 0.5, show.legend = FALSE) +
  scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
  scale_ysidex_discrete()+
  ggside::theme_ggside_void()  +
  scale_fill_manual(values = c("#4280fc", "#ffb452","#f7170a"))+
  ggtitle("Beta Diversity based on Bray-Curitis")


#why sample 1 is so different from the others of the "same zone"






#let's check the composition of the communities with a barplot


#we have to transform the counts to %
physeq_perc=transform_sample_counts(phyloseq_obj_css,function(x) 100 * x/sum(x))



#Turn all ASVs into genus counts


glom <- tax_glom(physeq_perc, taxrank = 'Genus')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object

data_glom$Zone <- as.character(data_glom$Zone) #convert to character
data_glom$Zone  <- factor(data_glom$Zone , levels =c("ANW","ANC","ANE"))






ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=Genus)) + 
  facet_grid(~Zone, scales = "free")+
  geom_bar(aes(), stat="identity", position="fill",width = 0.7) +
  scale_fill_manual(values = mycolors)+
  guides(fill=guide_legend(ncol=6)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none")


# can you see any difference ???




#let's do a RDA to see how




#calculate the ordination based on the effect of richness and lon
RDA = ordinate(phyloseq_obj_css ~ Richness + Longitude, "RDA")

# plot it
p2 = plot_ordination(phyloseq_obj_css, RDA, color = "Zone", )+ 
  theme_bw() +geom_text_repel(aes(label = SampleName, color = Zone), max.overlaps = Inf, show.legend = FALSE) + geom_point(size = 4)+
  scale_color_manual(values = c("#4280fc", "#ffb452", "#f7170a"))

p2



# Now add the environmental variables as arrows
arrowmat = vegan::scores(RDA, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = RDA1, yend = RDA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.7 * RDA1, y = 1.7 * RDA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.05, "npc"))
p3 = p2 + geom_segment(arrow_map, size = 1, data = arrowdf, color = "black", 
                       arrow = arrowhead) + geom_text(label_map, size = 4, data = arrowdf)
p3













