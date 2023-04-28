# SPATIAL ANALYSIS OF SARS-CoV2 qPCR+ SAMPLES ACROSS THE HUANAN MARKET BASED ON BINOMIAL GLMM FITS USING DATA OF LIU ET AL. NATURE 2023, https://www.nature.com/articles/s41586-023-06043-2 ####
# PLUS COMPOSITIONAL HEATMAP & COMPOSITIONAL PCA BIPLOT OF MITOCHONDRIAL METAGENOMIC READS OF LIU ET AL. 2023 #### 
# T. Wenseleers, 28 April 2023

library(devtools)
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(brms)
# devtools::install_github("m-clark/visibly")
library(visibly)
library(tidybayes)
library(tidyr)
library(lubridate)
library(mgcv)
library(gamm4)
library(nlme)
library(gstat)
library(splines)
library(splines2)
library(afex)
library(car)
library(ICglm)
library(MASS)
library(emmeans)
if (!require(marginaleffects)) devtools::install_github("tomwenseleers/marginaleffects") # fork of marginaleffects with some extra provisions for multinomial models
library(marginaleffects)
library(dplyr)
library(sf)
library(sp)
library(spatstat)
library(rgeos)
library(dsm)
# install_version("gpclib", version = "1.5-6", repos = "http://cran.r-project.org")
library(gpclib)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(export)
library(compositions)
library(ggforce)
library(ggrepel)
library(concaveman)
library(ggthemes)
# install_github("vqv/ggbiplot", force=TRUE)
library(ggbiplot)
library(pheatmap)
source("./scripts/soap_check.R") # from https://github.com/dill/soap_checker, currently not used
source("./scripts/autocrunch.R")
cols = c("#FFFFFF", brewer.pal(5, "YlOrRd")) # heatmap gradient colours


# 1. LOAD UTILITY FUNCTIONS ####

# convert longitudes & latidues to Northings and Eastings 
# required for soap film smoother, as it is isotropic so treats 1 unit change 
# in either dimension as equal, which isn’t true for lat/long
utm_zone = 50 # appropriate for location of Huanan market
targetproj = sprintf("+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs", utm_zone)

toNorthingsEastings = function(longitudes, latitudes, targetproj=targetproj) {
  if (length(longitudes) != length(latitudes)) {
    stop("Error: Longitude and latitude vectors must have the same length.")
  }
  
  data_points = data.frame(longitude = longitudes, latitude = latitudes)
  
  # Remove rows with NA values
  na_rows = is.na(longitudes) | is.na(latitudes)
  data_points_no_na = data_points[!na_rows, ]
  
  # Convert non-NA coordinates to Northings and Eastings
  coords_sf = st_as_sf(data_points_no_na, coords = c("longitude", "latitude"), crs = 4326)
  coords_utm = st_transform(coords_sf, targetproj)
  coords_utm_df = st_coordinates(coords_utm)
  
  # Add NA values back to the coordinates
  coords_utm_df_with_na = data.frame(x = rep(NA, nrow(data_points)), y = rep(NA, nrow(data_points)))
  coords_utm_df_with_na[!na_rows, ] = coords_utm_df
  
  return(coords_utm_df_with_na)
}

# function to return row-wise z score of matrix x, cf pheatmap:::scale_mat & pheatmap:::scale_rows
row_zscore = function (x) {  m = apply(x, 1, mean, na.rm = T)
                             s = apply(x, 1, sd, na.rm = T)
                             row_z_score = (x - m)/s
                             return (row_z_score) }
# function to return column-wise z score of matrix x, cf pheatmap:::scale_mat & pheatmap:::scale_rows
col_zscore = function (x) {  col_z_score = t(row_zscore(t(x)))
                             return (col_z_score) }



# 2. LOAD DATA ####

## load table S1 & S5 from Liu et al. Nature 2023 ####
tableS1 = read.csv("./tables/Liu-etal_2023_TableS1_samples.csv") %>%
  dplyr::rename("Sold_aquatic" = "Aquatic..type.of.vendor.sold.product.",
                "Sold_seafood" = "Seafood..type.of.vendor.sold.product.",
                "Sold_poultry" = "Poultry..type.of.vendor.sold.product.",
                "Sold_livestock" = "Livestock..type.of.vendor.sold.product.",
                "Sold_wildlife" = "Wildlife..type.of.vendor.sold.product.",
                "Sold_vegetable" = "Vegetable..type.of.vendor.sold.product.",
                "Sold_coldchain" = "Cold.chain..type.of.vendor.sold.product.") %>%
  mutate_at(vars(Sold_aquatic, Sold_seafood, Sold_poultry,
                 Sold_livestock, Sold_wildlife, Sold_vegetable, Sold_coldchain), 
            ~ as.factor(replace_na(., "unknown")))
tableS1$Sampling.location = gsub("Wine", "Wing", tableS1$Sampling.location)
tableS1$Sampling.location = gsub("wine", "wing", tableS1$Sampling.location)
tableS1$Sampling.location = gsub("Sewerage", "Sewage", tableS1$Sampling.location)
tableS1$Sampling.location = gsub("Sewers or sewerage wells", "Sewage well", tableS1$Sampling.location)
tableS1$date = as.Date(tableS1$Sampling.date)

tableS5 = read.csv("./tables/Liu-etal_2023_TableS5_pos_8782_28144.csv")
tableS5$group[tableS5$group==""] = NA


## data of SARS-CoV2 qPCR sample positivity with latitude & longitude (version April 1, based on Table S1 but with lats & longs) ####
data = read.csv("./tables/Liu-etal_2023_market_samples_Apr1.csv")
data$date = as.Date(data$Sampling.date, tryFormats="%d/%m/%Y")
data$SARS.CoV.2.qPCR.result = factor(data$SARS.CoV.2.qPCR.result, levels=c("Negative", "Positive"))
data$SARS.CoV.2.qPCR.result.numeric = as.numeric(data$SARS.CoV.2.qPCR.result)-1
# unique(data$Sample.type)
# "Environmental swab" "Sewage"             "Water drain"
data$Sample.type = factor(data$Sample.type, levels=c("Environmental swab", "Water drain", "Sewage"))
data$Sample_location = factor(data$Sample_location, levels=c("West","East","Warehouse"))
unique(data$Sample.type[!is.na(data$Latitude)]) # data for stalls & water drain samples (Environmental swab & Water drain)
unique(data$Sample.type[!is.na(data$Stall_merged_latitude)]) # data only for stalls (Environmental swab)
# latitude & longitude for stalls (Environmental swabs) or water drain samples
data$lat = data$Stall_merged_latitude
data$lat[is.na(data$lat)&data$Sample.type=="Water drain"] = data$Latitude[is.na(data$lat)&data$Sample.type=="Water drain"]
unique(data[!is.na(data$lat),"Sample.type"]) # Environmental swab Water drain 
data$long = data$Stall_merged_longitude
data$long[is.na(data$long)&data$Sample.type=="Water drain"] = data$Longitude[is.na(data$long)&data$Sample.type=="Water drain"]
unique(data[!is.na(data$long),"Sample.type"]) # Environmental swab Water drain 

data$stall = factor(data$Stall_corrected_merged) # stall, to include as random factor in models
# add types of products sold from Table S1 in Liu et al 2023
data = left_join(data, 
                 tableS1[,c("Sample.ID",
                            "Sold_aquatic",
                            "Sold_seafood",
                            "Sold_poultry",
                            "Sold_livestock",
                            "Sold_wildlife",
                            "Sold_vegetable",
                            "Sold_coldchain")],
                 by="Sample.ID", suffix=c("",""))
# add phylogenetic group of virus (A/B) found at stall based on NGS data (table S5)
data$virus_group = factor(replace_na(left_join(data, tableS5, by="Sample.ID")$group, 
                                     "unknown"),
                          levels=c("unknown", "A", "B"))
# NOTE: A & B group are now sometimes defined based on a very few mapping reads (e.g. 3 in some cases)
# not sure how reliable this is

# note: samples taken in other markets & swerage wells in surrounding areas & Animal swabs excluded
# tableS1[which(sapply(tableS1$Sample.ID, function (id) id %in% data$Sample.ID)==FALSE)[1:90],]

# longitudes & latitudes in Northings & Eastings coordinates
# (required for GAM )
xy = toNorthingsEastings(longitudes=data$long, 
                         latitudes=data$lat,
                         targetproj)
data$x = xy$x
data$y = xy$y
data_sf = st_as_sf(data[!is.na(data$x), ], coords=c("x", "y"), crs=targetproj)
data_sf$Sample.type = droplevels(data_sf$Sample.type)

# TODO: add raw nr of reads & centered logits of prop of reads mapping to particular host species from
# Table S8 (BOLD analyses) (cf Table S8 below)

# TODO: add which stalls were associated with human cases/hosps (+from which group) & add estimated date of onset
# (from map below, group HumanCase, based on closest distance or WHO report or data Marion Koopmans,
# cf https://gvn.org/SARS-CoV-2-response-efforts/gvn-forefront-of-virology-covid-19-webinar-series-featuring-dr-marion-koopmans/?)

# pos & neg environmental swab & water drain samples
positives_env = st_as_sf(data_sf[data_sf$Sample.type=="Environmental swab"&data_sf$SARS.CoV.2.qPCR.result=="Positive", ], coords=c("x", "y"), crs=targetproj)
negatives_env = st_as_sf(data_sf[data_sf$Sample.type=="Environmental swab"&data_sf$SARS.CoV.2.qPCR.result=="Negative", ], coords=c("x", "y"), crs=targetproj)
positives_wdrain = st_as_sf(data_sf[data_sf$Sample.type=="Water drain"&data_sf$SARS.CoV.2.qPCR.result=="Positive", ], coords=c("x", "y"), crs=targetproj)
negatives_wdrain = st_as_sf(data_sf[data_sf$Sample.type=="Water drain"&data_sf$SARS.CoV.2.qPCR.result=="Negative", ], coords=c("x", "y"), crs=targetproj)

# plot(data$x, data$y) # in Northings & Eastings coords
# plot(data$Longitude, data$Latitude) # lats & longs



# # dataset with repeated samples of the same stall aggregated

data_aggregated_per_stall = data %>%
  dplyr::filter(!is.na(x)) %>%
  group_by(Sample.type, stall) %>%
  dplyr::summarise(lat = mean(lat),
                   long = mean(long),
                   x = mean(x),
                   y = mean(y),
                   SARS.CoV.2.qPCR.result.numeric = mean(SARS.CoV.2.qPCR.result.numeric),
                   SARS.CoV.2.qPCR.positive = mean(SARS.CoV.2.qPCR.result.numeric)*n(),
                   SARS.CoV.2.qPCR.negative = n() - mean(SARS.CoV.2.qPCR.result.numeric)*n(),
                   nsamples = n() ) %>%
  as.data.frame()
data_aggregated_per_stall$SARS.CoV.2.qPCR.positive = as.integer(data_aggregated_per_stall$SARS.CoV.2.qPCR.positive)
data_aggregated_per_stall$SARS.CoV.2.qPCR.negative = as.integer(data_aggregated_per_stall$SARS.CoV.2.qPCR.negative)
data_aggregated_per_stall$nsamples = as.integer(data_aggregated_per_stall$nsamples)
data_aggregated_per_stall = data_aggregated_per_stall[data_aggregated_per_stall$Sample.type != "Sewage",]
data_aggregated_per_stall$Sample.type = droplevels(data_aggregated_per_stall$Sample.type)
data_aggregated_per_stall$stall=droplevels(data_aggregated_per_stall$stall)
data_aggregated_per_stall$rowid = 1:nrow(data_aggregated_per_stall)


data_aggregated_per_stall_sf = st_as_sf(data_aggregated_per_stall, coords=c("x", "y"), crs=targetproj)





## remaining suppl tables from Liu et al. Nature 2023 ####

tableS2 = read.csv("./tables/Liu-etal_2023_TableS2_pos_samples.csv")
tableS2$date = as.Date(tableS2$Sampling.date)
tableS3 = read.csv("./tables/Liu-etal_2023_TableS3_NGS_samples.csv")
tableS4 = read.csv("./tables/Liu-etal_2023_TableS4_whole_virus_genomes.csv")
tableS4$date = as.Date(tableS4$Collection.date)
tableS6 = read.csv("./tables/Liu-etal_2023_TableS6_kraken.csv")
tableS7 = read.csv("./tables/Liu-etal_2023_TableS7_DEseq2.csv")
tableS8 = read.csv("./tables/Liu-etal_2023_TableS8_mtDNA_BOLD.csv")
rownames(tableS8) = tableS8$Genus
tableS8$Genus = NULL
tableS8 = data.frame(t(tableS8))
# genus with max nr of reads in each sample :
most_common_genus = colnames(tableS8)[max.col(tableS8, ties.method="first")]
sort(unique(most_common_genus))
# [1] "Alcelaphus"   "Amaurornis"   "Anas"         "Anourosorex"  "Aves"         "Bos"          "Canis"        "Capra"       
# [9] "Casuarius"    "Columba"      "Coturnix"     "Felis"        "Gallinula"    "Gallus"       "Homo"         "Lariscus"    
# [17] "Liomys"       "Marmota"      "Microryzomys" "Nephelomys"   "Nyctereutes"  "Oryctolagus"  "Ovis"         "Phasianus"   
# [25] "Rattus"       "Rhizomys"     "Rupicapra"    "Spilopelia"   "Sus"
# TODO: check a few questionable ones, 
# e.g. Alcelaphus (hartebeest), 
# Microryzomys, Liomys & Nephelomys (S American rodents),
# Rupicapra (chamois),
# Casuarius (Casuary)

tableS8$total_reads = rowSums(tableS8)


# TODO: also needs checking:
# Axis (Axis deer),
# Odocoileus (white-tailed deer),
# (Sika deer, Cervus nippon, were reportedly sold in the WHS, cf WHO report, 
# but haven't read anything about Axis or Odocoileus, unless they were sold 
# under the wrong name)
# Hylobates (gibbon)

Sample.ID = rownames(tableS8)
tableS8$Sample.ID = Sample.ID
tableS8$SARS.CoV.2.qPCR.result = factor(tableS1$SARS.CoV.2.qPCR.result[match(tableS8$Sample.ID, tableS1$Sample.ID)]) # Positive / Negative
tableS8$group = as.character(tableS8$SARS.CoV.2.qPCR.result) # more detailed breakdown in A & B group
tableS8$group[rownames(tableS8) %in% tableS5$Sample.ID[which(tableS5$group=="A")]] = "Positive (A)" # "Env_0020" "Env_0033" "Env_0061"
tableS8$group[rownames(tableS8) %in% tableS5$Sample.ID[which(tableS5$group=="B")]] = "Positive (B)" # "Env_0087" "Env_0088" "Env_0126" "Env_0313" "Env_0346" "Env_0354" "Env_0398"
tableS8$group = factor(tableS8$group, levels=c("Negative", "Positive (A)", "Positive (B)", "Positive"),
                                      labels=c("Negative", "Positive (A)", "Positive (B)", "Positive (unknown group)"))
tableS8$most_common_genus = most_common_genus
data$most_common_genus = tableS8$most_common_genus[match(data$Sample.ID, tableS8$Sample.ID)]


## table with species of interest ####
# based on Xiao et al. Sci. Reports 2021 & SARS-CoV2 susceptibility & Liu et al's extended data Fig. 4 
# + some common species + unique genera that per sample came out with maximum reads
# Table needs some further work & checking
# also should indicate 
# - which ones could conceivably have been sold live
# - which species are susceptible to SARS-CoV2
# - which ones can transmit SARS-CoV2
species = read.csv("./tables/species_of_concern.csv")

# tableS8 subsetted to species/genera of concern + humans
tableS8_subset = cbind(Sample.ID=tableS8$Sample.ID, # host read counts together with qPCR result & sample ID
                       SARS.CoV.2.qPCR.result=tableS8$SARS.CoV.2.qPCR.result,
                       group=tableS8$group,
                       most_common_genus=tableS8$most_common_genus,
                       tableS8[, which(sapply(colnames(tableS8), 
                                        function(genus) genus %in% 
                                          c(unique(species$genus), "total_reads")))]) # unique(species$genus[species$possible_concern=="yes"])))])
genera_toremove = (colSums(tableS8_subset[,-c(1,2,3,4)])==0) # no reads in any sample
minreads = 50 # we only consider samples with >=50 host reads
samples_toremove = (tableS8_subset$total_reads<minreads) 
# sum(samples_toremove) : 27
tableS8_subset = tableS8_subset[-which(samples_toremove),
                                -(which(genera_toremove)+4)]

# tableS8_subset$Other = tableS8_subset$total_reads-rowSums(tableS8_subset[,-c(1,2,3,4,ncol(tableS8_subset))])
tableS8_subset$total_reads = NULL
metadata = left_join(data.frame(Sample.ID=rownames(tableS8_subset)), 
                     data[,c("Sample.ID",
                             "Sampling.date",
                             # "Sample_location",
                             "stall",
                             "Sample.type",
                             "lat",
                             "long",
                             "x",
                             "y",
                             "most_common_genus",
                             "Sold_aquatic",
                             "Sold_seafood",
                             "Sold_poultry",
                             "Sold_livestock",
                             "Sold_wildlife",
                             "Sold_vegetable",
                             "Sold_coldchain")],
                     by="Sample.ID", suffix=c("",""))
rownames(tableS8_subset) = paste0(metadata$Sample.ID, "|", metadata$stall, "|", as.character(metadata$Sample.type))



# 3. HEATMAP OF METAGENOMIC DATA OF LIU ET AL 2023 (READS MAPPING TO mtDNA OF PARTICULAR HOS SPECIES) IN RELATION TO SARS-CoV2 qPCR POSITIVITY ####
# mostly exploratory, better analyses to come
# see e.g. tutorial here https://microbiome.github.io/OMA/microbiome-community.html

most_common_reads = tableS8_subset$most_common_genus
annotation_row = data.frame(group = tableS8_subset$group,
                            most_common_reads = most_common_reads,
                            Sold_livestock = metadata$Sold_livestock, # add types of goods sold
                            Sold_wildlife = metadata$Sold_wildlife,
                            Sold_poultry = metadata$Sold_poultry,
                            Sold_coldchain = metadata$Sold_coldchain,
                            Sold_aquatic = metadata$Sold_aquatic #,
                            #Sold_vegetable = metadata$Sold_vegetable
                            ) 
rownames(annotation_row) = rownames(tableS8_subset)
annotation_col = data.frame(category = species$category[match(colnames(tableS8_subset)[-c(1:4)], 
                                                              species$genus)] ) 
rownames(annotation_col) = colnames(tableS8_subset)[-c(1:4)]
library(randomcoloR)
most_common_reads_unique = sort(unique(most_common_reads))
set.seed(1)
most_common_reads_unique_cols = distinctColorPalette(length(most_common_reads_unique))
names(most_common_reads_unique_cols) = most_common_reads_unique
annotation_colors = list(group = c("Negative"="grey",
                                        "Positive (A)"="magenta",
                                        "Positive (B)"="red",
                                        "Positive (unknown group)"="red3"),
                         Sold_livestock = c("no"="blue",
                                            "unknown"="grey",
                                            "yes"="red") ,
                         Sold_wildlife = c("no"="blue",
                                            "unknown"="grey",
                                            "yes"="red"),
                         Sold_poultry = c("no"="blue",
                                            "unknown"="grey",
                                            "yes"="red"),
                         Sold_coldchain = c("no"="blue",
                                            "unknown"="grey",
                                            "yes"="red"),
                         Sold_aquatic = c("no"="blue",
                                            "unknown"="grey",
                                            "yes"="red"),
                         category = c("human"="green",
                                      "livestock"="magenta",
                                      "poultry & ornamental birds"="red3",
                                      "ornamental birds"="red",
                                      "wildlife"="orange",
                                      "local wildlife"="orange3",
                                      "game meat"="purple",
                                      "pets"="cyan"
                                      ),
                         most_common_reads = most_common_reads_unique_cols)

# heatmap calculated using column-wise z scores of Aitchison centered logratio transformed props of reads mapping to particular host genera
# using 1-correlation as distance metric & UPGMA for clustering
reads = tableS8_subset[,-c(1,2,3,4)] # raw nr of reads mapping to particular host genus
# small pseudocount added to counts to avoid zeros (there is more fancy ways to deal with zeros)
pseudocount=1
reads_prop = sweep(reads+pseudocount,1, rowSums(reads+pseudocount), FUN="/")  # prop of reads mapping to particular host genus, 1E-2 added to avoid zeros
reads_clr = as.data.frame(clr(reads_prop)) # Aitchison centered logratio transformed prop of reads mapping to particular host genera, appropriate for compositional data
reads_clr_colzscore = col_zscore(reads_clr) # Aitchison centered logratio transformed data expressed as column z scores, to make sure that all genera get equal weight

maxabsv = max(abs(reads_clr_colzscore))

pheatmap(reads_clr_colzscore,
         scale="none", # column-wise z scores already pre-calculated, same as using reads_clr as input with scale="column"
         breaks = seq(-3, 3, length.out = 100), # or from -maxabsv to maxabsv
         col = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows="correlation", # using correlation-based distance metric for clustering
         clustering_distance_cols="correlation",
         clustering.method="average", # = UPGMA for clustering
         fontsize_row=5,
         annotation_row=annotation_row,
         annotation_col=annotation_col,
         annotation_colors=annotation_colors)
graph2pdf(file='./plots/heatmap.pdf', width=13, height=14)
graph2png(file='./plots/heatmap.png', width=13, height=14)

# TODO : convert tableS8 to long format & check if sample positivity is predictive
# of enrichment by particular host DNA using multinomial nnet::multinom or
# mclogit::mblogit fit (the latter can take into account overdispersion)
# or test if centered logits of reads mapping on particular host species are
# predictive of SARS-CoV2 qPCR positivity using binomial or multinomial GLM
# (maybe using some regularisation, e.g. using L0glm or glmnet)



# 4. COMPOSITIONAL PCA & BIPLOT OF METAGENOMIC DATA (READS MAPPING TO HOST mtDNA OF PARTICULAR SPECIES AS PROVIDED BY LIU ET AL 2023) IN RELATION TO SARS-COV2 qPCR POSITIVITY ####
# to check association between infection status and presence of particular host species reads
# very crude & mostly exploratory - better analyses to come
pca1 = prcomp(reads_clr, scale=FALSE)

# biplot
scores = as.data.frame(pca1$x)
loadings = as.data.frame(pca1$rotation)
scores_df = data.frame(PC1 = scores[, 1], PC2 = scores[, 2])
scores_df$group = tableS8_subset$group
loadings_df = data.frame(Variable = row.names(loadings), 
                         PC1 = 30*loadings[, 1], PC2 = 30*loadings[, 2])
ggplot(scores_df, aes(x = PC1, y = PC2, colour = group, fill = NA)) +
  geom_point() +
  geom_mark_hull(aes(colour=group, fill=group),
                 expand = unit(0, "mm"),
                 radius = unit(0, "mm"),
                 concavity=100) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.75, color = "blue") +
  geom_text(data = loadings_df, aes(x = PC1, y = PC2, label = Variable), size = 3, color = "blue") +
  theme_few() +
  scale_colour_manual("", values = c("blue", "magenta", "red", "red3")) +
  scale_fill_manual("", values = c("blue", "magenta", "red", "red3")) +
  labs(title = "", x = "PC1", y = "PC2") +
  theme(legend.position = c(0.25,0.17),
        legend.background = element_rect(fill="transparent", colour="transparent")) 

graph2pdf(file='./plots/biplot_compositional PCA.pdf', width=5, height=5)
graph2png(file='./plots/biplot_compositional PCA.png', width=5, height=5)
# some degree of separation between negative & pos B + pos A + pos (unknown group)
# SARS-CoV2 pos associated with humans+livestock, including Axis (Axis deer)+Odocoileus (white-tailed deer), 
# which are maybe of interest?
# based on Table S5, sample Env_0061 (from stall West|7|20-22-24, which was selling livestock & poultry),
# might tentatively contain A group SARS-CoV2 (though this is based on 3 mapping reads only)
# and it also has some reads supposedly mapping onto Axis & Odocoileus (Table S8) 
# (though pig & sheep dominant) (but check if these are not false positives & reliable)
# (other 2 stalls with possible group A SARS-CoV2 are also all from street 7 in West Wing,
# including Env_0020 from stall West|7|15-17 (selling seafood), linked to a human case, and
# Env_0033 (garbage cart)


# alt biplot using ggbiplot
# ggbiplot::ggbiplot(pca1, choices=c(1,2), groups=group,
#                    ellipse=F, circle=T, varname.size = 4) + 
#   geom_mark_hull(aes(colour=group, fill=group),
#                  expand = unit(0, "mm"),
#                  radius = unit(0, "mm"),
#                  concavity=11) +
#   scale_fill_manual("", values=c("blue","red")) +
#   scale_colour_manual(guide=F, values=c("blue","red")) +
#   coord_cartesian(xlim=c(-4,4), ylim=c(-4,4)) +
#   theme_few() 







# 5. LOAD MARKET MAP ####

# load geojson with map of market
map = st_read("./maps/geojson/huanan-market-internal.geojson")
map = st_transform(map, targetproj) # map in Northings and Eastings coords
unique(map$group) 
# "Env-Pos"        "StreetNumber"   "WildlifeVendor" "UnknownMeat"    "MarketMap"      "Env-Neg"        "HumanCase"  "Boundary"

market_map = map[map$group=='MarketMap',]
market_boundary = map[map$group=='Boundary',]
human_cases = map[map$group=='HumanCase',]
human_cases_df = as.data.frame(as(human_cases, 'Spatial')) %>% # as dataframe
  dplyr::rename(x = coords.x1, y = coords.x2)
# Before / After = before or after Dec 20 2019, but WHO report has more detailed date breakdown
# TODO: add which stalls were associated with human cases/hosps & add estimated date of onset
# (based on closest distance, cf https://stackoverflow.com/questions/34242154/merging-two-data-frames-both-with-coordinates-based-on-the-closest-location?)


ggplot() +
  geom_sf(data=market_boundary, 
          color="grey20", fill="white") +
  # geom_sf(data=preds_sf, aes(color = estimate, fill = estimate, pch=I(15), size=I(2))) +
  geom_sf(data=market_map, 
          fill=NA, color="black", lwd=0.01)  +
  geom_sf(data=negatives_env, color="blue") + # positive & negative Environmental swab samples
  geom_sf(data=positives_env, color="red") +
  geom_sf(data=negatives_wdrain, color="blue", pch=17) + # positive & negative Water drain samples
  geom_sf(data=positives_wdrain, color="red", pch=17) +
  # geom_sf(data=human_cases, color="magenta", pch=1, size=I(2)) +
  geom_sf_text(data=map[map$group=='StreetNumber',], aes(label=title, size=I(1.5)),
               color="blue") +
  scale_colour_gradientn("", colors=cols, na.value = "grey90") +
  scale_fill_gradientn("", colors=cols, na.value = "grey90") +
  labs(title="SARS-CoV2 qPCR positivity (%)",
       subtitle="binomial GAM soap film smooth", 
       fill="SARS-CoV2+ (%)", 
       colour="SARS-CoV2+ (%)") +
  theme_bw() +
  theme(legend.position = c(0.08, 0.2),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0))) +
  coord_sf(crs = "+proj=omerc +lonc=90 +lat_0=40 +gamma=28 +alpha=0") + # oblique Mercator projection to get view with market aligned vertically, https://stackoverflow.com/questions/66889558/rotate-ggplot2-maps-to-arbitrary-angles
  xlab("") + ylab("")

# TODO: check if this is consistent with Figs in Liu et al. Nature 2023
# (I believe there are some inconsistencies, so check geocoding & I think some samples may be missing
# if I compare to the figure of Liu et al. - maybe geocoding is not complete in input file)



# positives & negatives as data frame
positives_env_df = as.data.frame(as(positives_env, 'Spatial')) %>%
  dplyr::rename(x = coords.x1, y = coords.x2)

negatives_env_df = as.data.frame(as(negatives_env, 'Spatial')) %>%
  dplyr::rename(x = coords.x1, y = coords.x2)



# Convert market boundaries to GAM-compatible format to use as boundary
# in GAM soap film smooth (in Northings and Eastings coords)
boundary_gam = suppressWarnings(as(as(market_boundary, "Spatial"), "gpc.poly"))
boundary_gam = lapply(boundary_gam, function (slot) { out = slot@pts[[1]]
                                                       out$hole = NULL
                                                       names(out) = c("x", "y")
                                                       return(out)
                                                     } )
# [[1]] is west side of the market, [[2]] is east side


# Make prediction grid for model predictions (falling within market boundary)
nbins = 200
predgrid = rbind(make.soapgrid(boundary_gam[[1]], nbins),
                 make.soapgrid(boundary_gam[[2]], nbins))
# prediction grid points too far from any data point :
# distance threshold dist to decide that prediction grid point is too far from any 
# datapoint to be considered relative to unit circle, so dist=0.05 is 1/20th of the market
dist = 0.05 
toofar_envswabs = exclude.too.far(predgrid$x, predgrid$y, 
                                  data$x[!is.na(data$x)&data$Sample.type=="Environmental swab"], 
                                  data$y[!is.na(data$y)&data$Sample.type=="Environmental swab"], 
                                  dist=dist) 
predgrid_envswabs = predgrid[-which(toofar_envswabs),] # subsetted prediction grid for environmental swab test pos predictions
dev.off()
plot(predgrid_envswabs$x, predgrid_envswabs$y, pch=15, cex=1, xlab="", ylab="") # prediction grid points that will be used (removed points too far from datapoints)

toofar_wdrain = exclude.too.far(predgrid$x, predgrid$y, 
                                  data$x[!is.na(data$x)&data$Sample.type=="Water drain"], 
                                  data$y[!is.na(data$y)&data$Sample.type=="Water drain"], 
                                  dist=dist) 
predgrid_wdrain = predgrid[-which(toofar_wdrain),] # subsetted prediction grid for water drain test pos predictions
dev.off()
plot(predgrid_wdrain$x, predgrid_wdrain$y, pch=15, cex=1, xlab="", ylab="") # prediction grid points that will be used (removed points too far from datapoints)


# alt way of doing this:
# predgrid = expand.grid(x = seq(min(data$x, na.rm=T)-diff(range(data$x, na.rm=TRUE))/3, max(data$x, na.rm=T)+diff(range(data$x, na.rm=TRUE))/3, length.out=nbins),
#                        y = seq(min(data$y, na.rm=T)-diff(range(data$y, na.rm=TRUE))/3, max(data$y, na.rm=T)+diff(range(data$y, na.rm=TRUE))/3, length.out=nbins))
# # only keep points falling within market boundary
# x = predgrid$x
# y = predgrid$y
# inside = inSide(boundary_gam, x, y)
# predgrid = predgrid[inside,]
# nrow(predgrid) # 10938
# plot(predgrid$x, predgrid$y, pch=15, cex=1)

# alt way of doing this:
# predgrid_sf = 
# predgrid_sf = st_intersection(predgrid_sf, market_boundary_sf)
# plot(predgrid_sf)
# predgrid = as.data.frame(st_coordinates(predgrid_sf)) %>%
#  rename(x = X, y = Y)



# 6. 2D BAYESIAN BRM/STAN BINOMIAL GLMM TENSOR SPLINE FIT ####

set_sum_contrasts()
set.seed(1)
fit_brm = brm(SARS.CoV.2.qPCR.positive | trials(nsamples) ~
                (1|stall) + # stall included as a random intercept to take into account sampling dependencies
                t2(x, y, by=Sample.type), # 2D tensor spline with interaction by Sample.type (better waic than model with additive effect for Sample.type)
                family=binomial(logit),
                data=data_aggregated_per_stall)
waic(fit_brm) # 170, SE 23
summary(fit_brm)
summary(fit_brm)$fixed
prior_summary(fit_brm)

# plot of exp(random effect coefficients for stall) as measure of enrichment in
# qPCR SARS-CoV2 positivity at particular stalls relative to background
stall_effects = as.data.frame(ranef(fit_brm, robust=FALSE, probs = c(0.05, 0.95))$stall) %>%
  mutate(stall = rownames(.)) %>%
  arrange(desc(Estimate.Intercept))
stall_effects$stall = factor(stall_effects$stall, 
                             levels=stall_effects$stall[order(stall_effects$Estimate.Intercept, decreasing=F)])
ggplot(data=stall_effects[stall_effects$Estimate.Intercept>0,], 
       aes(y = stall, x = exp(Estimate.Intercept))) + 
  geom_pointrange(aes(xmax = exp(Q95.Intercept), xmin = exp(Q5.Intercept)), color = "darkblue") +
  scale_y_discrete("") +
  geom_vline(xintercept = 1, color = "red") + 
  theme_bw() + 
  theme(text = element_text(size=10)) +
  ylab(NULL) + xlab("Fold enrichment in qPCR positivity\n(mean plus 5% & 95% percentiles)") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512), 
                     labels = c("1", "2", "4", "8", "16", "32", "64", "128", "256", "512")) +
  labs(title="Fold enrichment in qPCR SARS-CoV2 positivity",
       subtitle="based on stall random intercept coefficients of\nbinomial GLMM tensor spline fit")
graph2pdf(file='./plots/stall enrichment positive samples.pdf', width=5, height=5)
graph2png(file='./plots/stall enrichment positive samples.png', width=5, height=5)

# plot(fit_brm)
# conditional_effects(fit_brm)
# plot(conditional_smooths(fit_brm), stype="raster") # logits, or stype="contour"
# pp_check(fit_brm, bdraws = 500)

# fitted values of SARS-Cov2 qPCR positivity per stall 
# (taking into account fixed effects + random intercept for stall) :
preds_brm = predictions(fit_brm, conf_level=0.95, type="link", newdata=data_aggregated_per_stall)
preds_brm_avg = preds_brm %>%  # avg estimates per sample type
                    group_by(Sample.type) %>% 
                    dplyr::summarise(estimate_avg = mean(predicted)) 
preds_brm = left_join(preds_brm, preds_brm_avg)
preds_brm = preds_brm  %>% 
  transform(estimate = plogis(predicted),
            conf.low = plogis(conf.low),
            conf.high = plogis(conf.high),
            estimate_avg = plogis(estimate_avg)) %>%
  arrange(desc(estimate)) %>%
  as.data.frame()
preds_brm$stall = factor(preds_brm$stall, 
                         levels=preds_brm$stall[order(preds_brm$estimate, decreasing=F)][!duplicated(preds_brm$stall[order(preds_brm$estimate, decreasing=F)])])
range(preds_brm$estimate)
data_aggregated_per_stall_with_preds_brm = as.data.frame(left_join(data_aggregated_per_stall, preds_brm))
data_aggregated_per_stall_with_preds_brm_sf = st_as_sf(data_aggregated_per_stall_with_preds_brm, coords=c("x","y"), crs=targetproj)

# plot of fitted qPCR SARS-CoV2 positivity per stall
ggplot(data=preds_brm[preds_brm$estimate >= preds_brm$estimate_avg,], 
       aes(y = stall, x = estimate)) + 
  facet_wrap(~ Sample.type) +
  geom_pointrange(aes(xmax = conf.high, xmin = conf.low), color = "darkblue") +
  scale_y_discrete("") +
  geom_vline(aes(xintercept = estimate_avg), color = "red") + 
  theme_bw() + 
  theme(text = element_text(size=10)) +
  ylab(NULL) + xlab("Modelled qPCR SARS-CoV2 positivity (%)\n(mean plus 95% CIs)") +
  labs(title="Modelled qPCR SARS-CoV2 positivity (%)",
       subtitle="based on binomial GLMM tensor spline fit")
graph2pdf(file='./plots/stall qPCR positivity.pdf', width=5, height=8)
graph2png(file='./plots/stall qPCR positivity.png', width=5, height=8)


# model predictions (global grand means, taking into account fixed effects, but leaving out stall random effect)
# cf https://vincentarelbundock.github.io/marginaleffects/articles/brms.html
# and https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/#different-kinds-of-average-predictions-with-multilevel-models
preds_brm_global = predictions(fit_brm, newdata=
                                 rbind(data.frame(predgrid_envswabs, Sample.type="Environmental swab", nsamples=1), 
                                       data.frame(predgrid_envswabs, Sample.type="Water drain", nsamples=1)),
                               re_formula = NA, # we leave out the random effect per stall part
                               conf_level=0.95, type="link") %>% 
  transform(estimate = plogis(predicted),
            conf.low = plogis(conf.low),
            conf.high = plogis(conf.high)) %>%
  as.data.frame()
range(preds_brm_global$estimate)

preds_brm_global_sf = st_as_sf(preds_brm_global, coords=c("x","y"), crs=targetproj)


# ggplot with model predictions for Environmental swabs
ggplot() +
  geom_sf(data=preds_brm_global_sf[preds_brm_global_sf$Sample.type=="Environmental swab",], # global grand mean fixed effect predictions
          aes(color=estimate*100, fill=estimate*100, pch=I(21), size=I(2))) +
  geom_sf(data=market_map, 
          fill=NA, color="black", lwd=0.01)  +
  geom_sf(data=data_aggregated_per_stall_with_preds_brm_sf[data_aggregated_per_stall_with_preds_brm_sf$Sample.type=="Environmental swab",], 
          aes(color=estimate*100, fill=estimate*100, pch=I(21), size=I(2))) +
  geom_sf(data=data_aggregated_per_stall_with_preds_brm_sf[data_aggregated_per_stall_with_preds_brm_sf$Sample.type=="Environmental swab",], 
          aes(pch=I(21), size=I(2)), color=I("black"), fill=alpha("white", 0)) +
  # geom_sf(data=negatives_env, color="blue") +
  # geom_sf(data=positives_env, color="red") +
  # geom_sf(data=negatives_wdrain, color="blue", pch=17) +
  # geom_sf(data=positives_wdrain, color="red", pch=17) +
  # geom_sf(data=human_cases, color="magenta", pch=1, size=I(2)) +
  geom_sf_text(data=map[map$group=='StreetNumber',], aes(label=title, size=I(1.5)),
               color="blue") +
  scale_colour_gradientn("", colors=cols, na.value = "grey90") +
  scale_fill_gradientn("", colors=cols, na.value = "grey90") +
  geom_sf(data=market_boundary, 
          color="grey20", fill="transparent") +
  labs(title="SARS-CoV2 qPCR positivity (%)",
       subtitle="Environmental swabs, binomial GLMM tensor spline fit\nbrm(SARS.CoV.2.qPCR.positive | trials(nsamples) ~ (1|stall) +\nt2(x, y, by=Sample.type), family=binomial(logit))", 
       fill="SARS-CoV2+ (%)", 
       colour="SARS-CoV2+ (%)") +
  theme_bw() +
  theme(legend.position = c(0.08, 0.2),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0))) +
  coord_sf(crs = "+proj=omerc +lonc=90 +lat_0=40 +gamma=28 +alpha=0") + # oblique Mercator projection to get view with market aligned vertically, https://stackoverflow.com/questions/66889558/rotate-ggplot2-maps-to-arbitrary-angles
  xlab("") + ylab("")
graph2pdf(file="./plots/SARSCoV2_positivity_market_envswabs_binomial_GLMM.pdf", width=6, height=6)
graph2png(file="./plots/SARSCoV2_positivity_market_envswabs_binomial_GLMM.png", width=6, height=6)


# ggplot with model predictions for Water drain samples
ggplot() +
  geom_sf(data=preds_brm_global_sf[preds_brm_global_sf$Sample.type=="Water drain",], # global grand mean fixed effect predictions
          aes(color=estimate*100, fill=estimate*100, pch=I(21), size=I(2))) +
  geom_sf(data=market_map, 
          fill=NA, color="black", lwd=0.01)  +
  geom_sf(data=data_aggregated_per_stall_with_preds_brm_sf[data_aggregated_per_stall_with_preds_brm_sf$Sample.type=="Water drain",], 
          aes(color=estimate*100, fill=estimate*100, pch=I(21), size=I(2))) +
  geom_sf(data=data_aggregated_per_stall_with_preds_brm_sf[data_aggregated_per_stall_with_preds_brm_sf$Sample.type=="Water drain",], 
          aes(pch=I(21), size=I(2)), color=I("black"), fill=alpha("white", 0)) +
  # geom_sf(data=negatives_env, color="blue") +
  # geom_sf(data=positives_env, color="red") +
  # geom_sf(data=negatives_wdrain, color="blue", pch=17) +
  # geom_sf(data=positives_wdrain, color="red", pch=17) +
  # geom_sf(data=human_cases, color="magenta", pch=1, size=I(2)) +
  geom_sf_text(data=map[map$group=='StreetNumber',], aes(label=title, size=I(1.5)),
               color="blue") +
  scale_colour_gradientn("", colors=cols, na.value = "grey90") +
  scale_fill_gradientn("", colors=cols, na.value = "grey90") +
  geom_sf(data=market_boundary, 
          color="grey20", fill="transparent") +
  labs(title="SARS-CoV2 qPCR positivity (%)",
       subtitle="Water drain samples, binomial GLMM tensor spline fit\nbrm(SARS.CoV.2.qPCR.positive | trials(nsamples) ~ (1|stall) +\nt2(x, y, by=Sample.type), family=binomial(logit))", 
       fill="SARS-CoV2+ (%)", 
       colour="SARS-CoV2+ (%)") +
  theme_bw() +
  theme(legend.position = c(0.08, 0.2),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0))) +
  coord_sf(crs = "+proj=omerc +lonc=90 +lat_0=40 +gamma=28 +alpha=0") + # oblique Mercator projection to get view with market aligned vertically, https://stackoverflow.com/questions/66889558/rotate-ggplot2-maps-to-arbitrary-angles
  xlab("") + ylab("")

graph2pdf(file="./plots/SARSCoV2_positivity_market_water drain_binomial_GLMM.pdf", width=6, height=6)
graph2png(file="./plots/SARSCoV2_positivity_market_water drain_binomial_GLMM.png", width=6, height=6)






# 7. TESTS WITH 2D BINOMIAL GAM SOAP FILM SMOOTHS ####

# BINOMIAL GAM SOAP FILM SMOOTH, TAKING INTO ACCOUNT MARKET BOUNDARY
# see tutorial https://fromthebottomoftheheap.net/2016/03/27/soap-film-smoothers/
# and ?soap
# NOTE: unfortunately, these smooths are not supported by brms for the moment

# NOTE2: these fits give a picture somewhat similar to the earlier published
# density plots, but right now they do not incorporate a stall random intercept
# or spatial autocorrelation in the residuals & hence are less good than the fits
# above, as they do not take into account sampling dependencies arising from
# repeated sampling from the same stall

nknots = 7 # gave best GCV score in binomial GAM soap film smooth below
knots = rbind(make.soapgrid(boundary_gam[[1]], nknots), # knots for soap film smooth, falling within market boundary
              make.soapgrid(boundary_gam[[2]], nknots))
dev.off()
plot(knots$x, knots$y, pch=16) # soap film spline knots (all within market boundary)

# fits on data aggregated per stall, with Sample.type acting additively or 
# including Sample.type interaction effect
fit_gam_soap1 = gam(cbind(SARS.CoV.2.qPCR.positive, SARS.CoV.2.qPCR.negative) ~ s(x, y, bs="so", k=5,  # soap film smooth, I chose k visually (also possible to use GCV score etc)
                                                                                  xt=list(bnd=boundary_gam, nmax=150)) +
                      Sample.type, # + # main effect for Sample.type
                    # s(stall, bs="re"), # TODO: add random intercept for stall to take into account sampling dependencies at stall level 
                    # note: a traditional tensor spline would be ~ ti(x)+ti(y)+te(x, y), but this would not take into account
                    # market boundaries, or in context of a binomial GLM ~ ns(x, df=XX):ns(y, df=XX):Sample.type
                    # predictions of all of these are quite similar
                    knots=knots,
                    data = data_aggregated_per_stall,
                    family = binomial(logit))
summary(fit_gam_soap1)
fit_gam_soap2 = gam(cbind(SARS.CoV.2.qPCR.positive, SARS.CoV.2.qPCR.negative) ~ s(x, y, bs="so", k=5,  # soap film smooth, I chose k visually (also possible to use cross-validation or based on GCV score etc)
                                                                                  xt=list(bnd=boundary_gam, nmax=150), 
                                                                                  by=Sample.type), # + # interaction by Sample.type 
                    # s(stall, bs="re"), # TODO: add random intercept for stall to take into account sampling dependencies at stall level 
                    # note: a traditional tensor spline would be ~ ti(x)+ti(y)+te(x, y), but this would not take into account
                    # market boundaries, or in context of a binomial GLM ~ ns(x, df=XX):ns(y, df=XX):Sample.type
                    # predictions of all of these are quite similar
                    knots=knots,
                    data = data_aggregated_per_stall,
                    family = binomial(logit))
summary(fit_gam_soap2)
GCV(fit_gam_soap1)
GCV(fit_gam_soap2) # model with interaction by Sample.type, fit_gam_soap2, fits better




# NOTE: for the moment I am taking into account spatial autocorrelation, but not 
# spatial autocorrelation in the residuals, to manually take into account 
# those in a GAM see
# https://stats.stackexchange.com/questions/459373/gam-with-binomial-distribution-and-with-spatial-autocorrelation-in-r
# or one could resort to a gamm with extra argument correlation=corExp(form = ~ x + y)
# https://stackoverflow.com/questions/57821185/gamm-with-spatial-auto-correlation-in-r 
# e.g.
# data_withnoise = data
# data_withnoise$x = jitter(data_withnoise$x, factor=1E-6) # add tiny bit of noise to coordinates to avoid zero distances
# data_withnoise$y = jitter(data_withnoise$y, factor=1E-6)
# fit_gamm = gamm(SARS.CoV.2.qPCR.result ~ ti(x)+ti(y)+te(x, y) + Sample.type, 
#                 # not this fit has downside that it would not take into account market boundaries, i.e. it uses a
#                 # regular tensor spline with additive effect for Sample.type
#                 correlation=corExp(form = ~ x + y), # spatially autocorrelated residuals
#                 family = binomial(logit),
#                 data = data_withnoise)
# below I will not consider possible spatial autocorrelation in the residuals
# and just use our GAM soap film smooth fits


# plot model predictions of soap film smooth model fit_gam_soap_env_aggr & fit_gam_soap_wdrain_aggr

# model predictions for environmental swab test positivity
fit = fit_gam_soap2
null_value_envswab = mean(predict(fit, newdata=data.frame(predgrid_envswabs,
                                                          Sample.type="Environmental swab"), 
                                  type="link")) # avg logit(prop positive) across market among Environmental swabs
plogis(null_value_envswab) 
# 0.04871769

preds_envswab_df = predictions(fit, newdata=data.frame(predgrid_envswabs,
                                                       Sample.type="Environmental swab"), 
                               vcov=T, conf_level=0.95, type="link") %>%
  transform(estimate = 100*plogis(predicted),
            conf.low = 100*plogis(conf.low),
            conf.high = 100*plogis(conf.high)
  )

preds_envswab_sf = st_as_sf(preds_envswab_df, coords=c("x","y"), crs=targetproj)



# ggplot with model predictions for Environmental swabs
ggplot() +
  geom_sf(data=market_boundary, 
          color="grey20", fill="white") +
  geom_sf(data=preds_envswab_sf, 
          aes(color = estimate, fill = estimate, pch=I(15), size=I(2))) +
  geom_sf(data=market_map, 
          fill=NA, color="black", lwd=0.01)  +
  geom_sf(data=negatives_env, color="blue") +
  geom_sf(data=positives_env, color="red") +
  # geom_sf(data=negatives_wdrain, color="blue", pch=17) +
  # geom_sf(data=positives_wdrain, color="red", pch=17) +
  # geom_sf(data=human_cases, color="magenta", pch=1, size=I(2)) +
  geom_sf_text(data=map[map$group=='StreetNumber',], aes(label=title, size=I(1.5)),
               color="blue") +
  scale_colour_gradientn("", colors=cols, na.value = "grey90") +
  scale_fill_gradientn("", colors=cols, na.value = "grey90") +
  labs(title="SARS-CoV2 qPCR positivity (%)",
       subtitle="Environmental swabs, binomial GAM soap film smooth", 
       fill="SARS-CoV2+ (%)", 
       colour="SARS-CoV2+ (%)") +
  theme_bw() +
  theme(legend.position = c(0.08, 0.2),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0))) +
  coord_sf(crs = "+proj=omerc +lonc=90 +lat_0=40 +gamma=28 +alpha=0") + # oblique Mercator projection to get view with market aligned vertically, https://stackoverflow.com/questions/66889558/rotate-ggplot2-maps-to-arbitrary-angles
  xlab("") + ylab("")
graph2pdf(file="./plots/SARSCoV2_positivity_market_envswabs_binomial_GAM_soap film smooth_no stall random effect.pdf", width=6, height=6)
graph2png(file="./plots/SARSCoV2_positivity_market_envswabs_binomial_GAM_soap film smooth_no stall random effect.png", width=6, height=6)


# model predictions for water drain test positivity
fit = fit_gam_soap2
null_value_wdrain = mean(predict(fit, newdata=data.frame(predgrid_wdrain,
                                                         Sample.type="Water drain"), 
                                 type="link")) # avg logit(prop positive) across market among Water drain samples
plogis(null_value_wdrain) 
# 0.02475341

preds_wdrain_df = predictions(fit, newdata=data.frame(predgrid_wdrain,
                                                      Sample.type="Water drain"), 
                              vcov=T, conf_level=0.95, type="link") %>%
  transform(estimate = 100*plogis(predicted),
            conf.low = 100*plogis(conf.low),
            conf.high = 100*plogis(conf.high)
  )
preds_wdrain_sf = st_as_sf(preds_wdrain_df, coords=c("x","y"), crs=targetproj)

# ggplot with model predictions for Water drain samples
ggplot() +
  geom_sf(data=market_boundary, 
          color="grey20", fill="white") +
  geom_sf(data=preds_wdrain_sf, aes(color = estimate, fill = estimate, pch=I(15), size=I(2))) +
  geom_sf(data=market_map, 
          fill=NA, color="black", lwd=0.01)  +
  # geom_sf(data=negatives_env, color="blue") +
  # geom_sf(data=positives_env, color="red") +
  geom_sf(data=negatives_wdrain, color="blue", pch=17) +
  geom_sf(data=positives_wdrain, color="red", pch=17) +
  # geom_sf(data=human_cases, color="magenta", pch=1, size=I(2)) +
  geom_sf_text(data=map[map$group=='StreetNumber',], aes(label=title, size=I(1.5)),
               color="blue") +
  scale_colour_gradientn("", colors=cols, na.value = "grey90") +
  scale_fill_gradientn("", colors=cols, na.value = "grey90") +
  labs(title="SARS-CoV2 qPCR positivity (%)",
       subtitle="binomial GAM soap film smooth", 
       fill="SARS-CoV2+ (%)", 
       colour="SARS-CoV2+ (%)") +
  theme_bw() +
  theme(legend.position = c(0.08, 0.2),
        legend.key.size = unit(0.2, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.background = element_rect(colour = alpha("white", 0), fill = alpha("white", 0))) +
  coord_sf(crs = "+proj=omerc +lonc=90 +lat_0=40 +gamma=28 +alpha=0") + # oblique Mercator projection to get view with market aligned vertically, https://stackoverflow.com/questions/66889558/rotate-ggplot2-maps-to-arbitrary-angles
  xlab("") + ylab("")

graph2pdf(file="./plots/SARSCoV2_positivity_market_wdrain_binomial_GAM_soap film smooth_no stall random effect.pdf", width=6, height=6)
graph2png(file="./plots/SARSCoV2_positivity_market_wdrain_binomial_GAM_soap film smooth_no stall random effect.png", width=6, height=6)
