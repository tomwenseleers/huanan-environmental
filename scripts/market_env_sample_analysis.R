# SPATIAL ANALYSIS OF SARS-CoV2 qPCR+ SAMPLES ACROSS THE HUANAN MARKET BASED ON BINOMIAL GAM FITS USING DATA OF LIU ET AL. NATURE 2023, https://www.nature.com/articles/s41586-023-06043-2 ####
# T. Wenseleers, 8 April 2023

library(lubridate)
library(mgcv)
library(gamm4)
library(nlme)
library(gstat)
library(splines)
library(splines2)
library(emmeans)
library(marginaleffects)
library(dplyr)
library(sf)
library(sp)
library(spatstat)
library(rgeos)
library(dsm)
library(devtools)
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


# some utility functions ####

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




# LOAD DATA ####

# tables Liu et al. Nature 2023
tableS1 = read.csv("./tables/Liu-etal_2023_TableS1_samples.csv")
tableS1$date = as.Date(tableS1$Sampling.date)
tableS2 = read.csv("./tables/Liu-etal_2023_TableS2_pos_samples.csv")
tableS2$date = as.Date(tableS2$Sampling.date)
tableS3 = read.csv("./tables/Liu-etal_2023_TableS3_NGS_samples.csv")
tableS4 = read.csv("./tables/Liu-etal_2023_TableS4_whole_virus_genomes.csv")
tableS4$date = as.Date(tableS4$Collection.date)
tableS5 = read.csv("./tables/Liu-etal_2023_TableS5_pos_8782_28144.csv")
tableS5$group[tableS5$group==""] = NA
tableS6 = read.csv("./tables/Liu-etal_2023_TableS6_kraken.csv")
tableS7 = read.csv("./tables/Liu-etal_2023_TableS7_DEseq2.csv")
tableS8 = read.csv("./tables/Liu-etal_2023_TableS8_mtDNA_BOLD.csv")


# data of SARS-CoV2 qPCR sample positivity (version April 1, based on Table S1 with lats & longs)
data = read.csv("./tables/Liu-etal_2023_market_samples_Apr1.csv")
data$date = as.Date(data$Sampling.date, tryFormats="%d/%m/%Y")
data$SARS.CoV.2.qPCR.result = factor(data$SARS.CoV.2.qPCR.result, levels=c("Negative", "Positive"))
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

xy = toNorthingsEastings(longitudes=data$long, # longitudes & latitudes in Northings & Eastings coordinates
                         latitudes=data$lat,
                         targetproj)
data$x = xy$x
data$y = xy$y
data_sf = st_as_sf(data[!is.na(data$x), ], coords=c("x", "y"), crs=targetproj)
data_sf$Sample.type = droplevels(data_sf$Sample.type)

# TO DO: add which stalls were associated with human cases/hosps & add estimated date of onset
# (from map below, group HumanCase, based on closest distance?)


# pos & neg environmental swab & water drain samples
positives_env = st_as_sf(data_sf[data_sf$Sample.type=="Environmental swab"&data_sf$SARS.CoV.2.qPCR.result=="Positive", ], coords=c("x", "y"), crs=targetproj)
negatives_env = st_as_sf(data_sf[data_sf$Sample.type=="Environmental swab"&data_sf$SARS.CoV.2.qPCR.result=="Negative", ], coords=c("x", "y"), crs=targetproj)
positives_wdrain = st_as_sf(data_sf[data_sf$Sample.type=="Water drain"&data_sf$SARS.CoV.2.qPCR.result=="Positive", ], coords=c("x", "y"), crs=targetproj)
negatives_wdrain = st_as_sf(data_sf[data_sf$Sample.type=="Water drain"&data_sf$SARS.CoV.2.qPCR.result=="Negative", ], coords=c("x", "y"), crs=targetproj)

# plot(data$x, data$y) # in Northings & Eastings coords
# plot(data$Longitude, data$Latitude) # lats & longs


# LOAD MARKET MAP ####

# load geojson with map of market
map = st_read("./maps/geojson/huanan-market-internal.geojson")
map = st_transform(map, targetproj) # map in Northings and Eastings coords
unique(map$group) 
# "Env-Pos"        "StreetNumber"   "WildlifeVendor" "UnknownMeat"    "MarketMap"      "Env-Neg"        "HumanCase"  "Boundary"

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

# TO DO: check if this is consistent with Figs in Liu et al. Nature 2023
# (I believe there are some inconsistencies, so check geocoding)


market_map = map[map$group=='MarketMap',]
market_boundary = map[map$group=='Boundary',]
human_cases = map[map$group=='HumanCase',]
human_cases_df = as.data.frame(as(human_cases, 'Spatial')) %>% # as dataframe
  dplyr::rename(x = coords.x1, y = coords.x2)
# Before / After = before or after Dec 20 2019, but WHO report has more detailed date breakdown
# TO DO: add which stalls were associated with human cases/hosps & add estimated date of onset
# (based on closest distance, cf https://stackoverflow.com/questions/34242154/merging-two-data-frames-both-with-coordinates-based-on-the-closest-location?)

# positives & negatives as data frame
positives_env_df = as.data.frame(as(positives_env, 'Spatial')) %>%
  dplyr::rename(x = coords.x1, y = coords.x2)

negatives_env_df = as.data.frame(as(negatives_env, 'Spatial')) %>%
  dplyr::rename(x = coords.x1, y = coords.x2)



# Convert market boundaries to GAM-compatible format to use as boundary
# in soap film smooth (in Northings and Eastings coords)
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



# 2D BINOMIAL GAM SPLINE FITS ####

# BINOMIAL GAM SOAP FILM SMOOTH, TAKING INTO ACCOUNT MARKET BOUNDARY ####
# see tutorial https://fromthebottomoftheheap.net/2016/03/27/soap-film-smoothers/
# and ?soap

nknots = 7 # gave best AIC in binomial GAM soap film smooth below
knots = rbind(make.soapgrid(boundary_gam[[1]], nknots), # knots for soap film smooth, falling within market boundary
              make.soapgrid(boundary_gam[[2]], nknots))
dev.off()
plot(knots$x, knots$y, pch=16) # soap film spline knots (all within market boundary)

fit_gam_soap1 = gam(SARS.CoV.2.qPCR.result ~ s(x, y, bs="so", k=7,  # soap film smooth, I chose k visually (also possible to use cross-validation or based on AIC or BIC etc)
                                              xt=list(bnd=boundary_gam, nmax=150)) +
                                             Sample.type, # + # main effect for Sample.type
                                             # s(stall, bs="re"), # random intercept for stall to take into account sampling dependencies at stall level 
                   # note: a traditional tensor spline would be ~ ti(x)+ti(y)+te(x, y), but this would not take into account
                   # market boundaries, or in context of a binomial GLM ~ ns(x, df=XX):ns(y, df=XX):Sample.type
                   # predictions of all of these are quite similar
                   knots=knots,
                   data = data,
                   family = binomial(logit))
fit_gam_soap2 = gam(SARS.CoV.2.qPCR.result ~ s(x, y, bs="so", k=7,  # soap film smooth, I chose k visually (also possible to use cross-validation or based on AIC or BIC etc)
                                               xt=list(bnd=boundary_gam, nmax=150), 
                                               by=Sample.type), # + # interaction by Sample.type 
                                             # s(stall, bs="re"), # random intercept for stall to take into account sampling dependencies at stall level 
                    # note: a traditional tensor spline would be ~ ti(x)+ti(y)+te(x, y), but this would not take into account
                    # market boundaries, or in context of a binomial GLM ~ ns(x, df=XX):ns(y, df=XX):Sample.type
                    # predictions of all of these are quite similar
                    knots=knots,
                    data = data,
                    family = binomial(logit))
AIC(fit_gam_soap1, fit_gam_soap2) # model with interaction by Sample.type, fit_gam_soap2, fits better

# NOTE: for the moment I am taking into account spatial autocorrelation, but not 
# possible spatial autocorrelation in the residuals, to manually take into account 
# those in a GAM see
# https://stats.stackexchange.com/questions/459373/gam-with-binomial-distribution-and-with-spatial-autocorrelation-in-r
# or one could resort to a gamm with extra argument correlation=corExp(form = ~ x + y)
# https://stackoverflow.com/questions/57821185/gamm-with-spatial-auto-correlation-in-r 
# e.g.
data_withnoise = data
data_withnoise$x = jitter(data_withnoise$x, factor=1E-6) # add tiny bit of noise to coordinates to avoid zero distances
data_withnoise$y = jitter(data_withnoise$y, factor=1E-6)
fit_gamm = gamm(SARS.CoV.2.qPCR.result ~ ti(x)+ti(y)+te(x, y) + Sample.type, 
                # not this fit has downside that it would not take into account market boundaries, i.e. it uses a
                # regular tensor spline with additive effect for Sample.type
                correlation=corExp(form = ~ x + y), # spatially autocorrelated residuals
                family = binomial(logit),
                data = data_withnoise)
# below I will not consider possible spatial autocorrelation in the residuals
# and just use our soap film smooth fit fit_gam_soap2


# plot model predictions of soap film smooth model fit_gam_soap2
fit = fit_gam_soap2
cols = c("#FFFFFF", brewer.pal(5, "YlOrRd")) # heatmap gradient colours

# model predictions for environmental swab test positivity
null_value_envswab = mean(predict(fit, newdata=data.frame(predgrid_envswabs,
                                                  Sample.type="Environmental swab"), 
                          type="link")) # avg logit(prop positive) across market among Environmental swabs
plogis(null_value_envswab) 
# 0.04713202

preds_envswab_df = predictions(fit, newdata=data.frame(predgrid_envswabs,
                                                       Sample.type="Environmental swab"), 
                               vcov=T, conf_level=0.95, type="link") %>%
  transform(estimate = 100*plogis(estimate),
            conf.low = 100*plogis(conf.low),
            conf.high = 100*plogis(conf.high)
  )
# z test to test for sign deviation from uniform constant prop across market
# still need to check this properly
preds_envswab_df$z = (preds_envswab_df$estimate-null_value_envswab)/preds_envswab_df$std.error
pthresh = 1E-10
preds_envswab_df$p.2tailed = p.adjust(2*pnorm(-abs(preds_envswab_df$z)), "bonferroni") # is modelled prop sign different from uniform?
preds_envswab_df$p.2tailed[preds_envswab_df$p.2tailed>pthresh] = NA
preds_envswab_df$p.lefttailed = p.adjust(pnorm(preds_envswab_df$z, lower.tail = T), "bonferroni") # is modelled prop sign lower than constant uniform value?
preds_envswab_df$p.lefttailed[preds_envswab_df$p.lefttailed>pthresh] = NA
preds_envswab_df$p.righttailed = p.adjust(pnorm(preds_envswab_df$z, lower.tail = F), "bonferroni") # is modelled prop sign higher than constant uniform value?
preds_envswab_df$p.righttailed[preds_envswab_df$p.righttailed>pthresh] = NA

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
graph2pdf(file="./plots/SARSCoV2_positivity_market_envswabs_binomial_GAM_soap.pdf", width=6, height=6)


# model predictions for water drain test positivity
null_value_wdrain = mean(predict(fit, newdata=data.frame(predgrid_wdrain,
                                                         Sample.type="Water drain"), 
                                 type="link")) # avg logit(prop positive) across market among Water drain samples
plogis(null_value_wdrain) 
# 0.01215741

preds_wdrain_df = predictions(fit, newdata=data.frame(predgrid_wdrain,
                                                      Sample.type="Water drain"), 
                              vcov=T, conf_level=0.95, type="link") %>%
  transform(estimate = 100*plogis(estimate),
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

graph2pdf(file="./plots/SARSCoV2_positivity_market_wdrain_binomial_GAM_soap.pdf", width=6, height=6)


# TO DO: maybe try to fit models including the centered logits of the prop of reads
# mapping to particular host species as covariates
# or
# maybe better: test for enrichment of particular host species among SARS-CoV2
# pos samples using multinomial models or multinomial mixed models
# (nnet::multinom or mclogit::mblogit - the latter also allows one to take
# into account overdispersion)




# PART BELOW SOME TESTS, DON'T LOOK AT THIS ####


# compositional data PCA

comp_data_transp = t(tableS8[,-1])
colnames(comp_data_transp) = tableS8[,1]
rownames(comp_data_transp) = colnames(tableS8)[-1]
comp_data_transp_clr = clr(comp_data_transp) # Aitchison centered log-ratio transformed
sample = rownames(comp_data_transp)
genus = colnames(comp_data_transp)
group = factor(tableS1$SARS.CoV.2.qPCR.result[match(sample, tableS1$Sample.ID)],
               levels=c("Negative","Positive"))

pca1 = prcomp(comp_data_transp_clr, scale=FALSE)

# Extract scores and loadings
scores = as.data.frame(pca1$x)
loadings = as.data.frame(pca1$rotation)
scores_df = data.frame(PC1 = scores[, 1], PC2 = scores[, 2])
scores_df$group = group
loadings_df = data.frame(Variable = row.names(loadings), 
                         PC1 = 50*loadings[, 1], PC2 = 50*loadings[, 2])

# Create the biplot using ggplot2
ggplot(scores_df, aes(x = PC1, y = PC2, colour = group)) +
  geom_point() +
  # geom_text_repel(aes(label = row.names(scores_df)), size = 3) +
  geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm")), alpha = 0.75, color = "blue") +
  geom_text(data = loadings_df, aes(x = PC1, y = PC2, label = Variable), size = 3, color = "blue") +
  theme_few() +
  scale_colour_manual("", values = c("blue", "red")) +
  labs(title = "", x = "PC1", y = "PC2") # +
# theme(legend.position = "none")


ggbiplot::ggbiplot(pca1, choices=c(2,3), groups=group,
                   obs.scale=1, var.scale=1, ellipse=F, circle=T, varname.size = 2.5
) + 
  geom_mark_hull(aes(colour=group, fill=group),
                 expand = unit(0, "mm"),
                 radius = unit(0, "mm"),
                 concavity=11) +
  scale_fill_manual("", values=c("blue","red")) +
  scale_colour_manual(guide=F, values=c("blue","red")) +
  ggtitle("Biplot") + 
  coord_cartesian(xlim=c(-25,25), ylim=c(-25,25)) 
# theme_few(base_size = 20) # +
# theme(legend.position = c(0.8, 0.2), legend.background=element_blank())

graph2ppt(file='biplot.pptx', width=10, height=14)


pheatmap(comp_data_transp_clr)

annotation = data.frame(group = group)
rownames(annotation) <- rownames(comp_data_transp_clr)
pheatmap(comp_data_transp_clr, 
         scale="none",
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         annotation_row=annotation)