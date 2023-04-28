# 2D BINOMIAL GAM SOAP FILM SMOOTHS ####

# BINOMIAL GAM SOAP FILM SMOOTH, TAKING INTO ACCOUNT MARKET BOUNDARY
# see tutorial https://fromthebottomoftheheap.net/2016/03/27/soap-film-smoothers/
# and ?soap
# NOTE: unfortunately, these smooths are not supported by brms for the moment

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
                    family = quasibinomial(logit))
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
                    family = quasibinomial(logit))
summary(fit_gam_soap2)
GCV(fit_gam_soap1)
GCV(fit_gam_soap2) # model with interaction by Sample.type, fit_gam_soap2, fits better



# fit on aggregated per-stall data (to remove sampling bias towards stalls with high SARS-CoV2 positivity)
# for Environmental swab samples
fit_gam_soap_env_aggr = gam(cbind(SARS.CoV.2.qPCR.positive, SARS.CoV.2.qPCR.negative) ~ s(x, y, bs="so", k=5,  # soap film smooth, I chose k visually (also possible to use GCV etc)
                                                                                          xt=list(bnd=boundary_gam, nmax=150)),
                            # note: a traditional tensor spline would be ~ ti(x)+ti(y)+te(x, y), but this would not take into account
                            # market boundaries, or in context of a binomial GLM ~ ns(x, df=XX):ns(y, df=XX)
                            # predictions of all of these are quite similar
                            knots=knots,
                            data = data_aggregated_per_stall[data_aggregated_per_stall$Sample.type=="Environmental swab",],
                            family = quasibinomial(logit)) # quasi to take into account overdispersion
summary(fit_gam_soap_env_aggr)

# for Water drain samples
fit_gam_soap_wdrain_aggr = gam(cbind(SARS.CoV.2.qPCR.positive, SARS.CoV.2.qPCR.negative) ~ s(x, y, bs="so", k=5,  # soap film smooth, I chose k visually (also possible to use GCV etc)
                                                                                             xt=list(bnd=boundary_gam, nmax=150)),
                               # note: a traditional tensor spline would be ~ ti(x)+ti(y)+te(x, y), but this would not take into account
                               # market boundaries, or in context of a binomial GLM ~ ns(x, df=XX):ns(y, df=XX)
                               # predictions of all of these are quite similar
                               knots=knots,
                               data = data_aggregated_per_stall[data_aggregated_per_stall$Sample.type=="Water drain",],
                               family = quasibinomial(logit)) # quasi to take into account overdispersion
summary(fit_gam_soap_wdrain_aggr)


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
fit = fit_gam_soap_env_aggr
null_value_envswab = mean(predict(fit, newdata=data.frame(predgrid_envswabs,
                                                          Sample.type="Environmental swab"), 
                                  type="link")) # avg logit(prop positive) across market among Environmental swabs
plogis(null_value_envswab) 
# 0.05830965

preds_envswab_df = predictions(fit, newdata=data.frame(predgrid_envswabs,
                                                       Sample.type="Environmental swab"), 
                               vcov=T, conf_level=0.95, type="link") %>%
  transform(estimate = 100*plogis(predicted),
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
fit = fit_gam_soap1
null_value_wdrain = mean(predict(fit, newdata=data.frame(predgrid_wdrain,
                                                         Sample.type="Water drain"), 
                                 type="link")) # avg logit(prop positive) across market among Water drain samples
plogis(null_value_wdrain) 
# 0.06057788

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

graph2pdf(file="./plots/SARSCoV2_positivity_market_wdrain_binomial_GAM_soap.pdf", width=6, height=6)






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