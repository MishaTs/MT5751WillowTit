library(statsecol)
library(unmarked)
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(viridisLite)

# read in data
data(willow)

willowNum <- willow %>% mutate(forestsq = forest^2) %>%
  mutate_all(as.numeric)

# figure out where the data is coming from
library(AHMBook)
data("SwissTits")
# lots of data cleaning for no reason
willowFilt <- data.frame(SwissTits$counts) %>% 
  select(ends_with("Willow.tit"))
willowCol <- colnames(willowFilt)
willowColNew <- str_sub(willowCol, end = 7)
willowFilt <- willowFilt %>% 
  rename_with(~ willowColNew, all_of(willowCol)) %>% 
  rownames_to_column("siteID") %>% 
  pivot_longer(
    cols = `X1.2004`:`X3.2013`, 
    names_to = "id",
    values_to = "count") %>% rowwise() %>% 
  mutate(y = paste0("y", str_sub(id, start = 2, end = 2)),
         yr = str_sub(id, start = 4)) %>% select(-id)
willowFin <- willowFilt %>%   pivot_wider(
  names_from = y, 
  values_from = count)
# data is from 2011: 237 observations
testObs <- willowFin %>% group_by(yr) %>% 
  summarise(numObs = sum(!(is.na(y1) & is.na(y2) & is.na(y3))))
# create data for plotting
willowOcc <- left_join(willowFin, SwissTits$sites, by = "siteID") %>% 
  mutate(y1 = as.numeric(y1 > 0),
         y2 = as.numeric(y2 > 0),
         y3 = as.numeric(y3 > 0)) %>% 
  filter(!if_all(c(y1, y2, y3), is.na))

# give up and try something else
# still doesn't work
siteID <- data.frame(SwissTits$sites) %>% mutate(forestP = forest/100,
                                                 elevZ = round(((elev - 1182.574)/646.333),6))
willowID <- willowNum %>% mutate(forestP = forest,
                                 rlength = length,
                                 elevZ = round(elev,6)) %>% 
  rownames_to_column("id")
willowPlot <- left_join(willowID, siteID, by = c("forestP","rlength"))



willowUnm <- unmarkedFrameOccu(
  y = willowNum[,c("y.1","y.2","y.3")], 
  # we are interested in breeding bird occupancy for the WHOLE breeding season, so days (which affects occupancy on a more granular level) are not good site covariates in this case
  siteCovs = data.frame(elev = willowNum$elev,
                        elev2 = willowNum$elevsq,
                        forest = willowNum$forest,
                        forest2 = willowNum$forestsq), 
  obsCovs = list(day = willowNum[,c("day1","day2","day3")], 
                 dur = willowNum[,c("dur1","dur2","dur3")],
                 intensity = willowNum[,c("intensity1","intensity2","intensity3")],
                 length = willowNum[,c("length","length","length")],
                 forest = willowNum[,c("forest","forest","forest")],
                 forest2 = willowNum[,c("forestsq","forestsq","forestsq")],
                 elev = willowNum[,c("elev","elev","elev")],
                 elev2 = willowNum[,c("elevsq","elevsq","elevsq")])
)
summary(willowUnm)

# null model
m0 <- occu(~1 ~1, data = willowUnm)

# full model
# intensity:length interaction is just dur and thus not included
mfull <- occu(formula = ~day + dur + intensity + length + day*dur + dur*intensity + day*length + day*intensity + dur*length  # p formula
              ~elev + elev2 + forest + forest2 + elev*forest + elev2*forest2 + elev*forest2 + elev2*forest, #psi formula
              data = willowUnm) #the data object

summary(mfull)
  
# all three agree on a few things
# elev2 and elev:forest are mutually exclusive
# p(day) is not worth including
# all p interaction terms are useless
mDredgeB <- dredge(mfull, rank = "BIC")
mDredgeA <- dredge(mfull, rank = "AIC")
mDredgeAc <- dredge(mfull, rank = "AICc")

mBIC <- occu(formula = ~dur # p formula
             ~elev + forest + forest^2 + elev*forest, #psi formula
             data = willowUnm) #the data object
summary(mBIC)

mfullAct <- occu(formula = ~day + dur + intensity + length + day*dur + dur*intensity + day*length + day*intensity + dur*length + forest + forest*dur + forest*day + forest*intensity + forest*length # p formula
              ~elev + elev2 + forest + forest2 + elev*forest + elev2*forest2 + elev*forest2 + elev2*forest, #psi formula
              data = willowUnm) #the data object

summary(mfullAct)

# all three agree on a few things
# elev2 and elev:forest are mutually exclusive
# p(day) is not worth including
# all p interaction terms are useless
mDredgeB2 <- dredge(mfullAct, rank = "BIC")
mDredgeA2 <- dredge(mfullAct, rank = "AIC")
mDredgeAc2 <- dredge(mfullAct, rank = "AICc")

mOptm <- occu(formula = ~day + dur + forest + forest*day # p formula
              ~elev + elev2 + forest + forest2 + elev2*forest + elev2*forest2 , #psi formula
              data = willowUnm) #the data object
summary(mOptm)

# test for VIF?
vif(mOptm, type = "state")
vif(mOptm, type = "det")

# GOF for best model....it's meh
gof.boot <- mb.gof.test(mOptm, nsim = 1000, ncores = 10)
ggplot() + 
  geom_histogram(data = data.frame(t.star = gof.boot$t.star),
                 aes(x=t.star), color="black", fill="#27813E", alpha = 0.3) +
  geom_vline(aes(xintercept = gof.boot$chi.square), linewidth = 0.8, color = "maroon") +
  xlab("Peasron Chi-square statistic") + theme_bw() +
  ggtitle(bquote("Bootstrapped"~chi^2~"fit statistic (1000 samples, p = 0.071)"))

# some plotting of results
data(Switzerland)
gelev <- ggplot(data = Switzerland, aes(x=x, y=y,fill=elevation)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "B") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="left") + 
  labs(x = "",
       y = "",
       fill = "Elevation (m)") + 
  coord_fixed()
gfor <- ggplot(data = Switzerland, aes(x=x, y=y,fill=forest)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "D") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="left") + 
  labs(x = "",
       y = "",
       fill = "% Forest") + 
  coord_fixed()

for_pred <- data.frame(elev = (Switzerland$elevation - 1182.574)/646.333,      # convert original m to z-score
                       elev2 = ((Switzerland$elevation - 1182.574)/646.333)^2, # convert original m to z-score
                       forest = Switzerland$forest/100,          #want prop not %
                       forest2 = Switzerland$forest/100,         #want prop not %
                       X = Switzerland$x,                        #keep the coordinates
                       Y = Switzerland$y)                        #keep the coordinates
cowplot::plot_grid(gelev,gfor,nrow=2)

willow_sdm <- predict(object = mOptm,     #the top model
                      type = "state",      #predict from state model
                      newdata = for_pred,  #spatially indexed data frame
                      append=TRUE)         #add data to predictions
gpred <- ggplot(data = willow_sdm, aes(x=X, y=Y,fill=Predicted)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "H") +
  # add actual observations if we have x,y data
  # geom_point(data = willowNum, aes(x=X, y=Y)) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="left") + 
  labs(x = "",
       y = "",
       fill = "P(Occupancy)") + 
  coord_fixed()
cowplot::plot_grid(gelev,gfor,gpred,nrow=3)



#---------------------------------------------------------------------------------------#
#psi ~ elev | mean(forest)
pred_psi_eleL <- data.frame(elev = seq(min(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      max(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      length = 30),
                           forest = quantile(probs = 0.20, willowUnm@siteCovs$forest, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
pred_psi_eleL <- predict(mOptm,type="state",newdata = pred_psi_eleL, append = TRUE) %>% 
  mutate(elevR = 1182.574 + elev*646.333) 
ggpsieleL <- ggplot(data = pred_psi_eleL, aes(x = elevR, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("Pr(Occupancy)") + xlab("Elevation (m)") + ylim(0,1) + theme_bw()

pred_psi_eleM <- data.frame(elev = seq(min(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      max(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      length = 30),
                           forest = median(willowUnm@siteCovs$forest, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
pred_psi_eleM <- predict(mOptm,type="state",newdata = pred_psi_eleM, append = TRUE) %>% 
  mutate(elevR = 1182.574 + elev*646.333) 
ggpsieleM <- ggplot(data = pred_psi_eleM, aes(x = elevR, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("Pr(Occupancy)") + xlab("Elevation (m)") + ylim(0,1) + theme_bw()

pred_psi_eleH <- data.frame(elev = seq(min(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      max(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      length = 30),
                           forest = quantile(probs = 0.80, willowUnm@siteCovs$forest, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
pred_psi_eleH <- predict(mOptm,type="state",newdata = pred_psi_eleH, append = TRUE) %>% 
  mutate(elevR = 1182.574 + elev*646.333) 
ggpsieleH <- ggplot(data = pred_psi_eleH, aes(x = elevR, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("Pr(Occupancy)") + xlab("Elevation (m)") + ylim(0,1) + theme_bw()

#---------------------------------------------------------------------------------------#
#psi ~ for | mean(elev)
pred_psi_forL <- data.frame(forest = seq(min(willowUnm@siteCovs$forest, na.rm=TRUE),
                                      max(willowUnm@siteCovs$forest, na.rm=TRUE),
                                      length = 30),
                           elev = quantile(probs = 0.20, willowUnm@siteCovs$elev, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
pred_psi_forL <- predict(mOptm,type="state",newdata = pred_psi_forL, append = TRUE)
ggpsiforL <- ggplot(data = pred_psi_forL, aes(x = forest, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("Pr(Occupancy)") + xlab("Forest Cover (%)") + ylim(0,1) + theme_bw()

pred_psi_forM <- data.frame(forest = seq(min(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        max(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        length = 30),
                           elev = median(willowUnm@siteCovs$elev, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
pred_psi_forM <- predict(mOptm,type="state",newdata = pred_psi_forM, append = TRUE)
ggpsiforM <- ggplot(data = pred_psi_forM, aes(x = forest, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("Pr(Occupancy)") + xlab("Forest Cover (%)") + ylim(0,1) + theme_bw()

pred_psi_forH <- data.frame(forest = seq(min(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        max(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        length = 30),
                           elev = quantile(probs = 0.80, willowUnm@siteCovs$elev, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
pred_psi_forH <- predict(mOptm,type="state",newdata = pred_psi_forH, append = TRUE)
ggpsiforH <- ggplot(data = pred_psi_forH, aes(x = forest, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("Pr(Occupancy)") + xlab("Forest Cover (%)") + ylim(0,1) + theme_bw()

cowplot::plot_grid(ggpsieleL, ggpsieleM, ggpsieleH, ggpsiforL, ggpsiforM, ggpsiforH, nrow=2)
