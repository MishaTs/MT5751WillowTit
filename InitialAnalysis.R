library(statsecol)
library(unmarked)
library(tidyverse)
library(MuMIn)
library(AICcmodavg)
library(viridisLite)
library(latex2exp)
# read in data
data(willow)
# prepare the data how I need it
willowNum <- willow %>% mutate(forestsq = forest^2,
                               iLength = 1/length) %>%
  mutate_all(as.numeric) %>% 
  rownames_to_column("id")

willowUnm <- unmarkedFrameOccu(
  y = willowNum[,c("y.1","y.2","y.3")], 
  # we are interested in breeding bird occupancy for the WHOLE breeding season, so days (which affects occupancy on a more granular level) are not good site covariates in this case
  siteCovs = data.frame(elev = willowNum$elev,
                        elev2 = willowNum$elevsq,
                        forest = willowNum$forest,
                        forest2 = willowNum$forestsq,
                        iLength = willowNum$iLength), 
  obsCovs = list(day = willowNum[,c("day1","day2","day3")], 
                 dur = willowNum[,c("dur1","dur2","dur3")],
                 intensity = willowNum[,c("intensity1","intensity2","intensity3")],
                 length = willowNum[,c("length","length","length")],
                 iLength = willowNum[,c("iLength","iLength","iLength")],
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

mAlt <- occu(formula = ~day + dur + forest  # p formula
              ~elev + forest, #psi formula
              data = willowUnm) #the data object

mAlt1 <- occu(formula = ~day + dur + forest # p formula
              ~elev + elev2 + forest + forest2, #psi formula
              data = willowUnm) #the data object

mAlt2 <- occu(formula = ~day + dur + forest + forest2 # p formula
              ~elev + elev2 + forest + forest2, #psi formula
              data = willowUnm) #the data object

mOptm_Alt <- occu(formula = ~day + dur + forest + forest*day # p formula
              ~elev, #psi formula
              data = willowUnm) #the data object

# recording visibility/leafing may be a good way to increase separability and avoid forest in both models
mOptm <- occu(formula = ~day + dur + forest + forest*day # p formula
              ~elev + elev2 + forest + forest2 + elev2*forest + elev2*forest2 , #psi formula
              data = willowUnm) #the data object

mAlt3 <- occu(formula = ~day + dur + forest + forest*day # p formula
              ~elev + elev2 + forest + forest2 + elev*forest + elev*forest2 + elev2*forest + elev2*forest2 , #psi formula
              data = willowUnm) #the data object

# combine some representative models for displaying
fl <- fitList(
  "p(.)                                   psi(.)"                                                                      = m0,
  "p(day + dur + forest)                  psi(elev + forest)"                                                          = mAlt,
  "p(day + dur + forest)                  psi(elev + elev^2 + forest + forest^2)"                                      = mAlt1,
  "p(day + dur + forest + forest^2)       psi(elev + elev^2 + forest + forest^2)"                                      = mAlt2,
  "p(day + dur + forest + forest*day)     psi(elev)"                                                                   = mOptm_Alt,  
  "p(day + dur + forest + forest*day)     psi(elev + elev^2*(forest + forest^2))"                                      = mOptm,
  "p(day + dur + forest + forest*day)     psi((elev + elev^2)*(forest + forest^2))"                                    = mAlt3)

# model output table to format
ms <- modSel(fl)

# actual model summary
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
        legend.position = "bottom") + 
  labs(x = "",
       y = "",
       fill = "Elevation (m)") + 
  guides(fill = guide_colorbar(# draw border around the legend
                               frame.colour = "black",
                               barwidth = 10)) + 
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
        legend.position="bottom") + 
  labs(x = "",
       y = "",
       fill = "Forest (%)")  + 
  guides(fill = guide_colorbar(# draw border around the legend
    frame.colour = "black",
    barwidth = 10)) + 
  coord_fixed()

for_pred <- data.frame(elev = (Switzerland$elevation - 1182.574)/646.333,      # convert original m to z-score
                       elev2 = ((Switzerland$elevation - 1182.574)/646.333)^2, # convert original m to z-score
                       forest = Switzerland$forest/100,          #want prop not %
                       forest2 = Switzerland$forest/100,         #want prop not %
                       X = Switzerland$x,                        #keep the coordinates
                       Y = Switzerland$y)                        #keep the coordinates
cowplot::plot_grid(gelev,gfor,nrow=2)

willowPredSDM <- modavgPred(list(mOptm), # top model
                            newdata = for_pred, #spatially indexed data frame
                            parm.type = "psi",  #predict from state model
                            c.hat = gof.boot$c.hat.est) #inflate SEs using Royle & Kery method
#add data to predictions manually
willow_sdm <- for_pred %>% mutate(Predicted = willowPredSDM$mod.avg.pred,
                                  SE = willowPredSDM$uncond.se,
                                  lower = willowPredSDM$lower.CL,
                                  upper = willowPredSDM$upper.CL)

gpredM_1 <- ggplot(data = willow_sdm, aes(x=X, y=Y,fill=Predicted)) +
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
        legend.position="bottom") + 
  labs(x = "",
       y = "",
       fill = TeX(r'(Estimated $\psi$)'))  + 
  guides(fill = guide_colorbar(# draw border around the legend
    frame.colour = "black",
    barwidth = 10)) + 
  coord_fixed()

gpredE <- ggplot(data = willow_sdm, aes(x=X, y=Y,fill=SE)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "H") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="bottom") + 
  labs(x = "",
       y = "",
       fill = TeX(r'(Estimated $\psi$ Error)'))  + 
  guides(fill = guide_colorbar(# draw border around the legend
    frame.colour = "black",
    barwidth = 10)) + 
  coord_fixed()

cowplot::plot_grid(gelev,gfor,gpredM_1,gpredE,nrow=2)

gpredL <- ggplot(data = willow_sdm, aes(x=X, y=Y,fill=lower)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "H", 
                       limits= c(0,1)) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="left") + 
  labs(x = "",
       y = "",
       fill = TeX(r'(Lower Bound $\psi$)')) + 
  guides(fill = guide_colorbar(# draw border around the legend
    frame.colour = "black")) + 
  coord_fixed()

gpredM_2 <- ggplot(data = willow_sdm, aes(x=X, y=Y,fill=Predicted)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "H", 
                       limits= c(0,1)) +
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
       fill = TeX(r'(Estimated $\psi$)'))  + 
  guides(fill = guide_colorbar(# draw border around the legend
    frame.colour = "black")) + 
  coord_fixed()

gpredH <- ggplot(data = willow_sdm, aes(x=X, y=Y,fill=upper)) +
  geom_raster() +
  scale_fill_viridis_c(direction = 1, 
                       option = "H", 
                       limits= c(0,1)) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position="left") + 
  labs(x = "",
       y = "",
       fill = TeX(r'(Upper Bound $\psi$)'))  + 
  guides(fill = guide_colorbar(# draw border around the legend
    frame.colour = "black")) + 
  coord_fixed()

cowplot::plot_grid(gpredL, gpredM_2, gpredH,nrow=3)



#---------------------------------------------------------------------------------------#
#psi ~ elev | mean(forest)
pred_psi_eleL <- data.frame(elev = seq(min(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      max(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      length = 30),
                           forest = quantile(probs = 0.25, willowUnm@siteCovs$forest, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
predPsiEleL <- modavgPred(list(mOptm), newdata = pred_psi_eleL, parm.type = "psi", c.hat = gof.boot$c.hat.est)
pred_psi_eleL <- pred_psi_eleL %>% mutate(Predicted = predPsiEleL$mod.avg.pred,
                                          SE = predPsiEleL$uncond.se,
                                          lower = predPsiEleL$lower.CL,
                                          upper = predPsiEleL$upper.CL,
                                          elevR = 1182.574 + elev*646.333)
ggpsieleL <- ggplot(data = pred_psi_eleL, aes(x = elevR, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("P(Occupied)") + xlab("Elevation (m)") + ylim(0,1) + theme_bw()

pred_psi_eleM <- data.frame(elev = seq(min(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      max(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      length = 30),
                           forest = median(willowUnm@siteCovs$forest, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
predPsiEleM <- modavgPred(list(mOptm), newdata = pred_psi_eleM, parm.type = "psi", c.hat = gof.boot$c.hat.est)
pred_psi_eleM <- pred_psi_eleM %>% mutate(Predicted = predPsiEleM$mod.avg.pred,
                                          SE = predPsiEleM$uncond.se,
                                          lower = predPsiEleM$lower.CL,
                                          upper = predPsiEleM$upper.CL,
                                          elevR = 1182.574 + elev*646.333)
ggpsieleM <- ggplot(data = pred_psi_eleM, aes(x = elevR, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("P(Occupied)") + xlab("Elevation (m)") + ylim(0,1) + theme_bw()

pred_psi_eleH <- data.frame(elev = seq(min(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      max(willowUnm@siteCovs$elev, na.rm=TRUE),
                                      length = 30),
                           forest = quantile(probs = 0.75, willowUnm@siteCovs$forest, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
predPsiEleH <- modavgPred(list(mOptm), newdata = pred_psi_eleH, parm.type = "psi", c.hat = gof.boot$c.hat.est)
pred_psi_eleH <- pred_psi_eleH %>% mutate(Predicted = predPsiEleH$mod.avg.pred,
                                          SE = predPsiEleH$uncond.se,
                                          lower = predPsiEleH$lower.CL,
                                          upper = predPsiEleH$upper.CL,
                                          elevR = 1182.574 + elev*646.333)
ggpsieleH <- ggplot(data = pred_psi_eleH, aes(x = elevR, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("P(Occupied)") + xlab("Elevation (m)") + ylim(0,1) + theme_bw()

#---------------------------------------------------------------------------------------#
#psi ~ for | mean(elev)
pred_psi_forL <- data.frame(forest = seq(min(willowUnm@siteCovs$forest, na.rm=TRUE),
                                      max(willowUnm@siteCovs$forest, na.rm=TRUE),
                                      length = 30),
                           elev = quantile(probs = 0.25, willowUnm@siteCovs$elev, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
predPsiForL <- modavgPred(list(mOptm), newdata = pred_psi_forL, parm.type = "psi", c.hat = gof.boot$c.hat.est)
pred_psi_forL <- pred_psi_forL %>% mutate(Predicted = predPsiForL$mod.avg.pred,
                                          SE = predPsiForL$uncond.se,
                                          lower = predPsiForL$lower.CL,
                                          upper = predPsiForL$upper.CL,
                                          forestP = forest*100)
ggpsiforL <- ggplot(data = pred_psi_forL, aes(x = forestP, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("P(Occupied)") + xlab("Forest Cover (%)") + ylim(0,1) + theme_bw()

pred_psi_forM <- data.frame(forest = seq(min(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        max(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        length = 30),
                           elev = median(willowUnm@siteCovs$elev, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
predPsiForM <- modavgPred(list(mOptm), newdata = pred_psi_forM, parm.type = "psi", c.hat = gof.boot$c.hat.est)
pred_psi_forM <- pred_psi_forM %>% mutate(Predicted = predPsiForM$mod.avg.pred,
                                          SE = predPsiForM$uncond.se,
                                          lower = predPsiForM$lower.CL,
                                          upper = predPsiForM$upper.CL,
                                          forestP = forest*100)
ggpsiforM <- ggplot(data = pred_psi_forM, aes(x = forestP, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("P(Occupied)") + xlab("Forest Cover (%)") + ylim(0,1) + theme_bw()

pred_psi_forH <- data.frame(forest = seq(min(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        max(willowUnm@siteCovs$forest, na.rm=TRUE),
                                        length = 30),
                           elev = quantile(probs = 0.75, willowUnm@siteCovs$elev, na.rm=TRUE)) %>% 
  mutate(elev2 = elev^2,
         forest2 = forest^2)
predPsiForH <- modavgPred(list(mOptm), newdata = pred_psi_forH, parm.type = "psi", c.hat = gof.boot$c.hat.est)
pred_psi_forH <- pred_psi_forH %>% mutate(Predicted = predPsiForH$mod.avg.pred,
                                    SE = predPsiForH$uncond.se,
                                    lower = predPsiForH$lower.CL,
                                    upper = predPsiForH$upper.CL,
                                    forestP = forest*100)
ggpsiforH <- ggplot(data = pred_psi_forH, aes(x = forestP, y = Predicted)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#440154", alpha=0.1) +
  geom_line(size=1,color="#440154") +
  ylab("P(Occupied)") + xlab("Forest Cover (%)") + ylim(0,1) + theme_bw()

cowplot::plot_grid(ggpsieleL, ggpsieleM, ggpsieleH, ggpsiforL, ggpsiforM, ggpsiforH, nrow=2)















# generate predictions for 4 quadrants of interest
willowPred <- willowNum %>% filter(id %in% c(25, 62, 150, 203)) %>% 
  rename(elev2 = elevsq,
         forest2 = forestsq)
# predicting occurrence from the state process
predQuads <- modavgPred(list(mOptm), newdata = willowPred, parm.type = "psi", c.hat = gof.boot$c.hat.est)
willowRes <- willowPred %>% mutate(Predicted = predQuads$mod.avg.pred,
                                   SE = predQuads$uncond.se,
                                   lower = predQuads$lower.CL,
                                   upper = predQuads$upper.CL,
                                   elev = round((646.333*elev+1182.574),0)) %>% 
  select(-c(elev2, forest2))
