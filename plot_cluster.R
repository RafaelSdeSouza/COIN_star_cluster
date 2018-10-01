require(dplyr)
require(magrittr)
library(factoextra)
require(mclust)
require('spatstat') ; library(spatstat)
require(VBmix)
require(MixAll)
require(classInt)
Star_dat <- read.csv("output_final.csv",header = TRUE)
#%>% 
  #select(.,-"proba") %>%
#  filter(.,phot_g_mean_mag < 20)  %>% 
#  na.omit() %>% select(c("ra","dec","parallax","pmra","pmdec","phot_g_mean_mag",
#                         "phot_bp_mean_mag","phot_rp_mean_mag")) %>%
#  mutate(.,bp.rp = phot_bp_mean_mag - phot_rp_mean_mag)

plot(Star_dat$l,Star_dat$b)

Van_Gogh <- c("#263C8B","#4E74A6", "#BDBF78", "#BFA524", "#2E231F")
colpal<-colorRampPalette(Van_Gogh)(100)
cutColor <- 100

colpa <- classIntervals(massinla$out[valid], n = cutColor, style= "fisher")

