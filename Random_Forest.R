# Random Forest
require("plot3D")
require(randomForest)
require(viridis)

Star_dat <- read.csv("ngc_2360_withproba.csv",header = TRUE) %>% 
  #select(.,-"proba") %>%
  filter(.,phot_g_mean_mag < 18)  %>% 
  na.omit() %>% select(c("ra","dec","parallax","pmra","pmdec","phot_g_mean_mag",
                         "phot_bp_mean_mag","phot_rp_mean_mag","proba")) %>%
  mutate(.,bp.rp = phot_bp_mean_mag - phot_rp_mean_mag)






x1 <- Star_dat$bp.rp
x2 <- Star_dat$phot_g_mean_mag
y <-  Star_dat$proba

# predict values on regular xy grid
#grid.lines = 100
#x1.pred <- seq(min(x1), max(x1), length.out = grid.lines)
#x2.pred <- seq(min(x2), max(x2), length.out = grid.lines)
#x1x2 <- expand.grid( x1 = x1.pred, x2 = x2.pred)

tree_model <- randomForest(y ~x1+x2,ntree=1000)
yr <- predict(tree_model)

Star_dat$yr <- yr
ggplot(Star_dat ,aes(x= bp.rp,y=phot_g_mean_mag,color = yr, size=yr, alpha=1)) + 
  geom_point(size=0.2) + 
  theme_bw() + theme() +
  #  geom_density2d() +
  scale_y_reverse() + scale_color_viridis_c()


partition.tree(tree_model, ordvars=c("x1","x2"), add=TRUE)

yt.pred <- 1 - matrix(predict(tree_model, newdata = x1x2), 
                    nrow = grid.lines, ncol = grid.lines)



scatter3D(x1, x2, y,  cex = 1, cex.lab=1.5,type="p",pch = 19,alpha=0.5,
          theta = 50, phi = 30, ticktype = "detailed",col="red",bty = "u",
          xlab="x1",
          ylab="x2",
          zlab="y", 
          expand =0.6, 
          colkey = FALSE,
          main = paste(i," Trees"),
          surf = list(x = x1.pred, y = x2.pred, z = yt.pred,col = viridis(200), shade = 0.35,
                      lwd=0.5,lty=1))


