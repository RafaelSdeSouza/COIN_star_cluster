require(dplyr)
require(magrittr)
library(factoextra)
require(mclust)
library(spatstat)
require(VBmix)
require(MixAll)
require(EMMIXcskew)
require(caret)

Star_dat <- read.csv("bad_fsr_0735_large_noproba.csv",header = TRUE) %>% 
  #select(.,-"proba") %>%
            filter(.,phot_g_mean_mag < 20)  %>% 
            na.omit() %>% select(c("ra","dec","parallax","pmra","pmdec","phot_g_mean_mag",
           "phot_bp_mean_mag","phot_rp_mean_mag")) %>%
  mutate(.,bp.rp = phot_bp_mean_mag - phot_rp_mean_mag) %>%
  mutate(.,Nparallax = spatialSign(parallax))

features <- c("ra","dec","pmra","pmdec","Nparallax")

Star_CL <- Mclust(Star_dat[,features],G=1:40)
Star_labels <- Star_CL$classification
Star_dat$labels <- as.factor(Star_labels)


mat <- matrix(c(iris[,1],iris[,2]), nrow=nrow(iris), byrow=TRUE)

Image3 <- readData(iris[,1:2])

radon(as.Image(iris[,1:2]))

rad = radon(t(im))$rData



ggplot(Star_dat,aes(x=pmra,y=pmdec,group=labels,color=labels,alpha=0.1)) + 
  geom_point() + 
#  facet_wrap(.~labels) +
  theme_bw() + theme(legend.position =  "none") 
#+
#  geom_density2d()  


plot_ly(Star_dat, x = ~pmra,  y = ~parallax, z = ~-pmdec,color = ~labels) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'pmra'),
                      yaxis = list(title = 'parallax [mas]'),
                      zaxis = list(title = 'pmdec')))



ggplot(Star_dat,aes(x= bp.rp,y=phot_g_mean_mag,group=labels,color = parallax, alpha=0.1)) + 
  geom_point(size=0.2) + facet_wrap(.~labels) +
  theme_bw() + theme() +
  #  geom_density2d() +
  scale_y_reverse() + scale_color_viridis_c()


ggplot(Star_dat,aes(x=pmra,y=pmdec,group=labels,color=labels,alpha=0.1)) + 
  geom_point() + 
  #  facet_wrap(.~labels) +
  theme_light() + theme(legend.position =  "none") +
  geom_voronoi_tile(aes(fill = labels), normalize = TRUE) +
  geom_voronoi_segment(normalize = TRUE) +
  scale_fill_viridis_d() +
  scale_color_viridis_d()


ggplot(Star_dat,aes(x=ra,y=dec,group=labels,color=labels,alpha=0.1)) + 
#  geom_point() +
  facet_wrap(.~labels) +
  theme_bw() + theme(legend.position =  "none") +
  geom_density2d()






ggplot(Star_dat,aes(x= bp.rp,y=phot_g_mean_mag)) + 
  geom_point(size=0.2) +
  theme_bw() + theme() +
  #  geom_density2d() +
  scale_y_reverse() + scale_color_viridis_c()



temp_dat <-  Star_dat
i <- 1
for (i in 1:10){
  Star_CL <- Mclust(temp_dat[,features],G=1:20)
  Star_labels <- Star_CL$classification
  temp_dat$labels <- as.factor(Star_labels)
  
  surv <- which(sqrt(Star_CL$parameters$variance$scale) <= 1 )
  temp_dat <- temp_dat[Star_CL$classification %in%surv,]
  i <- i  + 1  
}


spatial_test <- function(dat){
#dat <- apply(dat, 2, rescale) 
star.win <- owin(range(dat[,1]), range(dat[,2]))  
star <- as.ppp(dat, star.win)
Jc <- Jest(star, 0.01)
tt <- max(Jc$km - Jc$theo,na.rm=TRUE)
return(tt < 0.1)
}


library(RandomFields)
# 1d time series
n <- 256
rf <- GaussRF(x = c(0,1, 1/n), model = "stable", 
              grid = TRUE, gridtriple = TRUE,
              param = c(mean=0, variance=1, nugget=0, scale=1, kappa=1))
par(mfrow=c(4,2))


dataX <- rnorm(1500, 0, 0.1)
dataY <- rnorm(1500, 0, 0.1)
Rxyc <- data.frame(dataX,dataY)
z <- kde2d(Rxyc$dataX,Rxyc$dataY)$z
fd.estim.squareincr(z, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)



dataX <- runif(1500, 0, 1)
dataY <- runif(1500, 0,1)
Rxyr <- data.frame(dataX,dataY)
z <- kde2d(Rxyr$dataX,Rxyr$dataY)$z
fd.estim.squareincr(z, p.index = 1, plot.loglog = TRUE, plot.allpoints = TRUE)



dat = Rxyr
data(redwood)
J <- Jest(redwood, 0.01)
plot(J, main="redwood data")
# values are below J= 1, indicating clustered pattern

X <- rpoispp(lambda = 50,win = owin(c(0, 10), c(0, 10)))
plot(X)
Kc <- Kest(X)
Jc <- Jest(X)

plot(Kc)
E <- envelope(X, Kest, nsim = 500)
plot(E, main = NULL)


plot(nndensity(X))



spatial_test(Rxyc)




dataX <- runif(15, 0, 1)
dataY <- runif(15, 0, 1)
Rxyr <- data.frame(dataX,dataY)
spatial_test(Rxyr)


dat <- Rxyc
star.win <- owin(range(dat[,1]), range(dat[,2]))  
star <- as.ppp(dat, star.win)

Js <- Jest(star,0.1)
plot(Js)
mad.test(star, Jest, nsim = 100,correction = "border")



data(cells)
G <- Gest(cells)
plot(G)

plot(G$km-G$theo)
# P-P style plot
plot(G, cbind(km,theo) ~ theo)


plot(nndensity(star))

mad.test(star, Kest, verbose=FALSE, nsim=100,
         scale=function(r) { r })


fd2d <- fd.estimate(star, methods="filter1",
                    window.size = 100, step.size=100, plot.loglog = TRUE)

E <- envelope(star, Kest, nsim = 500)
plot(E, main = NULL)


dataX <- runif(15, 0, 1)
dataY <- runif(15, 0, 1)
Rxyr <- data.frame(dataX,dataY)

#spatial_test(Rxyr)

dat <- Rxyr
star.win <- owin(range(dat[,1]), range(dat[,2]))  
star <- as.ppp(dat, star.win)
dclf.test(star, Kest, nsim = 100)
ppm.bei <- ppm(star)
plot(predict.ppm(ppm.bei, type = "lambda"))


plot(dmixpois(1:100, 10, 5, invlink = I))


hopskel(redwood)
hopskel.test(star, alternative="clustered",method="MonteCarlo")

spatial_test(Rxyc)


dat <- apply(Rxyc, 2, rescale) 
star.win <- owin(range(dat[,1]), range(dat[,2]))  
star <- as.ppp(dat, star.win)
dclf.test(star, Jest, nsim = 100)

0.00016386/0.00036546


mad.test(star, Kest, nsim = 100)




spatial_test(Rxyr)




spatial_test(Rxyc)

tt <- NULL
for(i in 1:100){
dataX <- runif(15, 0, 1)
dataY <- runif(15, 0, 1)
Rxyr <- data.frame(dataX,dataY)
ss <- spatial_test(Rxyr)
tt <- rbind(tt,ss$statistic$mad)
}

tt <- NULL
for(i in 1:100){
dataX <- rnorm(15, 0, 1)
dataY <- rnorm(15, 0, 1)
Rxyc <- data.frame(dataX,dataY)
ss <- spatial_test(Rxyc)
tt <- rbind(tt,ss$statistic$mad)
}
stc <- spatial_test(Rxyc)



gal.win <- owin(range(Rxyc[,1]), range(Rxyc[,2]))  
gal <- as.ppp(Rxyc, gal.win)
#plot(Gest(gal))




E <- envelope(gal, Kest, nsim = 500)
cor(E$r,E$theo)
ad.test(E$theo,E$r)

plot(E, main = NULL)

#Significance level of pointwise Monte Carlo test: 2/101 = 0.0198



# Envelope of L function under CSR
  L(r) = sqrt(K(r)/pi)
## Not run: 
E <- envelope(X, Kest)
plot(E, sqrt(./pi) ~ r)

str <- spatial_test(Rxyr)


spatial_test(iris[,1:2])

centroid.owin(gal.win) ; area.owin(gal.win)
gal <- as.ppp(shap.lo[,c(1,2,4)], gal.win)  # ppp = spatstat's planar point pattern
summary(gal)
plot(density(gal, adjust=0.5), col=topo.colors(20), main='')# kernel density estimator
plot(gal, lwd=2, add=TRUE) # symbol size scaled to velocity marks

# Hypothesis tests for complete spatial randomness

clarkevans.test(gal, correction='donnelly', alternative='clustered') # Clark-Evans test
hopskel.test(gal, alternative='clustered') # Hopkins-Skellam test
mad.test(gal, interpolate=TRUE)   # median absolute deviation test
dclf.test(gal, interpolate=TRUE)  # Diggle-Cressie-Loosmore-Ford test
dg.test(gal, nsim=50, use.theory=TRUE, interpolate=TRUE)   #  Dao-Genton test
plot(dg.envelope(gal, Lest, nsim=50))  # Simulation envelope for DG test
