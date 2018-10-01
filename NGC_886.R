require(misc3d)
library(plotly)
require(caret)
data <- read.csv("Zoombie.csv") 



x <- data$ra
y <- data$dec
z <- data$parallax

p <- plot_ly(data, x = ~bp_rp,  y = ~parallax, z = ~-g,color = ~parallax) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'BP-RP'),
                      yaxis = list(title = 'parallax [mas]'),
                      zaxis = list(autorange = "reversed",title = 'G')))




plot_ly(data, x = ~pmra,  y = ~parallax, z = ~pmdec,color = ~parallax) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'pmra'),
                      yaxis = list(title = 'parallax [mas]'),
                      zaxis = list(title = 'pmdec')))



hist(data$parallax)
hist(spatialSign(data$parallax))