require(ggpubr)
require(reshape2)
require(ggthemes)
require(dplyr)
require(ggplot2)
require(rpart)
require(rpart.plot)

train <- read.csv("../data/bazin_all_girz.csv",header = T) %>%
  filter(sample == "train") %>%
  select("g_A","g_t0","g_tfall","g_trise","g_B","i_A","i_t0","i_tfall",
         "i_trise","i_B","r_A","r_t0","r_tfall","r_trise","r_B",
         "z_A","z_t0","z_tfall","z_trise","z_B","type")  %>%
  mutate(type = replace(type,type == 0,"Ia")) %>%
  mutate(type = replace(type,type %in%
                          c(2,3,4,12,15,17,19,20,21,24,25,26,27,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44),"II")) %>%
  mutate(type = replace(type,type %in%
                          c(1,5,6,7,8,9,10,11,13,14,16,18,22,23,29,45,28),"Ibc"))%>%
  filter(type %in% c("Ia","II")) 

DT.model <- rpart(type ~., data = train)
pdf("decision_tree.pdf",height = 8,width = 15)
rpart.plot(DT.model,tweak=1.5,under = T)
dev.off()


tree <- rpart(Kyphosis ~ ., data = kyphosis)
prp(tree, extra = 7)                               # left graph
library(plotmo)
plotmo(DT.model, type = "prob",nresponse = "Ia") # middle graph
# type = "prob" is passed to predict()
plotmo(DT.model, type = "prob", nresponse = "Ia", # right graph
       type2 = "image", ngrid2 = 200,              # type2 = "image" for an image plot
       pt.col = ifelse(train$type == "Ia", "red", "lightblue"))