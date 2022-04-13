library(tidyverse)

cheese <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/cheese.csv")

# Process data to have a log of price instead
cheese$logP <- log(cheese$price)
cheese$disp <- as.factor(cheese$disp)

# Scatter plot with log P on x access and volume on y, with display as color
cheese %>% 
  ggplot(aes(x=logP, y=vol, color=disp)) +
  geom_point(size=0.5, alpha=0.5)
ggsave("./Desktop/Sp22/modeling/SDS383c/ex05/cheese_eda.png")

# Looks like with display is heteroskedastic and with display is homoskedastic
# Just going to use heteroskedastic generally then






