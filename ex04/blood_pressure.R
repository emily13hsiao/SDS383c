library(dplyr)

bp <- read.csv("./Desktop/Sp22/modeling/SDS383D-master/data/bloodpressure.csv")

# (A) naive t-test

treatment <- bp %>% filter(treatment == 1)
control <- bp %>% filter(treatment == 2)

t.test(treatment$systolic, control$systolic)

# (B) mean t-test

avg_bp <- bp %>% 
  group_by(subject, .drop=FALSE) %>% 
  summarise(measurements = n(), mean_bp = mean(systolic), treatment=first(treatment))

treatment_mean <- avg_bp %>% filter(treatment==1)
control_mean <- avg_bp %>% filter(treatment==2)

t.test(treatment_mean$mean_bp, control_mean$mean_bp)




