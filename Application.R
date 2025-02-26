# EM on Classification Model
# Work on HPV trajectories
# Code by Dr Maxime Bonjour
# Helped by Dr Hadrien Charvat & Damien Georges

# Based on : http://tinyheero.github.io/2016/01/03/gmm-em.html

# last version : 20/03/2024

# setwd()

## Libraries
library("tidyverse")
library("lme4")
library("LaplacesDemon")
library("statmod")
library("numDeriv")
library("RColorBrewer")
library("cowplot")

## Functions
source("HPV.function.classi3.R")

## Import the data
temp1 <- readRDS('data.application.rds')

## Trajectory Classification 
### with fixed effect (without heterogenity):
#### 2 groups
Traj1.2 <- MyTrajEM.fixed(data = temp1, formula=~ 1+time+time2,num="value",denom="denom",clust="id",nb.grp=2,inter.grp=c("time","time2"),max.iter=200,simul=500,seed=200421)
save(Traj1.2, file = "Traj1.2.RData")

#### 3 groups
Traj1.3 <- MyTrajEM.fixed(data = temp1, formula=~1+time+time2,num="value",denom="denom",clust="id",nb.grp=3,inter.grp=c("time","time2"),max.iter=200,simul=500,seed=200421)
save(Traj1.3, file = "Traj1.3.RData")

#### 4 groups
Traj1.4 <- MyTrajEM.fixed(data = temp1, formula=~1+time+time2,num="value",denom="denom",clust="id",nb.grp=4,inter.grp=c("time","time2"),max.iter=200,simul=500,seed=200421)
save(Traj1.4, file = "Traj1.4.RData")


### with mixed effect (with heterogenity):
#### 2 groups
Traj2.2 <- MyTrajEM.mixed(data = temp1, formula=~1+time+time2,num="value",denom="denom",clust="id",nb.grp=2,inter.grp=c("time","time2"),max.iter=200,simul=50,seed=200421,rdm.eff="time",nb.aghq=1)
save(Traj2.2, file = "Traj2.2.RData")

#### 3 groups
Traj2.3 <- MyTrajEM.mixed(data = temp1, formula=~1+time+time2,num="value",denom="denom",clust="id",nb.grp=3,inter.grp=c("time","time2"),max.iter=200,simul=50,seed=200421,rdm.eff="time",nb.aghq=1)
save(Traj2.3, file = "Traj2.3.RData")

#### 4 groups
Traj2.4 <- MyTrajEM.mixed(data = temp1, formula=~1+time+time2,num="value",denom="denom",clust="id",nb.grp=4,inter.grp=c("time","time2"),max.iter=200,simul=50,seed=200421,rdm.eff="time",nb.aghq=1)
save(Traj2.4, file = "Traj2.4.RData")


# Or then load the data

# load(file = "Traj1.2.RData")
# load(file = "Traj1.3.RData")
# load(file = "Traj1.4.RData")

# load(file = "Traj2.2.RData")
# load(file = "Traj2.3.RData")
# load(file = "Traj2.4.RData")



## Models characteristics
### with fixed effect:
#### 2 groups
Traj1.2$final.prop.group     # proportion of trajectories in each cluster
Traj1.2$ICL.BIC              #ICL_BIC indicator
Traj1.2$glm1 %>% BIC         #BIC indicator

#### 3 groups
Traj1.3$final.prop.group
Traj1.3$ICL.BIC
Traj1.3$glm1 %>% BIC

#### 4 groups
Traj1.4$final.prop.group
Traj1.4$ICL.BIC
Traj1.4$glm1 %>% BIC

### with mixed effect:
#### 2 groups
Traj2.2$final.prop.group
Traj2.2$ICL.BIC
Traj2.2$glmer1 %>% BIC

#### 3 groups
Traj2.3$final.prop.group
Traj2.3$ICL.BIC
Traj2.3$glmer1 %>% BIC

#### 4 groups
Traj2.4$final.prop.group
Traj2.4$ICL.BIC
Traj2.4$glmer1 %>% BIC


## Models figures

Graphic representation of the clusters trajectories identified

age.intercept <- 30

### With fixed effects
a <- graph.result.1(Traj1.2, age.intercept = age.intercept)
b <- graph.result.1(Traj1.3, age.intercept = age.intercept)
c <- graph.result.1(Traj1.4, age.intercept = age.intercept)

plot_grid(a,b,c, nrow = 1)

### With mixed effects

d <- graph.result.1(Traj2.2, age.intercept = age.intercept)
e <- graph.result.1(Traj2.3, age.intercept = age.intercept)
f <- graph.result.1(Traj2.4, age.intercept = age.intercept)

plot_grid(d,e,f, nrow = 1)

