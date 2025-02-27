# EM on Classification Model
# Code by Maxime Bonjour
# & Dr Hadrien Charvat

# Based on : http://tinyheero.github.io/2016/01/03/gmm-em.html

# version 1: 20/04/2020
# version 2: 20/12/2023
# version 3: 20/03/2024

## library(nlme)
library("tidyverse")
library("lme4")
library("LaplacesDemon")
library("statmod")
library("numDeriv")
library("gghighlight")
require("maps")
library("readxl")


source("LogitAGHQ2.R")


## Load the df.sampling object
## TODO save this object as .rds
load("df.sampling.RData")

#### Algo Fixed

df.clust <- df.sampling %>%
  filter(mid.age >= 18, mid.age <= 70) %>%
  transmute(
    id = ID,
    time = mid.age,
    time2 = mid.age * mid.age,
    time3 = mid.age * mid.age * mid.age,
    value = num_HR,
    denom
  )

############ test
# temp <- df.clust
# data = temp
# formula=~1+time+time2+time3
# num="value"
# denom="denom"
# clust="id"
# nb.grp=4
# inter.grp=c("time","time2","time3")
# max.iter=50
# simul=10
# seed=200421

############
MyTrajEM.fixed <- function(data, formula = ~ 1 + time + time2, num = "value", denom = "denom", clust = "id", nb.grp = 3, inter.grp = c("time", "time2"), max.iter = 50, simul = 10, seed = 200421) {
  ## Data preparation (get cluster-specific information)
  data[, "Denum"] <- data[, denom] - data[, num]
  Clust <- unique(data[, clust])
  NbClust <- length(Clust$id)
  Num <- data[, num]
  Denom <- data[, denom]
  IdxClust <- sapply(1:NbClust, function(x) {
    which(data[, clust][[1]] == Clust[x, 1][[1]])
  }, simplify = FALSE)

  ## Model considered (partie fixe)
  ## Group-independent
  Design <- model.matrix(formula, data = data)
  fix.var <- colnames(Design)
  le.fix.var <- length(fix.var)

  ## Number of groups
  grp.names <- paste0("group.aug", 2:nb.grp)

  ## Variables in interaction with group (must be
  ## chosen among the variables in fixvar)
  inter.var <- c("(Intercept)", inter.grp)
  formul.inter <- paste0(inter.grp, collapse = "+")
  le.inter.var <- length(inter.var)
  idx.inter <- which(fix.var %in% inter.var)
  inter.var.grp <-
    matrix(
      paste0(c("", paste0(inter.var[-1], ":")), rep(grp.names, each = le.inter.var)),
      (nb.grp - 1),
      le.inter.var,
      byrow = TRUE
    )

  ## Create a "shell" table to store the group-specific
  ## cluster-specific llik
  table.vraiss <- matrix(NA, NbClust, nb.grp)

  ClustLogLik.fix <- function(i, beta) {
    Temp <- Design[IdxClust[[i]], ]
    TempN <- Num[[1]][IdxClust[[i]]]
    TempD <- Denom[[1]][IdxClust[[i]]]
    res <- sum(log(dbinom(x = TempN, size = TempD, prob = invlogit(Temp %*% beta))))
    return(res)
  }


  ## E-Step
  e_step <- function(glm.m, pi.pre) {
    ## effets fixes + calcul de la -log.like par groupe
    test.beta <- summary(glm.m)$coefficient[, 1]
    beta.vec.1 <- test.beta[fix.var]
    table.vraiss[, 1] <- sapply(1:NbClust, ClustLogLik.fix, beta = beta.vec.1)
    for (i in 2:nb.grp) {
      beta.vec.i <- beta.vec.1
      beta.vec.i[idx.inter] <- beta.vec.i[idx.inter] + test.beta[inter.var.grp[(i - 1), ]]
      table.vraiss[, i] <- sapply(1:NbClust, ClustLogLik.fix, beta = beta.vec.i)
    }

    ## selection de la -log.lkike la plus petite, attribution des nouveaux groupes
    if (length(t(table.vraiss)) %% length(pi.pre) != 0) {
      browser()
    }
    vraiss.group <- t(t(table.vraiss) + log(c(pi.pre)))
    vraiss <- apply(vraiss.group, 1, max)
    group.aug.new <- apply(vraiss.group, 1, which.max)

    # Calcul de l'ICL-BIC
    sum.vraiss <- apply(exp(table.vraiss), 1, sum)
    log.proba.appart <- vraiss.group - log(sum.vraiss)
    sum.log.proba.appart <- sum(log.proba.appart * exp(log.proba.appart))

    ## calcul de la loglikelihood totale du modele
    sum.of.comps.sum <- sum(vraiss)

    list(
      "lik" = sum.of.comps.sum,
      "group.aug.new" = as.character(group.aug.new),
      "sum.log.proba.appart" = sum.log.proba.appart
    )
  }

  ## M-Step
  m_step <- function(temp, sum.log.proba.appart) {
    temp$group.aug <- as.character(temp$group.aug)
    eval(parse(text = paste0("comp.model  <-  try(glm(cbind(", num, ",Denum) ~", as.character(formula)[2], "+(", formul.inter, ")*group.aug, data = temp, family = \"binomial\"),silent=TRUE)")))

    pi.pre <- prop.table(table(temp$group.aug))

    ICL.BIC <- BIC(comp.model) - sum.log.proba.appart

    list(
      "glm" = comp.model,
      "pi.pre" = pi.pre,
      "ICL.BIC" = ICL.BIC
    )
  }

  nb.group <- NA
  final.likelihood <- NA
  nb.iteration <- NA
  final.prop.group <- list(NA)
  final.glm <- list(NA)
  alloc.group <- list(NA)
  ICL.BIC <- NA

  ## Launch of simulations
  set.seed(seed)

  t0 <- proc.time()[3]
  for (j in 1:simul) {
    temp <- data
    grp.alea <- sample(1:nb.grp, size = NbClust, replace = TRUE)
    while (length(unique(grp.alea)) != nb.grp) {
      grp.alea <- sample(1:nb.grp, size = NbClust, replace = TRUE)
    }
    eval(parse(text = paste0("data.alea <- data.frame(", clust, "= Clust, group.aug = grp.alea)")))
    temp <- left_join(temp, data.alea, by = clust)
    temp$group.aug <- as.character(temp$group.aug)
    pi.pre.init <- table(as.character(grp.alea)) / NbClust

    ## temp <- temp %>% group_by(id) %>% summarise(group.aug = factor(sample (c(1:nb.grp),size=1, replace=TRUE), levels = c(1:nb.grp))) %>% right_join(temp, by = "id")

    ## pi.pre.init <- prop.table(table((temp %>% group_by(id) %>% slice(1) %>% select(group.aug))[,2]))  #voir pour enlever le "Adding missing groupin variables"

    eval(parse(text = paste0("glm.m.init  <-  try(glm(cbind(", num, ",Denum) ~", as.character(formula)[2], "+(", formul.inter, ")*group.aug, data = temp, family = \"binomial\"),silent=TRUE)")))

    if (inherits(glm.m.init, "try-error")) {
      nb.group.aug <- 1
    } else {
      e.step <- e_step(glm.m.init, pi.pre.init)

      eval(parse(text = paste0("data.temp <- data.frame(", clust, " = Clust, group.aug = e.step[[\"group.aug.new\"]])")))
      temp <- left_join(subset(temp, select = -group.aug), data.temp, by = clust)
      nb.group.aug <- length(unique(temp$group.aug))

      m.step <- m_step(temp, e.step[["sum.log.proba.appart"]])
      cur.lik <- e.step[["lik"]]
      lik.vector <- e.step[["lik"]]
    }


    for (i in 1:max.iter) {
      if ((nb.group.aug == nb.grp) & (!inherits(m.step[["glm"]], "try-error"))) {
        ## Repeat E and M steps till convergence
        e.step <- e_step(m.step[["glm"]], m.step[["pi.pre"]])

        eval(parse(text = paste0("data.temp <- data.frame(", clust, " = Clust, group.aug = e.step[[\"group.aug.new\"]])")))
        temp <- left_join(subset(temp, select = -group.aug), data.temp, by = clust)
        nb.group.aug <- length(unique(temp$group.aug))

        m.step <- m_step(temp, e.step[["sum.log.proba.appart"]])
        lik.vector <- c(lik.vector, e.step[["lik"]])

        if ((nb.group.aug != nb.grp) | (inherits(m.step[["glm"]], "try-error"))) {
          break
        }

        lik.diff <- abs(cur.lik - e.step[["lik"]])

        if (is.nan(lik.diff)) {

        } else {
          if (lik.diff < 0.001) {
            break
          } else {
            cur.lik <- e.step[["lik"]]
          }
        }
      } else {
        i <- 1
        break
      }
    }

    nb.group[[j]] <- length(unique(temp$group.aug))
    final.likelihood[[j]] <- max(lik.vector)
    nb.iteration[[j]] <- i
    final.prop.group[[j]] <- prop.table(table(temp$group.aug))
    final.glm[[j]] <- m.step$glm
    alloc.group[[j]] <- temp$group.aug
    ICL.BIC[[j]] <- m.step[["ICL.BIC"]]

    print(paste0("Iteration n. ", j, " - ", round(proc.time()[3] - t0, 2)))
  }

  ## Final
  IdxOK <- which(nb.group == nb.grp)
  if (length(IdxOK) > 0) {
    max.lik <- max(final.likelihood[IdxOK])
    num.final <- which(final.likelihood == max.lik)[1]

    print(final.prop.group[[num.final]])
    print(final.glm[[num.final]])

    glm1 <- final.glm[[num.final]]
    fix.effect <- summary(glm1)$coefficients[, 1]

    group.alloc <- alloc.group[[num.final]]

    final.ICL.BIC <- ICL.BIC[[num.final]]
  }

  ## Plot des resultats

  ## recuperation des effets fixes du modele selectionne et creation des
  ## reponses aux effets fixes a chaque temps
  MyTime <- seq(min(data$time), max(data$time), by = 1)
  # graph <- cbind(1,MyTime,MyTime^2)   # Hadrien

  # Moi 05.08.2024
  n.exp <- str_count(paste0(formula[2]), "time")
  graph <- matrix(
    nrow = length(MyTime),
    ncol = n.exp + 1
  )
  for (l in 0:n.exp) {
    graph[, l + 1] <- MyTime^l
  }

  MyResults <- matrix(0, length(MyTime), nb.grp)
  MyBeta <- fix.effect[fix.var]
  Temp <- exp(graph %*% MyBeta)
  MyResults[, 1] <- Temp / (1 + Temp)
  for (i in 2:nb.grp) {
    MyBeta.i <- MyBeta
    MyBeta.i[idx.inter] <- MyBeta.i[idx.inter] + fix.effect[inter.var.grp[(i - 1), ]]
    Temp <- exp(graph %*% MyBeta.i)
    MyResults[, i] <- Temp / (1 + Temp)
  }

  data$group.final <- group.alloc

  ## Iteration max atteinte
  ite.max <- sum(nb.iteration == max.iter)

  ## Bon nombre de groupe
  bon.group <- sum(nb.group == nb.grp) / simul


  print(proc.time()[3] - t0)
  return(list(data = data, final.prop.group = final.prop.group[[num.final]], glm1 = glm1, graph = cbind(MyTime, MyResults), ite.max = ite.max, bon.group = bon.group, ICL.BIC = final.ICL.BIC))
}



###### test
# temp <- df1
# data = temp1
# formula=~1+time+time2
# num="value"
# denom="denom"
# clust="id"
# nb.grp=4
# inter.grp=c("time","time2")
# max.iter=50
# simul=10
# seed=200421
# rdm.eff="time"
# nb.aghq=1


##### Algo mixed
MyTrajEM.mixed <- function(data, formula = ~ 1 + time + time2, num = "value", denom = "denom", rdm.eff = "time", clust = "id", nb.grp = 3, inter.grp = c("time", "time2"), nb.aghq = 1, max.iter = 50, simul = 10, seed = 200421, group.random = T, traces = F) {
  ## Data preparation (get cluster-specific information)
  data[, "Denum"] <- data[, denom] - data[, num]
  Clust <- unique(data[, clust])
  NbClust <- length(Clust$id)
  Num <- data[, num]
  Denom <- data[, denom]
  IdxClust <- sapply(1:NbClust, function(x) {
    which(data[, clust][[1]] == Clust[x, 1][[1]])
  }, simplify = FALSE)

  ## Model considered (partie fixe)
  ## Group-independent
  Design <- model.matrix(formula, data = data)
  fix.var <- colnames(Design)
  le.fix.var <- length(fix.var)

  ## Variables with random effect
  rdm.var <- c("(Intercept)", rdm.eff)
  formul.rdm <- paste0(rdm.eff, collapse = "+")
  n.rand <- length(rdm.var)
  pos <- which(fix.var %in% rdm.var)

  ## Number of groups
  grp.names <- paste0("group.aug", 2:nb.grp)

  ## Variables in interaction with group (must be
  ## chosen among the variables in fixvar)
  inter.var <- c("(Intercept)", inter.grp)
  formul.inter <- paste0(inter.grp, collapse = "+")
  le.inter.var <- length(inter.var)
  idx.inter <- which(fix.var %in% inter.var)
  # inter.var.grp <- matrix(paste0(c("",paste0(inter.var[-1],":")),rep(grp.names,each=nb.grp)),(nb.grp-1),le.inter.var,byrow=TRUE)
  inter.var.grp <-
    matrix(
      paste0(c("", paste0(inter.var[-1], ":")), rep(grp.names, each = le.inter.var)),
      (nb.grp - 1),
      le.inter.var,
      byrow = TRUE
    )

  ## Create a "shell" table to store the group-specific
  ## cluster-specific llik
  table.vraiss <- matrix(NA, NbClust, nb.grp)

  ClustLogLik <- function(i, beta, pos, Sigma, nbN = nb.aghq) {
    Temp <- Design[IdxClust[[i]], ]
    TempN <- Num[[1]][IdxClust[[i]]]
    TempD <- Denom[[1]][IdxClust[[i]]]
    res <- Clust.LogLik.M(beta, Temp, TempN, TempD, pos = pos, nbN = nbN, Sigma)
    return(res)
  }

  ## E-Step
  # glmer.m <- glmer.m.init
  # pi.pre <- pi.pre.init

  e_step <- function(glmer.m, pi.pre) {
    ## Matrix of random effects
    B <- matrix(bdiag(VarCorr(glmer.m)), n.rand, n.rand)

    ## effets fixes + calcul de la log.like par groupe
    test.beta <- summary(glmer.m)$coefficient[, 1]
    beta.vec.1 <- test.beta[fix.var]
    table.vraiss[, 1] <- sapply(1:NbClust, ClustLogLik, beta = beta.vec.1, pos = pos, Sigma = B)
    for (i in 2:nb.grp) {
      beta.vec.i <- beta.vec.1
      beta.vec.i[idx.inter] <- beta.vec.i[idx.inter] + test.beta[inter.var.grp[(i - 1), ]]
      table.vraiss[, i] <- sapply(1:NbClust, ClustLogLik, beta = beta.vec.i, pos = pos, Sigma = B)
    }

    ## selection de la log.lkike la plus petite, attribution des nouveaux groupes
    vraiss.group <- t(t(table.vraiss) - log(c(pi.pre)))
    vraiss <- apply(vraiss.group, 1, min)
    group.aug.new <- apply(vraiss.group, 1, which.min)


    # Calcul de l'ICL-BIC
    sum.vraiss <- apply(exp(-table.vraiss), 1, sum)
    log.proba.appart <- -vraiss.group - log(sum.vraiss)
    log.proba.appart[log.proba.appart == (Inf)] <- 0
    sum.log.proba.appart <- sum(log.proba.appart * exp(log.proba.appart))


    ## calcul de la loglikelihood totale du modele
    sum.of.comps.sum <- -sum(vraiss)

    list(
      "lik" = sum.of.comps.sum,
      "group.aug.new" = as.character(group.aug.new),
      "sum.log.proba.appart" = sum.log.proba.appart
    )
  }

  ## M-Step
  m_step <- function(temp, sum.log.proba.appart) {
    temp$group.aug <- as.character(temp$group.aug)
    eval(parse(text = paste0("comp.model  <-  try(glmer(cbind(", num, ",Denum) ~", as.character(formula)[2], "+(", formul.inter, ")*group.aug + (", formul.rdm, "|", clust, "), data = temp, family = \"binomial\"),silent=TRUE)")))

    pi.pre <- prop.table(table(temp$group.aug))

    if (!inherits(comp.model, "try-error")) {
      ICL.BIC <- BIC(comp.model) - sum.log.proba.appart
    } else {
      ICL.BIC <- NA
    }


    # if(inherits(ICL.BIC, 'try-error')){browser()}

    list(
      "glmer" = comp.model,
      "pi.pre" = pi.pre,
      "ICL.BIC" = ICL.BIC
    )
  }

  nb.group <- NA
  final.likelihood <- NA
  nb.iteration <- NA
  final.prop.group <- list(NA)
  final.glmer <- list(NA)
  alloc.group <- list(NA)
  ICL.BIC <- NA

  ## Launch of simulations
  set.seed(seed)

  t0 <- proc.time()[3]
  for (j in 1:simul) {
    temp <- data
    grp.alea <- sample(1:nb.grp, size = NbClust, replace = TRUE)

    eval(parse(text = paste0("data.alea <- data.frame(", clust, "= Clust, group.aug = grp.alea)")))
    temp <- left_join(temp, data.alea, by = clust)
    # save(temp, file = "temp.data.debug.RData")
    # load("temp.data.debug.RData")
    temp$group.aug <- as.character(temp$group.aug)
    pi.pre.init <- table(as.character(grp.alea)) / NbClust


    eval(parse(text = paste0("glmer.m.init  <-  try(glmer(cbind(", num, ",Denum) ~", as.character(formula)[2], "+(", formul.inter, ")*group.aug + (", formul.rdm, "|", clust, "), data = temp, family = \"binomial\"),silent=TRUE)")))

    if (inherits(glmer.m.init, "try-error")) {
      nb.group.aug <- 1
    } else {
      e.step <- e_step(glmer.m.init, pi.pre.init)

      ## AJOUT MF ##
      if (traces == T) {
        cat("\n", "__________________________", "\n")
        cat("\n", "glmer initital", "\n")
        print(glmer.m.init)

        cat("\n", "pi pre initital", "\n")
        print(pi.pre.init)
      }
      ## AJOUT MF ##

      eval(parse(text = paste0("data.temp <- data.frame(", clust, " = Clust, group.aug = e.step[[\"group.aug.new\"]])")))
      temp <- left_join(subset(temp, select = -group.aug), data.temp, by = clust)
      nb.group.aug <- length(unique(temp$group.aug))

      m.step <- m_step(temp, e.step[["sum.log.proba.appart"]])
      cur.lik <- e.step[["lik"]]
      lik.vector <- e.step[["lik"]]
    }

    if ((nb.group.aug == nb.grp) & (!inherits(m.step[["glmer"]], "try-error"))) {
      for (i in 1:max.iter) {
        ## Repeat E and M steps till convergence
        e.step <- e_step(m.step[["glmer"]], m.step[["pi.pre"]])

        eval(parse(text = paste0("data.temp <- data.frame(", clust, " = Clust, group.aug = e.step[[\"group.aug.new\"]])")))
        temp <- left_join(subset(temp, select = -group.aug), data.temp, by = clust)
        nb.group.aug <- length(unique(temp$group.aug))

        m.step <- m_step(temp, e.step[["sum.log.proba.appart"]])
        lik.vector <- c(lik.vector, e.step[["lik"]])

        if ((nb.group.aug != nb.grp) | (inherits(m.step[["glmer"]], "try-error"))) {
          break
        }

        lik.diff <- abs(cur.lik - e.step[["lik"]])

        ## AJOUT MF ##
        if (traces == T) {
          cat("\n", "__________________________", "\n")
          cat("\n", "iter =", i, "\n")


          cat("\n", "current lik", "\n")
          print(cur.lik)

          cat("\n", "e.step[['lik']]", "\n")
          print(e.step[["lik"]])
        }


        if (is.nan(lik.diff)) {

        } else { ## AJOUT MF ##
          if (lik.diff < 0.001) {
            break
          } else {
            cur.lik <- e.step[["lik"]]
          }
        }
      }
    } else {
      i <- 1
    }

    nb.group[[j]] <- length(unique(temp$group.aug))
    final.likelihood[[j]] <- max(lik.vector)
    nb.iteration[[j]] <- i
    final.prop.group[[j]] <- prop.table(table(temp$group.aug))
    final.glmer[[j]] <- m.step$glmer
    alloc.group[[j]] <- temp$group.aug
    ICL.BIC[[j]] <- m.step[["ICL.BIC"]]

    print(paste0("Iteration n. ", j, " - ", round(proc.time()[3] - t0, 2)))
  }

  ## Final
  IdxOK <- which(nb.group == nb.grp)
  if (length(IdxOK) > 0) {
    max.lik <- max(final.likelihood[IdxOK])
    num.final <- which(final.likelihood == max.lik)[1]

    print(final.prop.group[[num.final]])
    print(final.glmer[[num.final]])

    glmer1 <- final.glmer[[num.final]]
    fix.effect <- summary(glmer1)$coefficients[, 1]

    group.alloc <- alloc.group[[num.final]]
    final.ICL.BIC <- ICL.BIC[[num.final]]
  }

  ## Plot des resultats

  ## recuperation des effets fixes du modele selectionne et creation des
  ## reponses aux effets fixes a chaque temps
  MyTime <- seq(min(data$time), max(data$time), by = 1)
  # graph <- cbind(1,MyTime,MyTime^2)  # Hadrien

  # Moi 05.08.2024
  n.exp <- str_count(paste0(formula[2]), "time")
  graph <- matrix(
    nrow = length(MyTime),
    ncol = n.exp + 1
  )
  for (l in 0:n.exp) {
    graph[, l + 1] <- MyTime^l
  }
  #

  MyResults <- matrix(0, length(MyTime), nb.grp)
  MyBeta <- fix.effect[fix.var]
  Temp <- exp(graph %*% MyBeta)
  MyResults[, 1] <- Temp / (1 + Temp)
  for (i in 2:nb.grp) {
    MyBeta.i <- MyBeta
    MyBeta.i[idx.inter] <- MyBeta.i[idx.inter] + fix.effect[inter.var.grp[(i - 1), ]]
    Temp <- exp(graph %*% MyBeta.i)
    MyResults[, i] <- Temp / (1 + Temp)
  }


  data$group.final <- group.alloc

  ## Iteration max atteinte
  ite.max <- sum(nb.iteration == max.iter)

  ## Nombre de groupes voulu
  bon.group <- sum(nb.group == nb.grp) / simul


  print(proc.time()[3] - t0)
  return(list(data = data, final.prop.group = final.prop.group[[num.final]], glmer1 = glmer1, graph = cbind(MyTime, MyResults), ite.max = ite.max, bon.group = bon.group, ICL.BIC = final.ICL.BIC))
}


#### Function Graph
graph.result.1 <- function(Traj, age.intercept = 15) {
  NbGrp <- Traj[["final.prop.group"]] %>% length()
  MyCol <- RColorBrewer::brewer.pal(9, "Set1")
  pays.group <- Traj[["data"]] %>%
    mutate(group.final = group.final) %>%
    group_by(id) %>%
    slice(1) %>%
    select(id, group.final)
  dat <- list()
  for (i in 1:NbGrp) {
    dat[[i]] <- data.frame(MyTime = Traj[["graph"]][, "MyTime"] + age.intercept, Result = Traj[["graph"]][, i + 1])
  }


  Traj$data$age_cat <- as.factor(cut(Traj$data$time + 30, seq(15, 60, by = 5), dig.lab = 0, include.lowest = TRUE))
  levels.age <- Traj$data$age_cat %>% levels()

  mean.age <- NA
  for (i in 1:length(levels.age)) {
    mean.age[i] <- mean(c(as.numeric(substr(levels.age[i], start = 2, stop = 3)), as.numeric(substr(levels.age[i], start = 5, stop = 6))))
  }
  Traj$data$age_cat <- as.numeric(as.character(factor(Traj$data$age_cat, labels = mean.age)))
  Traj.graph <- Traj$data %>%
    group_by(id, age_cat) %>%
    summarise(num = sum(value), denom = sum(denom), group.final = unique(group.final))




  p <- ggplot() +
    geom_line(data = Traj.graph, aes(x = age_cat, y = num / denom, group = id, col = group.final)) +
    scale_color_manual(values = MyCol) +
    theme_bw() +
    labs(x = "Age", y = "HPV Prevalence")
  i <- 1
  while (i <= NbGrp) {
    df <- dat[[i]]
    p <- p + geom_line(data = df, aes(x = MyTime, y = Result), col = MyCol[i], linewidth = 1.5)
    i <- i + 1
  }



  temp2 <- Traj[["graph"]] %>% as.data.frame()

  if (NbGrp == 5) {
    names(temp2) <- c("MyTime", "1", "2", "3", "4", "5")
    temp2 <- temp2 %>% gather("group.final", "value", 2:6)
  }

  if (NbGrp == 4) {
    names(temp2) <- c("MyTime", "1", "2", "3", "4")
    temp2 <- temp2 %>% gather("group.final", "value", 2:5)
  }

  if (NbGrp == 3) {
    names(temp2) <- c("MyTime", "1", "2", "3")
    temp2 <- temp2 %>% gather("group.final", "value", 2:4)
  }

  if (NbGrp == 2) {
    names(temp2) <- c("MyTime", "1", "2")
    temp2 <- temp2 %>% gather("group.final", "value", 2:3)
  }

  p <- p + theme(legend.position = "none")

  q <- ggplot() +
    geom_line(data = Traj.graph, aes(x = age_cat, y = num / denom, group = id, col = group.final)) +
    geom_line(data = temp2, aes(x = MyTime + 30, y = value, group = group.final), linewidth = 1, col = "black", linetype = "solid") +
    theme_bw() +
    scale_color_manual(values = MyCol) +
    theme(legend.position = "none") +
    labs(x = "Age", y = "Prevalence", col = "Groups allocation") +
    facet_grid(. ~ group.final) +
    gghighlight(use_direct_label = F)


  print(p)
}
