if (!require(agricolae)) install.packages('agricolae',dependencies = T)
if (!require(nlme)) install.packages('nlme',dependencies = T)
if (!require(multcomp)) install.packages('mulcomp',dependencies = T)
if (!require(emmeans)) install.packages('emmeans',dependencies = T)
if (!require(ContaminatedMixt)) install.packages('ContaminatedMixt',dependencies=T)
if (!require(EnvStats)) install.packages('EnvStats',dependencies = T)
if (!require(ggplot2)) install.packages('ggplot2',dependencies = T)

mult.t.test <- function(dataC, level = 0.05) {
  
  groups <- levels(as.factor(dataC$group))
  lb <- length(unique(dataC$block))
  comb <- combn(x = groups, m = 2)
  
  Truth <- cont <- c()
  
  for (i in seq(ncol(comb))) {
    
    dataG <- subset(dataC, group %in% comb[,i])
    
    pval <- c()
    for (b in unique(dataC$block)) {
      
      dataB <- subset(dataG, block == b)
      data.norm <- shapiro.test(dataB$yield)$p.value>0.05
      if (data.norm) s <- t.test(yield~group, data = dataB) else s <- wilcox.test(yield~group, data = dataB)
      # s <- t.test(yield~group, data = dataB)
      pval <- c(pval, s$p.value<level)
      
    }
    Truth <- c(Truth, ifelse(mean(pval)>0.5,TRUE,FALSE))
    cont <- c(cont, paste(comb[1,i],comb[2,i],sep=' - '))
    
  }
  
  table.test <- cbind.data.frame(signficance = Truth)
  rownames(table.test) <- cont
  return(table.test)
}

extract.decision <- function(dataC, post.hoc, level = 0.05) {
  # list.post <- c("bh","holm","by","hoch","hom","snk", 
  #                "hsd","bonf","sidak","regwq","lsd","dunc","schef","dunnet")
  # kind.post <- c(rep('stepwise', 6), rep('simultaneous', 7), 'control')
  post.emmeans <- c('sidak', 'hommel', 'hoch','by', 'bh', 'holm', 'bonf', 'tukey')
  post.agricolae <- c('snk', 'lsd', 'duncan', 'scheffe', 'regwq')
  post.glht <- 'dunnet'
  
  model <- nlme::lme(yield ~ group,random =~1|block, data = dataC)
  lme.aov <- aov(yield ~ group+block, data = dataC)
  if (post.hoc %in% post.emmeans) {
    
    
    p <- summary(emmeans(model, pairwise ~ group, adjust = post.hoc))[['contrasts']][,c(1,6)]
    post.decision <- cbind.data.frame(decision = p$p.value<level, row.names = p$contrast)
    
  } else if (post.hoc == 'snk') {
    
    p <- SNK.test(lme.aov, "group", group = F, console = F)[['comparison']]['pvalue']
    post.decision <- cbind.data.frame(decision = p<level)
    
  } else if (post.hoc == 'lsd') {
    
    p <- LSD.test(lme.aov, "group", group = F, console = F)[['comparison']]['pvalue']
    post.decision <- cbind.data.frame(decision = p<level)
    
  } else if (post.hoc == 'scheffe') {
    
    p <- scheffe.test(lme.aov, "group", group = F, console = F)[['comparison']]['pvalue']
    post.decision <- cbind.data.frame(decision = p<level)
    
  } else if (post.hoc == 'duncan') {
    
    p <- duncan.test(lme.aov, "group", group = F, console = F)[['comparison']]['pvalue']
    post.decision <- cbind.data.frame(decision = p<level)
    
  } else if (post.hoc == 'regwq') {
    
    p <- regwq(lme.aov, "group", group = F, console = F)[['comparison']]['pvalue']
    post.decision <- cbind.data.frame(decision = p<level)
    
  } else {
    
    p <- summary(glht(model, linfct = mcp(group = "Dunnett"),
                      test = adjusted('none')))[['test']][['pvalues']]
    post.decision <- cbind.data.frame(decision = p<level)
    groups <- levels(as.factor(dataC$group))
    comb <- combn(x = groups, m = 2); cont <- c()
    for (i in 1:(length(groups)-1)) cont <- c(cont, paste(comb[1,i],comb[2,i],sep=' - '))
    rownames(post.decision) <- cont
  } 
  return(post.decision)
  
}

comp.pcp <- function(dataC, pcp, level = 0.05){
  
  postR <- extract.decision(dataC = dataC, post.hoc = pcp,level = level)
  goodS <- mult.t.test(dataC = dataC, level = level)
  if (pcp == 'dunnet') goodS <- goodS[1:nrow(postR),1];postR <- postR[,1]
  return(c(sum(postR), sum(postR[goodS==T])))
}


# simulate (sim) databases using parameters

sDataFrame <- function(t = 4,          # The number of group
                       ni= rep(5,4),   # The group size vector or scalar
                       d = "N",        # The distribution (N, E, T, S)
                       m = 2+1:4*0.5,  # The group means vectors 
                       v = rep(0.83*t+0.25, t) # The group variance vectors
) {
  # dummy variables
  x.dummy <- function(x = c('G2','G3'), fact.level = paste('G',1:4,sep='')) {
    l <- length(fact.level)
    m <- length(x)
    dummy <- c()
    for (j in seq(m)){
      for (i in seq(l)){
        dummy <- c(dummy, ifelse(test = x[j] == fact.level[i], yes = 1, no = 0))
      }
    }
    return(matrix(dummy, m, l, byrow = T))
  }
  # random values form distribution specified
  random <- function(n, dist = 'N', p1 = 0, p2 = 0.05){
    
    if (dist == "T"){
      
      r <- EnvStats::rnormTrunc(n, p1, p2, min = -1, max = 2)
      
    } else if (dist == "E") {
      
      r <- rexp(n, rate = 200)
      
    } else if (dist == "C") {
      
      r <- as.numeric(ContaminatedMixt::rCN(n, p1, p2, .95, 2))
      
      
    } else {
      
      r <- rnorm(n, p1, p2)
      
    }
    r
  }
  
  group <- as.factor(rep(paste("G", rep(1:t,ni), sep = ""),5))
  lev.group <- levels(group)
  blEff <- rnorm(5,0,1.5)
  block <- as.factor(rep(paste("B",1:5, sep =""), each = sum(ni)))
  lev.block <- levels(block)
  GroupEff <- as.vector(mvtnorm::rmvnorm(1,m, diag(rep(0.005,t))))
  x <- x.dummy(x = group, fact.level = lev.group)
  z <- x.dummy(x = block, fact.level = lev.block)
  epsi <- c()
  for (b in 1:5) {
    for (g in 1:t) {
      epsi <- c(epsi, random(ni[g], d,p1 = 0, p2 = v[g]))
    }
  }
  yield <- x%*%GroupEff + z%*%blEff + epsi
  
  xbase <- data.frame(block = block,
                      group = group,
                      yield = yield)
}

combine <- function (size, times = 0) {
  recombine <- function(size) {
    for (i in 1:length(size)) {
      if(i == 1) {
        comb <- paste(size[i], size[i+1], sep = "/")
      } else if (i == 2){
        next
      } else {
        comb <- paste(comb, size[i], sep ="/")
      }
    }
    comb
  }
  
  if (times == 0) {
    comb <- recombine(size)
  } else {
    comb <- c()
    for (i in size){
      comb <- c(comb, recombine(rep(i,times)))
    }
  }
  return(comb)
}

# Data bases generation ---------------------------------------------------

# This function create a list of all configurations of parameters

post.hoc.perf <- function(ng = c(3:5),                 # number of groups
                          sim = 10,                    # number of simulations
                          gs = c(3,5,10),              # group sizes vector
                          dist = c("N", "E", "T", "C"), # distribution
                          level = 0.05
){
  totsim <- length(ng)*((length(gs)+1)*(6+4))*sim
  
  postname <- c("GN","Bal","GS","DT","Mean","Variance","Sim","nb","hsd","xhsd","bonf",
                "xbonf","bh","xbh","holm","xholm","by","xby","hoch","xhoch",
                "hom","xhom","sidak","xsidak","regwq","xregwq","snk","xsnk",
                "lsd","xlsd","dunc","xdunc","schef","xschef","dunnet","xdunnet")
  
  NormRes <- data.frame(matrix(NA,nrow = totsim, ncol =  length(postname)))
  names(NormRes) <- postname
  source('devPosthoc.R')
  options(warn = -1)
  
  pb <- txtProgressBar(min = 0, max = totsim, style = 3)
  d <- 0;BG <- c("Bg", "Ubg");EM <- c("EM", "UM");EV <- c("EV", "UV")
  vsim <- paste("Sim", 1:sim)
  for (i in ng) {
    
    for (j in 1:2) {
      
      if (j==1) siz <- gs else siz <- 1
      
      for (k in siz) {
        
        if (j == 2) {
          numb<- sample(gs, i, replace = ifelse(i<=length(gs),FALSE,TRUE)); Gsiz <- combine(numb)
        } else {
          numb <- rep(k, i); Gsiz <- combine(k,i)
        }
        
        for (l in dist) {
          
          for (m in 1:2) {
            
            if (m == 1) mg <- rep(2, i) else mg <- 2+1:i*1
            
            for (n in 1:2) {
              
              if (l %in% c('C', 'E', 'T') & n == 2) next
              if (l == "N" & n == 2) {
                vg <- sample(seq(0.05,2,0.05)*i,i,replace = ifelse(i<=20,FALSE,TRUE))
              } else {vg <- rep(0.83*i+0.25,i)}
              va <- ifelse(test = l == 'N', yes = EV[n], no = "ND")
              for (s in 1:sim) {
                d <- d+1
                dataC <- sDataFrame(t = i,
                                    ni = numb,
                                    d = l,
                                    m = mg,
                                    v = vg)
                
                model <- nlme::lme(yield ~ group,random =~1|block, data = dataC)
                lme.aov <- aov(yield ~ group+block, data = dataC)
                NormRes[d,] <- c(i,BG[j],Gsiz,l,sum(mult.t.test(dataC)),va,vsim[s],sum(numb),
                                 comp.pcp(dataC = dataC, 'tukey'),
                                 comp.pcp(dataC = dataC, 'bonf'),
                                 comp.pcp(dataC = dataC, 'bh'),
                                 comp.pcp(dataC = dataC, 'holm'),
                                 comp.pcp(dataC = dataC, 'by'),
                                 comp.pcp(dataC = dataC, 'hoch'),
                                 comp.pcp(dataC = dataC, 'hommel'),
                                 comp.pcp(dataC = dataC, 'sidak'),
                                 comp.pcp(dataC = dataC, 'regwq'),
                                 comp.pcp(dataC = dataC, 'snk'),
                                 comp.pcp(dataC = dataC, 'lsd'),
                                 comp.pcp(dataC = dataC, 'duncan'),
                                 comp.pcp(dataC = dataC, 'scheffe'),
                                 comp.pcp(dataC = dataC, 'dunnet')
                )
                
                
              }
              setTxtProgressBar(pb, d)
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  close(pb)
  for (i in names(NormRes)[c(1,5,8:36)]) {
    NormRes[,i] <- as.numeric(NormRes[,i])
  }
  return(NormRes)
}
memory.limit(10.030e+12)
options(warn = 0)
# system.time(dset <- post.hoc.perf(ng = sort(c(3, 5, 7, 9, 11, 13, 15, 20, 25, 30),T),
#                                   gs = sort(c(5, 10, 15, 30, 50, 70, 100, 200),T),
#                                   sim = 100))


# Measures estimates ---------------------------------------------


criteria <- function (dset) {
  
  measures <- function(dset, pcp = 'snk') {
    d.col.name <- names(dset)
    data.pcp <- dset[,d.col.name[c(1:8,which(d.col.name%in%c(pcp,paste('x',pcp,sep=''))))]]
    RH <- data.pcp[,pcp]>0
    AH <- !RH
    H0 <- data.pcp$Mean == 0
    H1 <- !H0
    CF <- data.pcp[,paste('x',pcp,sep='')]>0
    all.pow <- sum(RH & H1)/sum(H1)
    any.pow <- sum(CF & RH & H1)/sum(data.pcp$Mean[H1])
    err.t1 <- sum(RH & H0)/sum(H0)
    fdr <- sum(RH & H0 | AH & H1)/sum(H0|H1)
    return(list(APP = all.pow, ANP = any.pow, ALPHA = err.t1, FDR = fdr))
  }
  
  info.df <- function(data) {
    gnvect <- levels(as.factor(data$GN))
    bvect <- levels(as.factor(data$Bal))
    gsvect <- levels(as.factor(data$GS))
    distvect <- levels(as.factor(data$DT))
    varvect <- levels(as.factor(data$Variance))
    sim <- length(levels(as.factor(data$Sim)))
    n <- nrow(data)
    post <- names(data)[9:36]
    return(list(gnvect=gnvect,bvect=bvect,gsvect=gsvect,distvect=distvect,
                varvect=varvect,sim=sim,post=post,n=n))
  }
  
  perf.info <- info.df(dset)
  estimates <- vector(mode = 'list', length = length(perf.info$distvect))
  names(estimates) <- perf.info$distvect
  for (i in names(estimates)) {
    # Separate data into different distributions
    estimates[[i]] <- vector(mode = 'list', length = 2)
    names(estimates[[i]]) <- c('data', 'criteria')
    
    estimates[[i]][['data']] <- na.omit(dset[dset$DT == i,])
    
    # --------------------
    d1 <- estimates[[i]][['data']]
    # if (i == 'N') {h1 <- H1.extract(d1);h0 <- Ho.extract(d1)} else {h1 <- nH1.extract(d1);
    # h0 <- nHo.extract(d1)}
    # h1 <- nH1.extract(d1);h0 <- nHo.extract(d1)
    # perf.info <- info.df(h1)
    power <- unique(d1[,c(1:4, 6,8)])
    pow <- data.frame(power,matrix(NA, nrow = nrow(power), 6))
    postType <- cbind.data.frame(post = c("hsd","schef","sidak","bonf","lsd",
                                          "dunc","regwq","bh",'by',"holm","hoch",
                                          "hom","snk","dunnet"),
                                 kind = rep(c('simultaneous','stepwise', 'control'),
                                            c(5, 8, 1)))
    c <- 0
    
    for (d in 1:nrow(power)) {
      for (p in perf.info$post[1:length(perf.info$post) %%2 == 1]) {
        c<-c+1
        dsetc <- d1[d1$GN==power$GN[d]&d1$Bal==power$Bal[d]&d1$GS==power$GS[d]
                    &d1$DT==power$DT[d]&d1$Variance==power$Variance[d],]
        me <- measures(dsetc,p)
        pow[c,] <- c(power[d, ], p,postType$kind[postType$post==p],
                     me$APP,me$ANP,me$ALPHA,me$FDR)
      }
    }
    
    names(pow) <- c(names(pow)[c(1:6)], 'post.hoc', 'kind', 'APP', 'ANP', 'ALPHA','FDR')
    rownames(pow) <- NULL
    
    estimates[[i]][['criteria']] <- pow
  }
  return(estimates)
}

dset <- post.hoc.perf(ng = sort(c(3, 5, 7, 9, 15),T),
                                    gs = sort(c(10, 15, 25, 30, 50),T),
                                    sim = 1000)

meas.estim <-criteria(dset = dset)

gen.gsize <- function(meas.estim) {
  
  g1 <- meas.estim$GS
  g2 <- meas.estim$GN
  
  lite.func <- function(g1) {
    g1 <- paste(unique(unlist(strsplit(g1,'/'))),collapse = '/')
  }
  g1 <- apply(as.data.frame(g1), 1, FUN = lite.func)
  return(paste(g2, '(', g1, ')', sep=''))
  
}

plotting <- function(meas.estim, distri = names(meas.estim), kind = levels(as.factor(meas.estim[[1]][[2]]$kind))){
  
  
  gsize <- sort(unique(gen.gsize(meas.estim$C$data)))[c(7,8,11,9,10,12,13:15,18,16,17,19:21,23,22,24:26,28,27,29:30,1:6)]
  save.plots <- function (pownorm, i = i, k = k, gsize = gsize){
    options(warn = -1)
    p <- ggplot(pownorm, aes(x = x, y = APP, group = post.hoc))
    p+geom_line(aes(linetype = post.hoc, color = post.hoc))+
      geom_point(aes(shape = post.hoc, color = post.hoc)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme(axis.text.x = element_text(face="bold", angle=90)) +
      scale_x_discrete(name ="Group size",limits=gsize) +
      geom_hline(yintercept=.8, linetype="dashed",
                 color = "red", size = 1)+
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
               ymax = .8, alpha = .1, fill = "red") +
      ggtitle(label = paste(i,"Distribution", "APP estimates", k, collapse = ' - '))+
      facet_grid(Variance~.)
    ggsave(plot = last_plot(), filename = paste('plots/',i,k,'APP.png',sep = ''),width = 7, height = 8)
    
    # ANPP
    
    p1 <- ggplot(pownorm, aes(x = x, y = ANP, group = post.hoc))
    p1+geom_line(aes(linetype = post.hoc, color = post.hoc))+
      geom_point(aes(shape = post.hoc, color = post.hoc)) +
      coord_cartesian(ylim = c(0, 1)) +
      theme(axis.text.x = element_text(face="bold", angle=90)) +
      scale_x_discrete(name ="Group size",limits=gsize) +
      geom_hline(yintercept=.8, linetype="dashed",
                 color = "red", size = 1)+
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
               ymax = .8, alpha = .1, fill = "red") +
      ggtitle(label = paste(i,"Distribution", "ANP estimates", k, collapse = ' - '))+
      facet_grid(Variance~.)
    ggsave(plot = last_plot(), filename = paste('plots/',i,k,'ANP.png',sep = ''),width = 7, height = 8)
    
    # ALPHA
    e <- ggplot(pownorm, aes(x = x, y = ALPHA, group = post.hoc))
    e +geom_line(aes(linetype = post.hoc, color = post.hoc))+
      geom_point(aes(shape = post.hoc, color = post.hoc)) +
      coord_cartesian(ylim = c(0, .5)) +
      theme(axis.text.x = element_text(face="bold", angle=90, size = 8)) +
      scale_x_discrete(name ="Group size",limits=gsize) +
      geom_hline(yintercept=.05, linetype="dashed",
                 color = "red", size = 1)+
      geom_hline(yintercept=0.0635, linetype="dashed",
                 color = "blue", size = 0.5)+
      geom_hline(yintercept=0.0365, linetype="dashed",
                 color = "blue", size = 0.5)+
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 0.0365,
               ymax = 0.0635, alpha = .2, fill = "green") +
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
               ymax = 0.0365, alpha = .3, fill = "yellow") +
      ggtitle(label = paste(i,"Distribution", "ALPHA estimates", k, collapse = ' - '))+
      facet_grid(Variance~.)
    ggsave(plot = last_plot(), filename = paste('plots/',i,k,'ALPHA.png',sep = ''),width = 7, height = 8)
    # False Discovery Error rates
    
    f <- ggplot(pownorm, aes(x = x, y = FDR, group = post.hoc))
    f +geom_line(aes(linetype = post.hoc, color = post.hoc), size = .5)+
      geom_point(aes(shape = post.hoc, color = post.hoc)) +
      coord_cartesian(ylim = c(0, .5)) +
      theme(axis.text.x = element_text(face="bold", angle=90, size = 8)) +
      scale_x_discrete(name ="Group size",limits=gsize) +
      geom_hline(yintercept=.05, linetype="dashed",
                 color = "red", size = 1)+
      geom_hline(yintercept=0.0635, linetype="dashed",
                 color = "blue", size = 0.5)+
      geom_hline(yintercept=0.0365, linetype="dashed",
                 color = "blue", size = 0.5)+
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
               ymax = 0.0635, alpha = .2, fill = "green") +
      ggtitle(label = paste(i,"Distribution", "FDR estimates", k, collapse = ' - '))+
      facet_grid(Variance~.)
    ggsave(plot = last_plot(), filename = paste('plots/',i,k,'FDR.png',sep = ''),width = 7, height = 8)
  }
  
  for (i in distri) {
    pownorm2 <- meas.estim[[i]][[2]][,c(1,3,5:12)]
    
    for (k in kind){
      pownorm <- pownorm2[pownorm2$kind == k,]
      pownorm$x <- as.factor(gen.gsize(pownorm))
      save.plots(pownorm, i, k,gsize)
    }
    
  }
  
}

plotting(meas.estim=meas.estim)

# Simultaneous and stepwise PCPs
N.step <- subset(meas.estim$N$criteria, kind == 'stepwise')
for (i in names(N.step)[c(2:8)]) {
  N.step[,i] <- as.factor(N.step[,i])
}
# Power
# ANP
summary(aov(ANP~GN+nb+Variance+post.hoc, data = N.step))

# APP
summary(aov(APP~Bal*GN*nb+Variance+post.hoc, data = N.step))

# ALPHA
summary(aov(ALPHA~Bal*GN*nb+Variance+post.hoc, data = N.step))

# FDR
summary(aov(FDR~GN*Bal*nb+Variance+post.hoc, data = N.step))


N.sim <- subset(meas.estim$N$criteria, kind == 'simultaneous')
for (i in names(N.sim)[c(1:8)]) {
  N.sim[,i] <- as.factor(N.sim[,i])
}
# Power
# ANP
summary(aov(ANP~GN+nb+Variance+post.hoc, data = N.sim))

# APP
summary(aov(APP~Bal*GN*nb+Variance+post.hoc, data = N.sim))

# ALPHA
summary(aov(ALPHA~Bal*GN*nb+Variance+post.hoc, data = N.sim))

# FDR
summary(aov(FDR~GN*Bal*nb+Variance+post.hoc, data = N.sim))




# Contaminated

C.sim <- subset(meas.estim$C$criteria, kind == 'simultaneous')
for (i in names(C.sim)[c(1:8)]) {
  C.sim[,i] <- as.factor(C.sim[,i])
}
summary(aov(ANP~Bal*GN*nb+post.hoc, data = C.ss.anp))

# Change box plot colors by groups
ggplot(N.ss.anp, aes(x=nb, y=ANP, fill=GN)) +
  geom_line()


# ALPHA
N.ss.alp <- subset(meas.estim$N$criteria[,c(1:8,11)], kind != 'control')
for (i in names(N.ss.alp)[c(1:8)]) {
  N.ss.alp[,i] <- as.factor(N.ss.alp[,i])
}
summary(aov(ALPHA~Bal*GN*nb+Variance+post.hoc, data = N.ss.alp))


# Change box plot colors by groups
ggplot(N.ss.alp, aes(x=as.factor(Bal), y=ALPHA, fill=as.factor(GN))) +
  geom_boxplot()

# FDR
N.ss.fdr <- subset(meas.estim$N$criteria[,c(1:8,12)], kind != 'control')
for (i in names(N.ss.fdr)[c(1:8)]) {
  N.ss.fdr[,i] <- as.factor(N.ss.fdr[,i])
}

# Change box plot colors by groups
ggplot(N.ss.fdr, aes(x=Bal, y=FDR, fill=GN)) +
  geom_boxplot()


# Control
N.con.anp <- subset(meas.estim$N$criteria[,c(1:8,10)], kind == 'control')
for (i in names(N.con.anp)[c(1:8)]) {
  N.con.anp[,i] <- as.factor(N.con.anp[,i])
}
summary(aov(ANP~Bal+Variance+GN+GS, data = N.con.anp))

# Change box plot colors by groups
ggplot(N.con.anp, aes(x=Variance, y=ANP, fill=GN)) +
  geom_boxplot()


N.con.app <- subset(meas.estim$N$criteria[,c(1:9)], kind == 'control')
for (i in names(N.con.app)[c(1:8)]) {
  N.con.app[,i] <- as.factor(N.con.app[,i])
}
summary(aov(APP~Bal+Variance+GN+GS, data = N.con.app))

# Change box plot colors by groups
ggplot(N.con.app, aes(x=Variance, y=APP, fill=GN)) +
  geom_boxplot()

N.con.alp <- subset(meas.estim$N$criteria[,c(1:8,11)], kind == 'control')
for (i in names(N.con.alp)[c(1:8)]) {
  N.con.alp[,i] <- as.factor(N.con.alp[,i])
}
summary(aov(ALPHA~Bal+Variance+GN+GS, data = N.con.alp))

# Change box plot colors by groups
ggplot(N.con.alp, aes(x=Bal, y=ALPHA, fill=GN)) +
  geom_boxplot()

N.con.fdr <- subset(meas.estim$N$criteria[,c(1:8,12)], kind == 'control')
for (i in names(N.con.fdr)[c(1:8)]) {
  N.con.fdr[,i] <- as.factor(N.con.fdr[,i])
}
summary(aov(FDR~Bal+Variance+GN+GS, data = N.con.fdr))

# Change box plot colors by groups
ggplot(N.con.fdr, aes(x=Variance, y=FDR, fill=GN)) +
  geom_boxplot()

C.con.app <- subset(meas.estim[[4]]$criteria[,c(1:9,10)], kind == 'control')
for (i in names(C.con.app)[c(1:8)]) {
  C.con.app[,i] <- as.factor(C.con.app[,i])
}
summary(aov(APP~GN*GS, data = C.con.app))

# Change box plot colors by groups
ggplot(C.con.app, aes(x=Bal, y=APP, fill=GN)) +
  geom_boxplot()
