######## Loading necessary libraries ########################################
library(strataG)
library(effects)
library(adegenet)
library(car)
library(multcomp)
library(MuMIn)
library(snow)
#library(viridis) #colour palette
library(Rfast)
options(max.print=10000)
library(extrafont)
library(extrafontdb)
setEPS()
######## FUNCTION Shannon.gen MODIFIED FROM SpadeR ########################
# Usage:
# Shannon.gen(gInd, estimator = c("ch","cs","mle","z"))
#
## Arguments:
## gInd         a genind object (the function was tested only on microsatellite data)
## estimator    a character string, a vector or a list of H estimators that
##              will be estimated:
##              - "ch" Chao et al 2013 (default),
##              - "cs" Chao and Shen 2003,
##              - "mle" original Shannon diversity index,
##              - "z" Zahl 1975.
##
## Details:
##  For the citations of the estimators see Konopiński MK (submitted)
##  The functions MLE, Z, CS and Ch were copied from Diversity function
##  from SpadeR package (https://cran.r-project.org/web/packages/SpadeR/)
##
## Values:
##  A list of data frames (one per each estimator) presenting selected
##  Shannon H estimators for all loci (loci in rows, populations in columns)
##
Shannon.gen = function(gInd, estimator = NULL) {
  MLE = function(X) {
    X = X[X > 0]
    n = sum(X)
    - sum(X / n * log(X / n))
  }
  Z = function(X){
    X=X[X>0]
    Y=X[X>1]
    n=sum(X)
    -n*sum(X/n*log(X/n))-(n-1)/n*sum((n-X)*(-X/(n-1)*log(X/(n-1))) )-(n-1)/n*sum(-Y*(Y-1)/(n-1)*log((Y-1)/(n-1)))
  }
  CS = function(X) {
    x = X
    x = x[x > 0]
    n = sum(x)
    f1 = sum(x == 1)
    C_head = 1 - f1 / n
    a = -sum(C_head * (x / n) * log(C_head * (x / n)) / (1 - (1 - C_head *
          (x / n)) ^ n))
    a
  }
  Ch = function(X) {
    x = X
    x = x[x > 0]
    n = sum(x)
    UE <- sum(x / n * (digamma(n) - digamma(x)))
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if (f1 > 0)
    {
      A <-
        1 - ifelse(f2 > 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2), (n - 1) * f1 /
                     ((n - 1) * f1 + 2))
      B = sum(x == 1) / n * (1 - A) ^ (-n + 1) * (-log(A) - sum(sapply(1:(n - 1), function(k) {
        1 / k * (1 - A) ^ k})))
    }
    if (f1 == 0) {
      B = 0
    }
    if (f1 == 1 & f2 == 0) {
      B = 0
    }
    UE + B
  }
  inputs <- lapply(seppop(gInd),
                   function (i)
                     sapply(seploc(i),
    function(j)
   colSums(j@tab,na.rm = TRUE)))
  out <- list()
  if (is.null(estimator))
    estimator <- "ch"
  if ("s" %in% estimator) {
    out[["Shannon_1949"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, MLE)))
  }
  if ("z" %in% estimator) {
    out[["Zahl_1975"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, Z)))
  }
  if ("cs" %in% estimator) {
    out[["Chao_Shen_2003"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, CS)))
  }
  if ("ch" %in% estimator) {
    out[["Chao_et_al_2013"]] <-
      as.data.frame(lapply(inputs, function (i)
        sapply(i, Ch)))
  }
  return(out)
}

######## FUNCTION SIMULATING SAMPLING OF THE POPULATION ########
## Arguments:
##  i - a genind object
##  testSize - numeric vector of sample sizes to be generated from the original
##  genind object - maximum value cannot exceed nuber of samples in the smalles
##  population in genind object
## values:
##  A 4-dimmensional array containing estimates of the four Shannon H estimators [,,x,]
##  in all loci [x,,,], populations [,x,,] and sample sizes [,,,x]
##
varTable <- function (i, testSize) {
  #creating array to store calculated indices
  # [x,,,] locus
  # [,x,,] population
  # [,,x,] index of genetic variation: (1) H' - Shannon's index of diversity, (2) Jacknife H' estimator,
  # (3) Chao & Shen 2003, (4) Chao et. al 2013
  # [,,,x] sample size
  everything <-
    array(
      data = NA,
      dim = c(nLoc(i), nPop(i), 4, length(testSize)),
      dimnames = list(
        c(locNames(i)),
        c(popNames(i)),
        c("H(MLE)", "H(Z)", "H(CS)", "H(Chao)"),
        testSize
      )
    )
  for (j in testSize) {
    resampled <- c()
    for (popul in popNames(i))
      resampled <- append(resampled,tryCatch(sample(which(i$pop==popul),j,replace = FALSE),
                                             error = function(x) print ("problem with sampling")))
    g <- i[resampled]
    Shannon_results  <-
      Shannon.gen(g,estimator = c("s", "z", "cs", "ch"))
    for (m in 1:length(Shannon_results)) {
      everything[, , m, paste0(j)] <- as.matrix(Shannon_results[[m]])
    }
  }
  return(everything)
}

######## FUNCTION CALCULATING GENETIC DIVERSITY IN THE WHOLE POPULATION #######
PopEstimates <- function (g,He=FALSE) {
  ## Arguments:
  ##  g - a genind object
  ##
  ## values:
  ##  A 3-dimmensional array containing estimates of the four Shannon H estimators [,,x]
  ##  in all loci [x,,] and populations [,x,]
  ##
  # creating array to store calculated estimators
  # [x,,] locus
  # [,x,] population
  # [,,x] index of genetic variation: (1) H' - Shannon's index of diversity, (2) Jacknife H' estimator,
  # (3) Chao & Shen 2003, (4) Chao et. al 2013,(5) He - expected heterozygosity,
  everything <- array(
    data = NA,
    dim = c(nLoc(g), nPop(g), 4),
    dimnames = list(
      c(locNames(g)),
      popNames(g),
      c("H(MLE)", "H(Z)", "H(CS)", "H(Chao)")
    )
  ) 
  if (He==TRUE){
    everything <- array(
      data = NA,
      dim = c(nLoc(g), nPop(g), 5),
      dimnames = list(
        c(locNames(g)),
        popNames(g),
        c("H(MLE)", "H(Z)", "H(CS)", "H(Chao)", "He")
      )
    ) 
    everything[,,5]  <- as.matrix(as.data.frame(lapply(seppop(g),function(i) unlist(lapply(seploc(i),Hs)))))
  }
  Shannon_results  <-
    Shannon.gen(g, estimator = c("s", "z", "cs", "ch"))
  for (m in 1:length(Shannon_results)) {
    everything[, , m] <- as.matrix(Shannon_results[[m]])
  }
  return(everything)
}


######## FUNCTION CALCULATING SAMPLING ERROR #######
SamplingError <- function (gind,ssizes = c(5,20,80),repl = 1000){
  npops <- nPop(gind)
  results <- array(data = NA, dim = c(npops,4,length(ssizes)),
                   dimnames = list(popNames(gind),c("Hmle","Hz","Hcs","Hchao"), ssizes))
  for (size in 1:length(ssizes)){
    samplVar <- array(data = NA, dim = c(npops,4,repl), 
                      dimnames = list(popNames(gind),c("Hmle","Hz","Hcs","Hchao"),c()))
    x <- for (rrr in 1:repl){
      resampled <- c()
      for (popul in popNames(gind)) resampled <- 
          append(resampled,sample(which(gind$pop==popul),ssizes[size],replace = FALSE))
      samplVar[,,rrr]<- apply(PopEstimates(gind[resampled]),2:3,mean)
    }
    results[,,size] <- apply(samplVar,1:2,sd)
  }
  return(results)
}

############ COALESCENT SIMULATIONS IN FASTSIMCOAL ############################
popSize <- 10000
# mutation rates in loci
mutRates  <- c(0.0001, 0.0002, 0.0005, 0.001)
# max no. of alleles in loci
maxRange <- c(3, 6, 9, 12, 15, 20)
numLoci <- length(mutRates) * length(maxRange)
# bottleneck sizes
botSize <- c(500, 50, 20)
botTime <- 20 # since t0
botLength <- 20 # since the beginning of the bottleneck
splitTime <- 50 # since t0
numPops <- length(botSize) + 1

# Setting the parameters for fastsimcoal
# Parameters of microsatellite loci
locus.params <-
  strataG::fscLocusParams(
    locus.type = "msat",
    num.loci = 1,
    mut.rate = rep(mutRates, each = length(maxRange)),
    #proportion of non-stepwise mutations
    gsm.param = 0.2,
    #maximum number of alleles at locus
    range.constraint = rep(maxRange, length(mutRates)),
    #diploid individuals
    ploidy = 2,
    #24 chromosomes (free recombination)
    chromosome = c(1:(numLoci))
  )
# Parameters of coalescent analysis
hist.ev <- strataG::fscHistEv(
  #dates (in generations) of historical events
  num.gen = c(
    rep(botTime, numPops - 1),
    rep(botTime + botLength, numPops - 1),
    rep(splitTime, numPops - 1)
  ),
  source.deme = as.numeric(c(seq(1:(
    numPops - 1
  )),
  seq(1:(
    numPops - 1
  )),
  seq(1:(
    numPops - 1
  )))),
  sink.deme = c(seq(1:(numPops - 1)),
                seq(1:(numPops - 1)),
                rep(0, numPops - 1)),
  prop.migrants = c(rep(1, 3 * (numPops - 1))),
  #simulating bottlenecks and recovery to previous popsize
  new.sink.size = c(as.numeric(rep(botSize / popSize)),
                    as.numeric(rep(popSize / botSize)),
                    rep(1, numPops - 1))
)
splSize <- rep(popSize, numPops)
pop.info <- strataG::fscPopInfo(
  pop.size = rep(popSize, numPops),
  sample.size = splSize,
  sample.times = c(rep(0, numPops))
)

##### fastsimcoal run
n_iter <- 1000
nsamples <- 2000
GIndS <- list()
PopData <- list()
ParamResults <- list()
SamplResults <- list()
# Starting simulations
SimStart <- Sys.time()
for (n in 1:n_iter) {
  SimTime <- Sys.time()
  simulated <- strataG::fastsimcoal(
    pop.info,
    locus.params,
    mig.rates = NULL,
    hist.ev,
    num.cores = 12,
    exec = paste0("fsc26"),
    delete.files = TRUE,
    quiet = TRUE
  )
  GIndS <- gtypes2genind(simulated)
  ParamResults[[n]] <- PopEstimates(GIndS, He=TRUE)
  numbers <- c()
  for (popul in popNames((GIndS)))
    numbers <-
    append(numbers, sample(which(GIndS$pop==popul),nsamples,replace = FALSE))
  PopData[[n]] <- GIndS[numbers]
  cat("###############", "\n")
  cat("##### Simulation", n, "out of", n_iter, "============", "\n")
  cat("##### It took: ",
      round(difftime(Sys.time(), SimTime, units = "secs")),
      "seconds",
      "================",
      "\n")
  cat("##### The run will end in",
      round(difftime(Sys.time(), SimStart, units = "mins") * (n_iter - n) / n),
      "minutes ====",
      "\n")
  cat("###############", "\n", "\n", "\n")
}

save.image()

# Simulating sampling

############ SAMPLING ERROR OF METRICES ###################################################

cl <- makeCluster(15)
testSamples <- c(5, 10, 20, 40, 80, 120, 200)
clusterExport(cl, list = c("testSamples","varTable","Shannon.gen","PopEstimates"))
SamplResults <- parLapply(cl,PopData, function (i) varTable(i,testSamples))
SplErr <- parLapply(cl,PopData, SamplingError)
stopCluster(cl)

vars <- dimnames(SamplResults[[1]])
ErrNames <- dimnames(SplErr[[1]])

HcsMiss <- array(data = NA, dim = c(3,4), dimnames = list(ErrNames[[3]],ErrNames[[1]]))
for (size in ErrNames[[3]]) for (popul in ErrNames[[1]]) {
  HcsMiss[size,popul] <- length(which(is.na(sapply(SplErr,function(i) i[popul,3,size]))))
}

loadfonts(device = "postscript")

botSizeTxt <- c(expression(italic("P"[20])),expression(italic("P"[50])),
                expression(italic("P"[500])),expression(italic("P"[C])))
popSizeTxt <- c(expression(italic("Ns = 5")),expression(italic("Ns = 20")),
                expression(italic("Ns = 80")))


pdf("Figure_1.pdf",width = 10,height = 15)
par(mfcol=c(4,3),oma = c(5, 5, 3, 3),mar = c(4.5, 2, 2, 2))
counter=0
# Drawing the plots
for (size in dimnames(SplErr[[1]])[[3]]) {
  for (popul in dimnames(SplErr[[1]])[[1]]) {
    counter = counter + 1
    resamplSD <- sapply(SplErr,function(i) i[popul,,size])
    d1 <- density(resamplSD[1,])
    d2 <- density(resamplSD[2,])
    d3 <- density(resamplSD[3,],na.rm = TRUE)
    d4 <- density(resamplSD[4,])
    plot(d1, main = "",
         xlim = c(range(resamplSD,na.rm = TRUE)),
         ylim = c(range(c(d1$y,d2$y,d3$y,d4$y),na.rm = TRUE)),
         lwd = 2, lty=1, cex = 2, cex.axis = 1.5, bty = "n",
         xlab = NA, ylab = NA)
    lines(d2, lwd = 2,lty=2)
    lines(d4, lwd = 2, lty=4)
    if (d3$n > 950) lines(d3, lwd = 2, lty=3)
    if (counter==9)
    legend("topright", lty = c(1:4),lwd=2, cex = 1.5,bty = "n", inset = c(-0.03,-0.03),
           c(expression(italic("H"[MLE])),expression(italic("H"[Z])),
             expression(italic("H"[CS])),expression(italic("H"[Chao]))), title = NULL)
    mtext(paste(" ",LETTERS[counter]), side = 3, outer = FALSE, line = -0.1, adj = 0, cex = 1.5)
    if ((counter+3) %% 4==0)
      mtext(popSizeTxt[(counter+3)%/%4], 
            side = 3, outer = FALSE, line = 1, cex = 1.2)
    if (counter > 8)
      mtext(botSizeTxt[counter-8], 
            side = 4, outer = FALSE, line = 1, cex = 1.2)
  }}
mtext('Standard deviation', side = 1, outer = TRUE, line = 2, cex = 1.5)
mtext('Density', side = 2, outer = TRUE, line = 2, cex = 1.5)
dev.off()


# Randomization test
niter=10000000
sdDif<- data.frame(test = integer(0), metrics = integer(0), "p-value" = integer(0))
for (size in ErrNames[[3]]) for (popul in ErrNames[[1]]) {
  j=0
  for (j in 1:3) {
    for (k in (j+1):4){
      test <- paste0("Ns ",size,", P: ",popul)
      metrics <- paste0(ErrNames[[2]][j]," vs ",ErrNames[[2]][k])
      a <- sapply(SplErr,function(i) i[popul,j,size])
      b <- sapply(SplErr,function(i) i[popul,k,size])
      if (max(a, na.rm = TRUE)<min(b,na.rm = TRUE)){
        n <- data.frame(test = test, metrics = metrics, "p-value" = 1/niter)
      }else{
        sa <- sample(a,niter,replace = TRUE)
        sb <- sample(b,niter,replace = TRUE)
        c <- sum(sa>sb,na.rm = TRUE)/(niter-sum(is.na(sa>sb)))
        if (c<=0.5){
          result <- c
        } else (result<- -(1-c))
        n <- data.frame(test = test, metrics  = metrics, "p-value" = round(result,6))
      }
      write.table(n,col.names = FALSE,quote = FALSE, sep = "\t")
      sdDif <- rbind(sdDif, n)
    }
    j=j+1
  }
}
citation("vegan")

MeanSDs <- data.frame(popul = integer(0), metric = integer(0), mean = integer(0))
for (popul in ErrNames[[1]]) for (metric in ErrNames[[2]]){
  MeanSDs <- rbind(MeanSDs,data.frame(popul = popul,metric = metric,
                  mean = round(mean(sapply(SplErr,function(i) i[popul,metric,1]),na.rm = TRUE),4)))
}
write.table(MeanSDs,file="SDmean.txt",quote = FALSE, sep = "\t")


############ SUMMARIZING SIMULATIONS ####
meanLociPop <-
  lapply(ParamResults, function (i)
    apply(i, c(2, 3), mean))

# Table of the estimators' mean values within samples and theirs SD
Distributions <-
  array(
    data = NA,
    dim = c((length(vars[[3]])), length(vars[[2]]), 4),
    dimnames = list(c(unlist(vars[[3]])), unlist(vars[[2]]), c("median", "min", "max","SD"))
  )

for (metric in 1:4){
  for (pop in 1:4) {
      data <- sapply(meanLociPop, function (i) i[pop, metric])
      Distributions[metric, pop, 1] <- round(median(data), 4)
      Distributions[metric, pop, 2] <- round(min(data), 4)
      Distributions[metric, pop, 3] <- round(max(data), 4)
      Distributions[metric, pop, 4] <- round(sd(data), 5)
      
}
}
write.table(file = "Sampling_Distributions.txt",
            sep = "\t",
            t(data.frame(Distributions)),
            quote = FALSE)


# Table of median and 95% confidence intervals of relative bias
# of Shannon H estimators averaged over loci
meanLociSample <-
  lapply(SamplResults, function (i)
    apply(i, c(2, 3, 4), mean, na.rm = TRUE))

rBias <-
  array(
    data = NA,
    dim = c(4, 7, 4),
    dimnames = list(unlist(vars[[3]]), unlist(vars[[4]]), c("median", "lowCI", "upCI", "SD"))
  )
for (metric in 1:length(vars[[3]])){
  for (size in 1:length(vars[[4]])) {
    data <- c()
    for (pop in 1:length(vars[[2]])) {
      data <-
        append(data,
               sapply(meanLociSample, function (i)
                 i[pop, metric, size]) - sapply(meanLociPop, function (i)
                   i[pop, metric]))
    }
    rBias[metric, size, 1] <- round(median(data), 4)
    rBias[metric, size, 2] <- round(Rfast::nth(data, 201, descending = FALSE), 4)
    rBias[metric, size, 3] <- round(Rfast::nth(data, 201, descending = TRUE), 4)
    rBias[metric, size, 4] <- round(sd(data), 5)
  }}
write.table(file = "rB_median_CIs.txt",
            sep = "\t",
            t(data.frame(rBias)),
            quote = FALSE)


############ Estimating bias and error of means over loci ####
bplotStd <-
  array(
    data = NA,
    dim = c(
      length(vars[[2]]),
      length(vars[[3]]),
      length(vars[[4]]),
      length(SamplResults)
    ),
    dimnames = list(vars[[2]], vars[[3]], vars[[4]], 1:1000)
  )
bpop <-
  array(
    data = NA,
    dim = c(length(vars[[2]]), length(vars[[3]])+1, length(SamplResults)),
    dimnames = list(vars[[2]], c(vars[[3]],"He"), 1:1000)
  )

for (pop in 1:length(vars[[2]])) {
  for (metric in 1:(length(vars[[3]])+1)) {
    bpop[pop, metric, ] <-
      sapply(ParamResults, function(i) mean(i[, pop, metric]))
  }
}

bplot <-
  array(data = NA,dim = c(length(vars[[2]]),length(vars[[3]]),length(vars[[4]]),
      length(SamplResults)),dimnames = list(vars[[2]], vars[[3]], vars[[4]], 1:1000)
  )

for (pop in 1:length(vars[[2]])) {
  for (metric in 1:length(vars[[3]])) {
    for (size in 1:length(vars[[4]])) {
      bplot[pop, metric, size, ] <- sapply(SamplResults, function(i) mean(i[, pop, metric, size]))
      bplotStd[pop, metric, size, ] <- (((bplot[pop, metric, size, ] - bpop[pop, metric, ])) / bpop[pop, metric, ])
    }
  }
}
loadfonts("postscript")

##### Plotting standardized bias for all sample sizes and metrics in simulated populations ###

Colours = c("#FFFFFF", "#DDDDDD", "#999999", "#666666")
Hnames = c(expression(italic("H"[MLE])),expression(italic("H"[Z])),
           expression(italic("H"[CS])),expression(italic("H"[Chao])))
pdf(file= "Figure_2.pdf", width = 10, height = 7,
             family = "sans", fonts = c("Arial"))
par(mfrow = c(2, 2),  oma = c(4, 4, 0, 0),mar = c(3, 1, 2, 1))
counter = 0
for (metric in c(1:length(vars[[3]]))) {
  counter = counter + 1
  boxplot(
    t(bplotStd[pop, metric, , ]),
    boxfill = NA, border = NA, frame = F,cex.axis = 1.2,
    las = 1,  yaxp = c(-0.3, 0.2, n = 5),bty = "n",
    ylim = c(range(bplotStd[, metric, , ], na.rm = TRUE))
  )
  abline(h = 0, lty = 2, lwd = 1)
  legend("bottomright",legend =Hnames[counter],bty = "n",cex = 1.2,adj = c(0.2,0.2))
  mtext(LETTERS[counter], side = 3, outer = FALSE, line = 0.5, adj = 0, cex = 1.5)
  for (pop in 1:length(vars[[2]])) {
    boxplot(
      t(bplotStd[pop, metric, , ]),
      add = TRUE, frame = F,
      boxfill = Colours[pop],
      xaxt = "n", yaxt = "n", boxwex = 0.15,
      range = 5,outline = FALSE,pch = 19, bty = "n",
      at = (1:ncol(t(bplotStd[pop, metric, , ])) + (pop / 5 - 0.5))
    )
  }
}
mtext('Number of samples', side = 1, outer = TRUE, line = 2, cex = 1.5)+
mtext(expression(italic("rB")), side = 2, outer = TRUE, line = 2, cex = 1.5)
dev.off()

############ GLM Analyses of populations' means ####
GLMpopError <- data.frame()
n = 0
for (pop in 1:length(vars[[2]])) {
  for (metric in 1:length(vars[[3]])) {
    for (size in 1:length(vars[[4]])) {
      range = c((n * 1000 + 1):(n * 1000 + 1000))
      n = n + 1
      GLMpopError[range, 1] <- abs(bplotStd[pop, metric, size, ])
      GLMpopError[range, 2] <- rep(vars[[3]][metric], 1000)
      GLMpopError[range, 3] <- rep(vars[[2]][pop], 1000)
      GLMpopError[range, 4] <- rep(vars[[4]][size], 1000)
      GLMpopError[range, 5] <- bpop[pop, 5, ]
    }
  }
}
names(GLMpopError) <- c("error", "metric", "pop", "size", "He")
GLMpopError$metric <- as.factor(GLMpopError$metric)
GLMpopError$pop <- as.factor(GLMpopError$pop)
GLMpopError$size <- as.integer(GLMpopError$size)
GLMpopError$size <- as.factor(GLMpopError$size)
GLMpopError200 <- GLMpopError
GLMpopError200$size <- relevel(GLMpopError200$size,ref=7)
which(is.nan(GLMpopError$error))

# GLM.gauss <- glm(error ~ metric * He * size,
#                  data = GLMpopError,
#                  family = gaussian(link = identity))
# GLM.gauss.log <- glm(error ~ metric * He * size,
#                      data = GLMpopError,
#                      family = gaussian(link = log))
# GLM.gauss.inverse <-  glm(error ~ metric * He * size,
#       data = GLMpopError, family = gaussian(link = inverse))
# GLM.gamma <- glm(error ~ metric * He * size, data = GLMpopError,
#                  family = Gamma())
GLM.gamma.log <- glm(error ~ metric * He * size, data = GLMpopError,
                 family = Gamma(link = log))

# MuMIn::model.sel(GLM.gamma,
#                  GLM.gauss,
#                  GLM.gauss.inverse,
#                  GLM.gauss.log,
#                  GLM.gamma.log)
# 

############ Printing Size and He Effect #####
linMod <-
  glm(error ~ metric * size,
      data = GLMpopError,
      family = Gamma(link = log))
car::Anova(linMod)
summary(linMod)
capture.output(summary(linMod),file = "TableS1.txt")

pdf(paste0("Figure_3.pdf"), width = 5, height = 5)
plot(
  effects::predictorEffect("size", linMod),
  axes = list(x = list(rug = FALSE),
              y = list(ticks=list(at=c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.15)))),
  lines = list(
    multiline = TRUE,
    lty = c(1:5),
    lwd = 2
  ),
  symbols = list(pch = c(15, 25, 17, 18)),
  colors = "black",
  main = NULL,
  ylab = list(expression(italic("rE"))),
  xlab = "Sample size",
  #  rescale.axis = FALSE,
  lattice = list(
    key.args = list(
      title = NULL,
      space = "right",
      columns = 1,
      border = FALSE,
      cex = 1,
      cex.title = 1.15,
      x=0.56, y=0.9, 
      padding.text=1,
      text = list(c(expression(italic(H[Chao])),expression(italic(H[CS])),
                    expression(italic(H[MLE])),expression(italic(H[Z]))))
    )
  )
)
dev.off()

linMod2 <-
  glm(error ~ He*metric*size,
      data = GLMpopError,
      family = Gamma(link = log))
pdf("Figure_4a.pdf", width =10, height = 6, family  =  "sans", fonts = "Arial")
plot(
  effects::predictorEffect("He", linMod2),
  axes = list(x = list(rug = FALSE), alternating = FALSE),
  lines = list(
    multiline = TRUE,
    lty = c(1:6, 2),
    lwd = c(rep(2, 6), 1.5)
  ),
  between = list(x = 0.5,y=0.5),
  colors = c(rep("#000000", 6), rep("#888888", 1)),
  main = NULL,
  ylab = list(expression(italic("rE"))),
  xlab = list(expression(italic("He"))),
  z.var = "size",
  lattice = list(
    strip = list(
      factor.names = FALSE),
    key.args = list(
      space = "right",
      title = "Sample size",
      columns = 1,
      border = FALSE,
      cex = 1,
      cex.title = 1.15
    )
  )
)
dev.off()



############ Tuckey HSD method #####
cld_res_pop <-
  array(dim = c(0, length(vars[[4]])))
for (metric in 1:length(vars[[3]])) {
  testData <- GLMpopError[GLMpopError$metric == vars[[3]][metric], ]
  mod <- glm(error ~ size, data = testData, family = Gamma(link = log))
  glht_metric <- glht(mod, linfct = mcp(size = "Tukey"))
  cldres <- cld(glht_metric)
  cld_res_pop <- rbind(cld_res_pop, cldres$mcletters$Letters)
  cld_res_pop <- rbind(cld_res_pop, round(glht_metric$coef, 5))
}

write.table(
  file = "Pop_TuckeySize.txt",
  cld_res_pop,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

cld_size <-
  array(dim = c(0, length(vars[[3]])))
for (size in vars[[4]]) {
  testData <- GLMpopError[GLMpopError$size == size, ]
  mod <- glm(error ~ metric, data = testData, family = Gamma(link = log))
  glht_size <-
    glht(mod,
         linfct = mcp(metric = "Tukey"))
  par(mai = c(1,2,1,1))
  print(plot(glht_size, main = paste(size)))
  cldres <- cld(glht_size)
  cld_size <- rbind(cld_size, size)
  cld_size <- rbind(cld_size, cldres$mcletters$Letters)
  cld_size <- rbind(cld_size, round(glth_size$coef, 5))
};cld_size
write.table(
  file = "Pop_TuckeyMetric.txt",
  cld_size,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)



########## LOCUS EFFECTS #####
bplotStdLoci <-
  array(
    data = NA,
    dim = c(
      length(vars[[1]]),
      length(vars[[2]]),
      length(vars[[3]]),
      length(vars[[4]]),
      length(SamplResults)
    ),
    dimnames = list(vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1:1000)
  )
bpopLoci <-
  array(
    data = NA,
    dim = c(
      length(vars[[1]]),
      length(vars[[2]]),
      (1 + length(vars[[3]])),
      length(ParamResults)
    ),
    dimnames = list(vars[[1]], vars[[2]], c(unlist(vars[[3]]), "He"), 1:1000)
  )
bplotLoci <-
  array(
    data = NA,
    dim = c(
      length(vars[[1]]),
      length(vars[[2]]),
      length(vars[[3]]),
      length(vars[[4]]),
      length(SamplResults)
    ),
    dimnames = list(vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1:1000)
  )
for (locus in 1:length(vars[[1]])) {
  for (pop in 1:length(vars[[2]])) {
    for (metric in 1:(length(vars[[3]]) + 1)) {
      bpopLoci[locus, pop, metric, ] <-
        sapply(ParamResults, function(i)
          i[locus, pop, metric])
    }
  }
}
for (locus in 1:length(vars[[1]])) {
  for (pop in 1:length(vars[[2]])) {
    for (metric in 1:(length(vars[[3]]))) {
      for (size in 1:length(vars[[4]])) {
        bplotLoci[locus, pop, metric, size, ] <-
          sapply(SamplResults, function(i)
            i[locus, pop, metric, size])
        bplotStdLoci[locus, pop, metric, size, ] <-
          (((bplotLoci[locus, pop, metric, size, ] - bpopLoci[locus, pop, metric, ])) /
             bpopLoci[locus, pop, metric, ])
      }
    }
  }
}

##### Detecting 'NaNs'######
NaNsloci <-
  data.frame(
    locus = integer(0),
    pop = character(0),
    metric = character(0),
    iteration = integer(0),
    size = integer(0)
  )
for (locus in 1:length(vars[[1]]))
  for (pop in 1:length(vars[[2]]))
    for (metric in 1:length(vars[[3]]))
      for (size in 1:length(vars[[4]]))
        for (iter in 1:1000) {
          if (is.nan(bplotLoci[locus, pop, metric, size, iter])) {
            NaNsloci <-
              rbind(
                NaNsloci,
                data.frame(
                  locus = vars[[1]][locus],
                  pop = vars[[2]][pop],
                  metric = vars[[3]][metric],
                  size = vars[[4]][size],
                  iteration = (1:1000)[iter]
                )
              )
          }
        }
tableNaN <- matrix(ncol = 2, nrow = 24)
colnames(tableNaN) <- vars[[2]][1:2]
rownames(tableNaN) <- vars[[1]]
for (locus in 1:24) {
  for (pop in 1:2) {
    tableNaN[locus, pop] <-
      sum((NaNsloci$locus == vars[[1]][locus]) * (NaNsloci$pop == vars[[2]][pop]))
  }
}
sum(colSums(tableNaN))

write.table(
  file = "NaNs.txt",
  cbind(tableNaN, locus.params$param.6, locus.params$param.4),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

save.image()
##### Preparing table for GLM for loci ####
GLMlocError <- data.frame()
n = 0
for (locus in 1:length(vars[[1]])) {
  for (pop in 1:length(vars[[2]])) {
    for (metric in 1:length(vars[[3]])) {
      for (size in 1:length(vars[[4]])) {
        range = c((n * 1000 + 1):(n * 1000 + 1000))
        n = n + 1
        GLMlocError[range, 1] <-
          abs(bplotStdLoci[locus, pop, metric, size, ])
        GLMlocError[range, 2] <- rep(vars[[1]][locus], 1000)
        GLMlocError[range, 3] <- rep(vars[[2]][pop], 1000)
        GLMlocError[range, 4] <- rep(vars[[3]][metric], 1000)
        GLMlocError[range, 5] <- rep(vars[[4]][size], 1000)
        GLMlocError[range, 6] <- bpopLoci[locus, pop, 5, ]
      }
    }
  }
}

names(GLMlocError) <-
  c("error", "locus", "pop", "metric", "size", "He")
GLMlocError$locus <- as.factor(GLMlocError$locus)
GLMlocError$pop <- as.factor(GLMlocError$pop)
GLMlocError$metric <- as.factor(GLMlocError$metric)
GLMlocError$size <- as.integer(GLMlocError$size)
GLMlocError$size <- as.factor(GLMlocError$size)

# Checking for errors = 0 (they are very unlikely,
# but cannot be included in gamma model)
Error0 <- which(GLMlocError$error == 0)
if (length(Error0) != 0) {
  write.table(file = "0_errors.txt",
              GLMlocError[Error0, ],
              row.names = FALSE,
              quote = FALSE)
  GLMlocError <- GLMlocError[GLMlocError$error != 0, ]
}

# Model choice
GLM.Loc.gauss <- glm(error ~ metric + locus,
                     data = GLMlocError,
                     family = gaussian(link = identity))
GLM.Loc.gauss.log <- glm(error ~ metric + locus,
                         data = GLMlocError,
                         family = gaussian(link = log))
GLM.Loc.gauss.inverse <-
  glm(error ~ metric + locus,
      data = GLMlocError,
      family = gaussian(link = inverse))
GLM.Loc.gamma <- glm(error ~ metric + locus, data = GLMlocError,
                     family = Gamma())
GLM.Loc.gamma.log <- glm(error ~ locus * metric,
                         data = GLMlocError,
                         family = Gamma(link = log))
MuMIn::model.sel(
  GLM.Loc.gauss,
  GLM.Loc.gauss.log,
  GLM.Loc.gauss.inverse,
  GLM.Loc.gamma,
  GLM.Loc.gamma.log
)


# Preparing Figure S1 - locus and size specific error level association with
# locus expected heterozygosity predictor for all metrics
# [sizes need to selected to due to computer power limitations]
sizes = c(1, 3, 5)
for (size in sizes) {
  testData <- GLMlocError[GLMlocError$size == vars[[4]][size], ]
  mod <-
    glm(error ~ locus * metric * He,
        data = testData,
        family = Gamma(link = log))
  # out <- capture.output(summary(mod))
  # out2 <- capture.output(car::Anova(mod, type = 3))
  # cat(
  #   out,
  #   file = paste0("Loc_", vars[[4]][size], "_summary.txt"),
  #   sep = "\n",
  #   append = FALSE
  # )
  # cat(
  #   out2,
  #   file = paste0("Loc_", vars[[4]][size], "_summary.txt"),
  #   sep = "\n",
  #   append = TRUE
  # )
  pdf(file = paste0("FigureS", which(sizes==size), ".pdf"), width = 12, height = 12)
  print(
    plot(
      effects::predictorEffect("He", mod),
      axes = list(x = list(rug = FALSE, rotate = 90,He=list(lab = list(expression(italic("He"))),
                                                            ticks=list(at=c(0.1,.5,0.9)))),
                  y = list(ticks=list(at=c(0.01,0.03,0.1,0.3,1,3,10)),lab = "ASE"),
                  alternating = FALSE),
      main = NULL,
      lines = list(multiline = TRUE,lty = c(2:5),lwd = c(1.5, 1.3, 1.1, 0.9)),
      colors = c("#222222"),
      las = 2,
      ylab = list(expression(italic("rE"))),
      lattice = list(
        key.args = list(
          title = NULL,
          space = "top",
          columns = 4,
          border = FALSE,
          cex = 1,
          padding.text=1,
          text = list(c(expression(italic(H[Chao])),expression(italic(H[CS])),
                        expression(italic(H[MLE])),expression(italic(H[Z]))))),
          strip = list(
            factor.names = FALSE),
          layout = c(6, 4))
      )[c(19:24, 13:18, 7:12, 1:6)]
    )
    dev.off()
}

for (size in sizes) {
  for (metric in vars[[3]]){
  testData <- GLMlocError[GLMlocError$size == vars[[4]][size], ]
  mod <-
    glm(error ~ locus*He,
        data = testData,
        family = Gamma(link = log))
  out <- capture.output(summary(mod))
  cat(
    out,
    file = paste0("Size_", vars[[4]][size],"_m_",metric,"_summary.txt"),
    sep = "\n",
    append = FALSE
  )
  out2 <- capture.output(car::Anova(mod, type = 3))
  cat(
    out2,
    file = paste0("Size_", vars[[4]][size],"_m_",metric, "_Anova.txt"),
    sep = "\n",
    append = TRUE
  )
}
}



# Preparing Figure 5 - locus specific error level for all metrics
testData <- GLMlocError[GLMlocError$size %in% vars[[4]][sizes], ]
modFig4 <- glm(error ~ metric * locus * size,
            data = testData,
            family = Gamma(link = log))
write.table(
  file = "Locus_size_metric_GLM.txt",
  capture.output(summary(modFig4)),
  quote = FALSE,
  row.names = FALSE
)
write.table(
  file = "Locus_size_metric_GLM_AOV.txt",
  capture.output(car::Anova(modFig4)),
  quote = FALSE,
  row.names = FALSE
)


cbbPalette <- c("black", "red", "gold", "blue")

##################################
threeSizes <- GLMlocError[GLMlocError$size %in% vars[[4]][c(1,3,5)], ]
mod3sizes <- glm(error ~ locus:metric + size:metric ,
                 data = threeSizes,
                 family = Gamma(link = log))
effect3sizes <- effects::predictorEffect(c("metric"), mod3sizes)

postscript("Figure_5.eps", fonts=c("serif", "Arial"))
plot (effect3sizes,
    symbols = list(pch = c(4,5,6,15,17,18), cex = 0.5),
    lines = list(
      multiline = TRUE,
      lty = c(1:6),
      lwd = 2,
      col = rep(cbbPalette,each = 6)
    ),
    axes = list(
      alternating = FALSE, 
      between = list(x = 0.5,y=1),
      
      grid = TRUE,
      x = list(rug = FALSE, 
               rotate = 90),
      y = list(rescale.axis = TRUE,
               type = "rescale")),
    ylab = list(expression(italic("rE"))),
    xlab = "Metric",
    z.var = "metric",
    main = NULL,
    text = list(fontfamily = "serif"),
    lattice = list(
      layout = c(3,1),
      key.args = list(
        fontfamily = "serif",
        space = "top",
        columns = 4,
        border = 0,
        cex = 1,
        cex.title = 1.15
      )
    )
)
dev.off()


# Preparing CLD table for Tukey HSD test
cld_res <-
  array(dim = c(0, length(vars[[1]])), dimnames = list(c(), unlist(vars[[1]])))
for (metric in 1:length(vars[[3]])) {
  for (size in sizes) {
    testData <- GLMlocError[GLMlocError$size == vars[[4]][size], ]
    testData <- testData[testData$metric == vars[[3]][metric], ]
    mod <- glm(error ~ locus,
               data = testData,
               family = Gamma(link = log))
    glht_locus <-  glht(mod, linfct = mcp(locus = "Tukey"))
    cldres <- cld(glht_locus)
    cld_res <- rbind(cld_res, cldres$mcletters$Letters)
    cld_res <- rbind(cld_res, round(glht_locus$coef, 6))
  }
}
write.table(
  file = "Loci_Tuckey.txt",
  cld_res,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
# save.image()