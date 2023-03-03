

library(readxl)
library(plyr)

# Select Population
Ethnicity <- "EUR"

getSEER_data <- function(Ethnicity="EUR") {
  ethnicity.select <- switch(Ethnicity,
                             "EUR"="NH White",
                             "AFR"="NH Black",
                             "EAS"="AIAN",
                             "LA"="Hispanic")
  
  #mortality 
  mort.SEER <- read.table(file="Data/Seer_PCa_asmr_byrace.csv", header=T, sep=",")
  mort.SEER <- mort.SEER[mort.SEER$Race.Ethnicity==ethnicity.select,]
  mort.SEER <- apply(mort.SEER, 1, FUN=function(v) {
    age.v <- unlist(strsplit(trimws(strsplit(v["Age"], "years")), "[-+]+"))
    age.v <- age.v[1]:age.v[length(age.v)]
    rate.v <- rep(as.numeric(v["Rate.per.100.000"])/100000, times=length(age.v))
    return(cbind(age.v, rate.v))
  })
  mort.SEER <- do.call(rbind, mort.SEER)
  seer_mu_c <- rep(0, max(mort.SEER[,1])+1)
  seer_mu_c[mort.SEER[,1]+1] <- mort.SEER[,2]
  
  #incidence
  inc.SEER <- read.table(file="Data/Seer_PCa_asir_byrace.csv", na.strings=c(" ", "^"), sep=",", header=T)
  inc.SEER <- inc.SEER[inc.SEER$stage_prostate=="All",]
  inc.SEER <- inc.SEER[inc.SEER$Race.Ethnicity==ethnicity.select,]
  inc.SEER <- apply(inc.SEER, 1, FUN=function(v) {
    age.v <- unlist(strsplit(trimws(strsplit(v["Age"], "years")), "[-+]+"))
    age.v <- age.v[1]:age.v[length(age.v)]
    rate.v <- as.numeric(v["Rate.per.100.000"])/100000
    #rate.v <- rep(ifelse(is.na(rate.v), 0, rate.v/100000), times=length(age.v))
    return(cbind(age.v, rate.v))
  })
  inc.SEER <- do.call(rbind, inc.SEER)
  seer_mu <- rep(0, max(inc.SEER[,1])+1)
  seer_mu[inc.SEER[,1]+1] <- inc.SEER[,2]
  return(list(mu=seer_mu, mu_c=seer_mu_c))
}

calc_AR <- function(mu=NULL, mu_c=NULL, betas.mu=NULL, PRS.categories=NULL){
  # setting up the data
  max.age <- length(mu)-1
  age_t <- 0:max.age # age in terms of time/years
  age_index <- age_t+1
  OR.cat <- exp(betas.mu)
  
  # risk-specific categories for incidence calculations 
  S_c <- exp(-cumsum(mu_c)) # probability of not dying from another cause of death by age t, based on the age- specific mortality rates mc_p.
  #f <- PRS.categories[2:length(PRS.categories)]- PRS.categories[1:(length(PRS.categories)-1)] # proportion or population in each PRS category
  f <- PRS.categories[2:length(PRS.categories)]- PRS.categories[1:(length(PRS.categories)-1)] # proportion or population in each PRS category
  K = length(f)
  lambda_0 <- rep(NA, max.age+1)
  S_k <- matrix(NA, nrow=(max.age+1), ncol=K) # age x risk k
  lambda_k <- matrix(NA, nrow=(max.age+1), ncol=K)  # age x risk k
  S_k[0+1,] <- 1 # survival is 1 at age 0
  lambda_k[0+1,] <- 0 # 0 incidence at age 0
  a <- 1 # age 1
  t <- a+1
  lambda_0[t] <- mu[t]* sum(f*S_k[(t-1),]) / sum(f*S_k[(t-1),]*OR.cat)  # initialize value at first time t=a+1
  lambda_k[t,] <- lambda_0[t]*OR.cat # initialize value at first time t=a+1
  for(k in 1:K) { # initialize for all categories
    S_k[t, k] <- exp(-sum(lambda_k[1:(t-1),k]))
  }
  for(a in 2:max(age_t)) {
    t <- a+1
    lambda_0[t] <- mu[t]* sum(f*S_k[(t-1),]) / sum(f*S_k[(t-1),]*OR.cat)
    lambda_k[t,] <- lambda_0[t]*OR.cat
    for(k in 1:K) {
      S_k[t, k] <- exp(-sum(lambda_k[1:(t-1),k]))
    }
  }
  
  # AR calculations
  AR <- matrix(NA, nrow=(max.age+1), ncol=K)
  a <- 0
  t <- a+1
  for(k in 1:K) {
    AR[t,k] <- lambda_k[t,k]*S_k[t,k]*S_c[t]
  }
  for(a in 1:max(age_t)) {
    t<- a+1
    for(k in 1:K) {
      AR[t,k] <- AR[(t-1),k]+lambda_k[t,k]*S_k[t,k]*S_c[t]
    }
  }
  return(AR)
}
 

# overall meta-analysis results used for AR calculations: 1 SD increase OR
OR <- switch(Ethnicity,
             "EUR"=2.32,
             "AFR"=2.04,
             "EAS"=2.15,
             "LA"=2.12)

# calculate an OR for each level of % risk
quant.n <- seq(from=.01, to=0.99, by=.01) 
v.n <- qnorm(quant.n)
risk.at.each.level <- v.n*log(OR)
OR.at.each.level <- exp(risk.at.each.level)

PRS.categories <- c(0, quant.n, 1)
PRS.categories <- PRS.categories[PRS.categories!=0.5]
PRS.labels <- paste(paste(100*PRS.categories[1:(length(PRS.categories)-1)], 100*PRS.categories[2:(length(PRS.categories))], sep="% - "), "%", sep="")
PRS.ref <- "49% - 51%"
PRS.labels.s <- PRS.labels[!PRS.labels==PRS.ref]
OR.cat <- OR.at.each.level 



d <- getSEER_data(Ethnicity=Ethnicity)
mu <- d$mu
mu_c <- d$mu_c
betas.mu <- log(OR.cat)

AR.r <- calc_AR(mu=mu, mu_c=mu_c, betas.mu=betas.mu, PRS.categories=PRS.categories)
save(AR.r, file=paste("AR", Ethnicity, "RData", sep="."))


pdf(file=paste("AR.figure.", Ethnicity, ".pdf", sep=""))
yrPal <- colorRampPalette(c('green','yellow', 'orange', 'red'))
plot.breaks <- seq(from=-0.1, to=plyr::round_any(max(AR.r), .1, ceiling), by=.05)
plot.col <- function(v) {
  yrPal(length(plot.breaks))[as.numeric(cut(v,breaks = plot.breaks))]
}

max.age <- length(mu)-1
age <- 50:max.age
age_t <- 0:max.age
index <- age+1 
K <- length(OR.cat)
k<- K
#plot.col <- rainbow(K, rev=T)
max.y <- .5
plot(age_t[index], AR.r[index,k], ylab="Absolute Risk", xlab="Age", 
     pch=16, col=plot.col(AR.r[index,k]), ylim=c(0,max.y), yaxt="n", main=paste("Absolute Risk for", Ethnicity, "with 1SD OR = ", OR))
axis(4,at=seq(from=0, to=max.y, by=.05))
for(k in (K-1):1) {
  points(age_t[index], AR.r[index,k], pch=16, col=plot.col(AR.r[index,k]))
}
PRS.lines <- c(.9, .8, .51, .1)
for(PRS.index in PRS.lines) {
  k <- match(PRS.index, round(PRS.categories,2))-1
  lines(age_t[index], AR.r[index,k], type="l", lwd=2, lty=match(PRS.index, PRS.lines), col=1)
}
legend("top", title="GRS Categories", legend=PRS.labels[match(PRS.lines, round(PRS.categories,2))-1], col=1, lty=1:length(PRS.lines), lwd=2)

AR.levels <- ifelse(plot.breaks<0, 0, round(plot.breaks,2))
AR.labels <- paste(AR.levels[1:(length(AR.levels)-1)], "< AR <", AR.levels[2:(length(AR.levels))])
AR.labels <- c(AR.labels[AR.labels!="0 < AR < 0"], paste("AR >", AR.levels[length(AR.levels)]))
legend("topleft", title="Absolute Risk Levels", legend=AR.labels, col=yrPal(length(plot.breaks))[as.numeric(cut(AR.levels,breaks = plot.breaks))], pch=16)

dev.off()

AR.r <- as.data.frame(AR.r)
names(AR.r) <- PRS.labels
row.names(AR.r) <- age_t
write.table(AR.r, file=paste("AR", Ethnicity, "csv", sep="."), sep=",", row.names=T)

write.table(cbind(quant.n, OR.at.each.level), file=paste("OR", Ethnicity, "csv", sep="."), sep=",", row.names=F)

