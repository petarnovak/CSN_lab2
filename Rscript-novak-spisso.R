# 1: Introduction
library(xtable)
write_table <- function(label,file,df) {
   degree_sequence = read.table(file, header = FALSE)
   df <- rbind(df,data.frame(label, length(degree_sequence$V1), max(degree_sequence$V1),
                                       sum(degree_sequence$V1)/length(degree_sequence$V1),
                                       length(degree_sequence$V1)/sum(degree_sequence$V1)))
   return(df)
   }

source = read.table("list_in.txt", 
         header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
         as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
        )

lang.df <- data.frame()

for (x in 1:nrow(source)) {
  lang.df <- write_table(source$language[x], source$file[x],lang.df)
}
colnames(lang.df) <- c("Language", "N", "Maximum degree", "M/N", "N/M")
#Table with all languages and some important values
lang.df
print(xtable(lang.df), file = "Table 1.tex")

# 2: Visualization
# Function to plot the degree sequence of a given language
degree_plot <- function(language,file){
  degree_sequence = read.table(file, header = FALSE)
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum, main = language, xlab = "degree", ylab = "Number of vertices", log = "y")
}

par(ask = TRUE) #We stop at each plot
for (x in 1:nrow(source)){
  degree_plot(source$language[x], source$file[x])
}

# 3: Log-likelihood function
require(stats4) # for MLE
require(VGAM) # for the Riemann-zeta function

# Geometric distribution
minus_log_like_geo <- function(p){
  -(sum(x)-length(x)) * log(1-p) - length(x) * log(p)
}

# Poisson distribution
minus_log_like_pois <- function(lambda){
  C <- 0
  for (i in 1:length(x)) {
    C = C + sum(log(2:x[i]))
  }
  - sum(x) * log(lambda) + length(x) * (lambda + log(1-exp(1)^(-lambda))) + C
}

# Zeta distribution
minus_log_like_zeta <- function(gamma){
  length(x) * log(zeta(gamma)) + gamma * sum(log(x))
}

# Zeta (gamma=2) distribution
minus_log_like_zeta2 <- function(){
  mle <- length(x) * log(pi^2/6) + 2 * sum(log(x))
  return(mle)
}

# Right-truncated zeta distribution
minus_log_like_zeta_trunc <- function(gamma){
  length(x) * log(sum((1:max(x))^(-gamma))) + gamma * sum(log(x)) 
}

# Function to calculate the mle useful parameters
mle_calc <- function(p,lambda,x){
  mle_geo <- mle(minus_log_like_geo,
                 start = list(p = p),
                 method = "L-BFGS-B",
                 lower = c(0.001),
                 upper = c(0.999))
  geo_par <- c(attributes(summary(mle_geo))$coef[1],attributes(summary(mle_geo))$m2logL)
  mle_pois <- mle(minus_log_like_pois,
                  start = list(lambda = lambda),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  pois_par <- c(attributes(summary(mle_pois))$coef[1],attributes(summary(mle_pois))$m2logL)
  mle_zeta <- mle(minus_log_like_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  zeta_par <- c(attributes(summary(mle_zeta))$coef[1],attributes(summary(mle_zeta))$m2logL)
  mle_zeta_trunc <- mle(minus_log_like_zeta_trunc,
                        start = list(gamma = 2),
                        method = "L-BFGS-B",
                        lower = c(1.0000001))
  zeta_trunc_par <- c(attributes(summary(mle_zeta_trunc))$coef[1],
                      attributes(summary(mle_zeta_trunc))$m2logL)
  df <- data.frame(pois_par,geo_par,zeta_par,zeta_trunc_par)
  colnames(df) <- NULL
  return(df)
}

# Functions to calculate the AIC
get_AIC <- function(m2logL,K,N){
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

AIC_calc <- function(label, m2logL.list,df,x){
  AIC.list <- sapply(m2logL.list,get_AIC,K=1,N=length(x))
  mle_zeta2 <- 2*minus_log_like_zeta2()
  AIC.list <- append(AIC.list,get_AIC(mle_zeta2,0,length(x)),after = 3)
  best.AIC <- min(AIC.list)
  AIC.list <- AIC.list-best.AIC
  df <- rbind(df,data.frame(label,t(AIC.list)))
  return(list(AICdf = df, AICbest = best.AIC))
}

param.df <- data.frame() #table with best parameters
AIC.df <- data.frame() #AIC table
bestAIC.list <- list()

for (i in 1:nrow(source)){
  x <- read.table(source$file[i], header = FALSE)$V1
  param.list <- mle_calc(lang.df$"N/M"[i],lang.df$"M/N"[i],x)
  param.df <- rbind(param.df,data.frame(source$language[i], param.list[1,],max(x)))
  AIC <- AIC_calc(source$language[i],param.list[2,],AIC.df,x)
  AIC.df <- AIC$AICdf
  bestAIC.list[i] <- AIC$AICbest
}

colnames(param.df) <- c("Language", "lambda", "p", "gamma 1","gamma 2","k max")
colnames(AIC.df) <- c("Language", "1", "2", "3", "4", "5")

param.df
print(xtable(param.df), file = "Table Parametri.tex")
AIC.df
print(xtable(AIC.df), file = "Table AIC.tex")

# Plots of distributions vs real data
# Probability functions
#Displaced geometric
geo_dist <- function(p,k){
  q <- (1-p)^(k-1)*p
  return(q)
}

#Displaced poisson
pois_dist <- function(lambda,k){
  p <- (lambda^k) * exp(-lambda) / (factorial(k) * (1 - exp(-lambda)))
}

#Zeta truncated
zetatrunc_dist <- function(k_max,gamma,k){
  k^(-gamma)/(sum((1:k_max)^(-gamma)))
}

full_plot <- function(i,label,file,df){
  degree_sequence = read.table(file, header = FALSE)
  degree_spectrum = table(degree_sequence) 
  barplot(degree_spectrum, main = label, xlab = "degree", 
          ylab = "Number of vertices", log = "y")
  #degree_spectrum = data.frame(table(degree_sequence))
  #plot(degree_spectrum$V1,degree_spectrum$Freq,main = label,type="l",
  #     xlab = "degree", ylab = "Number of vertices",log="y")
  
  x <- 1:max(degree_sequence)
  
  geo_prob <- sapply(x, geo_dist,p=df$p[i])
  lines(x,geo_prob*nrow(degree_sequence),type="l",col = "blue",lwd = 3)
  
  pois_prob <- sapply(x, pois_dist,lambda=df$lambda[i])
  lines(x,pois_prob*nrow(degree_sequence),type="l",col = "green",lwd = 3)
  
  zeta_prob <- sapply(x, dzeta, shape = df$`gamma 1`[i])
  lines(x,zeta_prob*nrow(degree_sequence),type="l",col = "orange",lwd = 3)
  
  zeta2_prob <- sapply(x, dzeta, shape = 2)
  lines(x,zeta2_prob*nrow(degree_sequence),type="l",col = "magenta",lwd = 3)
  
  zetatrunc_prob <- sapply(x, zetatrunc_dist, gamma = df$`gamma 2`[i], 
                           k_max = df$`k max`[i])
  lines(x,zetatrunc_prob*nrow(degree_sequence),type="l",col = "red",lwd = 3)
  legend("topright", legend = c("Geometric", "Poisson", "Zeta", "Zeta (Lambda=2)", "Zeta truncated"), 
         col = c("blue","green","red","magenta","orange"), lty = 1, lwd = 2)
}

for (x in 1:nrow(source)){
  full_plot(x,source$language[x], source$file[x],param.df)
}


# 6: Samples from discrete distribution
folder_path <- "./samples_from_discrete_distributions/data"
files <- list.files(path = folder_path, full.names = TRUE)
prob_list <- c("geo 0.05","geo 0.1","geo 0.2","geo 0.4","geo 0.8","zeta 1.5","zeta 2.5",
          "zeta 2","zeta 3.5","zeta 3")

# Introduction
prob.df <- data.frame()

for (x in 1:length(files)) {
  prob.df <- write_table(prob_list[x], files[x],prob.df)
}

colnames(prob.df) <- c("Distribution", "N", "Maximum degree", "M/N", "N/M")
prob.df

# Visualization
for (i in 1:length(files)){
  degree_sequence = read.table(files[i], header = FALSE)
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum, main = paste("Distribution :",prob_list[i]), xlab = "degree", ylab = "number of vertices", log = "y")
}

param.sample.df <- data.frame() #table with best parameters
AIC.sample.df <- data.frame() #AIC table

for (i in 1:length(files)){
  x <- read.table(files[i], header = FALSE)$V1
  param.list <- mle_calc(prob.df$"N/M"[i],prob.df$"M/N"[i],x)
  param.sample.df <- rbind(param.sample.df,data.frame(prob.df$Distribution[i], param.list[1,],max(x)))
  AIC.sample.df <- AIC_calc(prob.df$Distribution[i],param.list[2,],AIC.sample.df,x)
}

colnames(param.sample.df) <- c("Distribution", "lambda", "p", "gamma 1","gamma 2","k max")
colnames(AIC.sample.df) <- c("Distribution", "1", "2", "3", "4","5")

param.sample.df
AIC.sample.df

# Plots of distributions vs sample data
# NOTE: the code requires a lot of time to plot gamma 1.5
for (x in 1:length(files)){
  full_plot(x,prob_list[x],files[x],param.sample.df)
}


# Additional work
# Altmann distribution
minus_log_like_altmann <- function(gamma,delta){
  gamma * sum(log(x)) + delta * sum(x) + length(x)*log(sum(((1:length(x))^(-gamma)*exp(-delta*(1:length(x))))))
}

altmann_dist <- function(gamma,delta,k,N){
  k^(-gamma) * exp(-delta * k) * (1/sum((1:N)^(-gamma) * exp(-delta * (1:N))))
}

altmann.param <- data.frame()

for (i in 1:nrow(source)){
  degree_sequence <- read.table(source$file[i], header = FALSE)
  x <- degree_sequence$V1
  mle_altmann <- mle(minus_log_like_altmann,
                     start = list(gamma = 2, delta = 0.01),
                     method = "L-BFGS-B",
                     lower = c(0.001,0.001))
  alt.gamma <- attributes(summary(mle_altmann))$coef[1]
  alt.delta <- attributes(summary(mle_altmann))$coef[2]
  altmann.AIC <- get_AIC(attributes(summary(mle_altmann))$m2logL,2,length(x)) - bestAIC.list[[i]]
  altmann.param <- rbind(altmann.param,data.frame(alt.gamma,alt.delta,altmann.AIC))
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum, main = source$language[i], xlab = "degree", 
          ylab = "Number of vertices", log = "y")
  
  z <- 1:max(degree_sequence)
  altmann_prob <- sapply(z, altmann_dist,gamma=alt.gamma,delta=alt.delta,N=max(x))
  lines(z,altmann_prob*nrow(degree_sequence),type="l",col = "red",lwd = 3)
}

colnames(altmann.param) <- c("gamma", "delta", "Alt.AIC-best.AIC")

altmann.param
print(xtable(altmann.param), file = "Table Altmann.tex")


plot(degree_spectrum$V1,degree_spectrum,main = "Geo 0.05", xlab = "degree", 
     ylab = "Number of vertices",log="y")

z <- 1:max(degree_sequence)

altmann_prob <- sapply(z, altmann_dist,gamma=0.5,delta=0.01,N=max(x))
lines(z,altmann_prob,type="l",col = "blue")
minus_log_like_altmann(0.5,0.01)