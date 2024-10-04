library(xtable)

# Function to create the initial summary table
write_table <- function(label,file,df) {
   degree_sequence = read.table(file, header = FALSE)
   df <- rbind(df,data.frame(label, length(degree_sequence$V1), max(degree_sequence$V1),
                                       sum(degree_sequence$V1)/length(degree_sequence$V1),
                                       length(degree_sequence$V1)/sum(degree_sequence$V1)))
   return(df)
   }

# LOG-LIKELIHOOD
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
                 lower = c(0.01),
                 upper = c(0.99))
  geo_par <- c(attributes(summary(mle_geo))$coef,attributes(summary(mle_geo))$m2logL)
  mle_pois <- mle(minus_log_like_pois,
                  start = list(lambda = lambda),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  pois_par <- c(attributes(summary(mle_pois))$coef,attributes(summary(mle_pois))$m2logL)
  mle_zeta <- mle(minus_log_like_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  zeta_par <- c(attributes(summary(mle_zeta))$coef,attributes(summary(mle_zeta))$m2logL)
  mle_zeta_trunc <- mle(minus_log_like_zeta_trunc,
                        start = list(gamma = 2),
                        method = "L-BFGS-B",
                        lower = c(1.0000001))
  zeta_trunc_par <- c(attributes(summary(mle_zeta_trunc))$coef,
                      attributes(summary(mle_zeta_trunc))$m2logL)
  df <- data.frame(pois_par,geo_par,zeta_par,zeta_trunc_par)
  colnames(df) <- NULL
  return(df)
}

# AIC
get_AIC <- function(m2logL,K,N){
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

# Function to create the AIC table and
# to store the best AIC
AIC_calc <- function(label, m2logL.list,df,x){
  AIC.list <- sapply(m2logL.list,get_AIC,K=1,N=length(x)) #row of the table
  mle_zeta2 <- 2*minus_log_like_zeta2()
  AIC.list <- append(AIC.list,get_AIC(mle_zeta2,0,length(x)),after = 3)
  best.AIC <- min(AIC.list) #best AIC
  AIC.list <- AIC.list-best.AIC 
  df <- rbind(df,data.frame(label,t(AIC.list)))
  return(list(AICdf = df, AICbest = best.AIC))
}

# PLOTS
# Probability functions, necessary to plot the distributions
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

# Function to create the plots
full_plot <- function(i,x,label,n.samples,deg_spec,df){
  degree <- as.numeric(names(deg_spec))  # X-axis
  vertices <- as.numeric(deg_spec)  # Y-axis
  plot(degree, vertices, type = "l", main = label, lwd = 2,
       xlab = "Degree", ylab = "Number of vertices", log = "xy") #plot of real data
  
  geo_prob <- sapply(x, geo_dist,p=df$p[i])
  lines(x,geo_prob*n.samples,col = "blue",lwd = 2) #plot of geometric distribution
  
  pois_prob <- sapply(x, pois_dist,lambda=df$lambda[i])
  lines(x,pois_prob*n.samples,col = "green",lwd = 2) #plot of Poisson distribution
  
  zeta_prob <- sapply(x, dzeta, shape = df$`gamma 1`[i])
  lines(x,zeta_prob*n.samples,col = "orange",lwd = 2) #plot of Zeta distribution
  
  zeta2_prob <- sapply(x, dzeta, shape = 2)
  lines(x,zeta2_prob*n.samples,col = "magenta",lwd = 2) #plot of Zeta (gamma = 2) distribution
  
  zetatrunc_prob <- sapply(x, zetatrunc_dist, gamma = df$`gamma 2`[i], 
                           k_max = df$`k max`[i]) 
  lines(x,zetatrunc_prob*n.samples,col = "red",lwd = 2) #plot of right-truncated Zeta distribution
  
  legend("topright", legend = c("Data","Geometric", "Poisson", "Zeta", 
                                "Zeta (Lambda=2)", "Zeta truncated","Confidence interval"), 
         col = c("black","blue","green","orange","magenta","red","red"), 
         lty = c(1,1,1,1,1,1,2), lwd = 2)
}

# Code for language networks
# Import the data and create the initial summary table
source = read.table("list_in.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

lang.df <- data.frame() #Summary table with all languages and important values

for (x in 1:nrow(source)) {
  lang.df <- write_table(source$language[x], source$file[x],lang.df)
}

colnames(lang.df) <- c("Language", "N", "Maximum degree", "M/N", "N/M")
lang.df 
print(xtable(lang.df), file = "Table languages.tex")


# Tables with estimated parameters and AIC
param.df <- data.frame() #table with best parameters
se.df <- data.frame() #table with the standard errors
AIC.df <- data.frame() #AIC table
bestAIC.list <- list() #list with best AIC for every language

for (i in 1:nrow(source)){
  x <- read.table(source$file[i], header = FALSE)$V1
  param.list <- mle_calc(lang.df$"N/M"[i],lang.df$"M/N"[i],x)
  param.df <- rbind(param.df,data.frame(source$language[i], param.list[1,],max(x)))
  se.df <- rbind(se.df,data.frame(source$language[i], param.list[2,]))
  AIC <- AIC_calc(source$language[i],param.list[3,],AIC.df,x)
  AIC.df <- AIC$AICdf
  bestAIC.list[i] <- AIC$AICbest
}

colnames(param.df) <- c("Language", "lambda", "p", "gamma 1","gamma 2","k max")
colnames(se.df) <- c("Language", "1", "2", "3", "5")
colnames(AIC.df) <- c("Language", "1", "2", "3", "4", "5")
param.df
se.df
AIC.df
print(xtable(param.df), file = "Table Parameters.tex")
print(xtable(AIC.df), file = "Table AIC.tex")


# Plot
for (i in 1:nrow(source)){
  degree_sequence = read.table(source$file[i], header = FALSE)
  degree_spectrum = table(degree_sequence) 
  x <- 1:max(degree_sequence)
  full_plot(i, x, source$language[i], nrow(degree_sequence), degree_spectrum, param.df)
  gamma_min = param.df$`gamma 2`[i] - sqrt(se.df$"5"[i])
  gamma_max = param.df$`gamma 2`[i] + sqrt(se.df$"5"[i])
  zetatrunc_prob <- sapply(x, zetatrunc_dist, gamma = gamma_min, 
                           k_max = param.df$`k max`[i])
  lines(x,zetatrunc_prob*nrow(degree_sequence),col = "red",lwd = 1, lty = 2)
  zetatrunc_prob <- sapply(x, zetatrunc_dist, gamma = gamma_max, 
                           k_max = param.df$`k max`[i])
  lines(x,zetatrunc_prob*nrow(degree_sequence),col = "red",lwd = 1, lty = 2)
  
}

# Samples from discrete distributions
# Import the data and create the initial summary table
folder_path <- "./samples_from_discrete_distributions/data"
files <- list.files(path = folder_path, full.names = TRUE)
prob_list <- c("geo 0.05","geo 0.1","geo 0.2","geo 0.4","geo 0.8","zeta 1.5","zeta 2.5",
          "zeta 2","zeta 3.5","zeta 3")

prob.df <- data.frame()

for (x in 1:length(files)) {
  prob.df <- write_table(prob_list[x], files[x],prob.df)
}

colnames(prob.df) <- c("Distribution", "N", "Maximum degree", "M/N", "N/M")
prob.df
print(xtable(prob.df), file = "Table sample.tex")

# LOG-LIKELIHOOD & AIC
param.sample.df <- data.frame() #table with best parameters
AIC.sample.df <- data.frame() #AIC table
bestAIC.sample.list <- list()#list with best AIC for every distribution

# (The cicle is quite slow, mainly because of the Gamma 1.5 distribution)
for (i in 1:length(files)){
  x <- read.table(files[i], header = FALSE)$V1
  param.list <- mle_calc(prob.df$"N/M"[i],prob.df$"M/N"[i],x)
  param.sample.df <- rbind(param.sample.df,data.frame(prob.df$Distribution[i], param.list[1,],max(x)))
  AIC <- AIC_calc(prob.df$Distribution[i],param.list[2,],AIC.sample.df,x)
  AIC.sample.df <- AIC$AICdf
  bestAIC.sample.list[i] <- AIC$AICbest
}

colnames(param.sample.df) <- c("Distribution", "lambda", "p", "gamma 1","gamma 2","k max")
colnames(AIC.sample.df) <- c("Distribution", "1", "2", "3", "4","5")

param.sample.df
AIC.sample.df
print(xtable(param.sample.df), file = "Table Parameters sample.tex")
print(xtable(AIC.sample.df), file = "Table AIC sample.tex")

# ADDITIONAL WORK
# Altmann distribution
minus_log_like_altmann <- function(gamma,delta){
  gamma * sum(log(x)) + delta * sum(x) + length(x)*log(sum(((1:length(x))^(-gamma)*exp(-delta*(1:length(x))))))
}

altmann_dist <- function(gamma,delta,k,N){
  k^(-gamma) * exp(-delta * k) * (1/sum((1:N)^(-gamma) * exp(-delta * (1:N))))
}

altmann.param <- data.frame()

for (i in 1:nrow(source)){
  degree_seq <- read.table(source$file[i], header = FALSE)
  degree_spec = table(degree_seq)
  x <- degree_seq$V1
  mle_altmann <- mle(minus_log_like_altmann,
                     start = list(gamma = 2, delta = 0.01),
                     method = "L-BFGS-B",
                     lower = c(0.00001,0.00001))
  alt.gamma <- attributes(summary(mle_altmann))$coef[1]
  alt.delta <- attributes(summary(mle_altmann))$coef[2]
  
  altmann.AIC <- get_AIC(attributes(summary(mle_altmann))$m2logL,2,length(x)) - bestAIC.list[[i]]
  altmann.param <- rbind(altmann.param,data.frame(alt.gamma,alt.delta,altmann.AIC))
  
  degree <- as.numeric(names(degree_spec))  # X-axis
  vertices <- as.numeric(degree_spec)  # Y-axis
  plot(degree, vertices, type = "l", main = source$language[i], lwd = 2,
       xlab = "Degree", ylab = "Number of vertices", log = "xy")
  
  z <- 1:max(degree_seq)
  altmann_prob <- sapply(z, altmann_dist,gamma=alt.gamma,delta=alt.delta,N=max(x))
  lines(z,altmann_prob*nrow(degree_seq),col = "green",lwd = 3)
  
  zetatrunc_prob <- sapply(z, zetatrunc_dist, gamma = param.df$`gamma 2`[i], 
                           k_max = param.df$`k max`[i])
  lines(z,zetatrunc_prob*nrow(degree_seq),col = "red",lwd = 3)
  legend("topright", legend = c("Data","Altmann","Zeta truncated"), 
         col = c("black","green","red"), lty = 1, lwd = 2)
}

colnames(altmann.param) <- c("gamma", "delta", "Alt.AIC-best.AIC")

altmann.param
print(xtable(altmann.param), file = "Table Altmann.tex", digits = c(0, 2, 3, 2))
