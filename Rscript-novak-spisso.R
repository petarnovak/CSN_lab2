# 1: Introduction
write_table <- function(label,file,df) {
   degree_sequence = read.table(file, header = FALSE)
   lang.df <- rbind(df,data.frame(label, length(degree_sequence$V1), max(degree_sequence$V1),
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

# 2: Visualization
# Function to plot the degree sequence of a given language
degree_plot <- function(language,file){
  degree_sequence = read.table(file, header = FALSE)
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum/sum(degree_spectrum), main = language, xlab = "degree", ylab = "percentage of vertices", log = "y")
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
                 lower = c(0.0000001),
                 upper = c(0.9999999))
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

AIC_calc <- function(label, m2logL.list,df){
  AIC.list <- sapply(m2logL.list,get_AIC,K=1,N=length(x))
  mle_zeta2 <- 2*minus_log_like_zeta2()
  AIC.list <- append(AIC.list,get_AIC(mle_zeta2,0,length(x)),after = 3)
  best.AIC <- min(AIC.list)
  AIC.list <- AIC.list-best.AIC
  df <- rbind(df,data.frame(label,t(AIC.list)))
  return(df)
}

param.df <- data.frame() #table with best parameters
AIC.df <- data.frame() #AIC table

for (i in 1:nrow(source)){
  x <- read.table(source$file[i], header = FALSE)$V1
  param.list <- mle_calc(lang.df$"N/M"[i],lang.df$"M/N"[i],x)
  param.df <- rbind(param.df,data.frame(source$language[i], param.list[1,],max(x)))
  AIC.df <- AIC_calc(source$language[i],param.list[2,],AIC.df)
}

colnames(param.df) <- c("Language", "lambda", "p", "gamma 1","gamma 2","k max")
colnames(AIC.df) <- c("Language", "1", "2", "3", "4", "5")

param.df
AIC.df

# Plots of distributions
full_plot <- function(i,language,file){
  degree_sequence = read.table(file, header = FALSE)
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum/sum(degree_spectrum), main = language, xlab = "degree", ylab = "number of vertices")
  
  x_vals <- 1:rownames(degree_spectrum)[nrow(degree_spectrum)]
  best_params <- param.df[i,]
  geo <- dgeom(x_vals, prob = best_params$p)
  lines(x_vals,geo)
  #pois <- dpois(x_vals, lambda = best_params$lambda)
}

for (x in 1:nrow(source)){
  full_plot(x,source$language[x], source$file[x])
}

#plot(x_vals, pmf_vals, type = "l", lwd = 2, col = "red",
#     main = "PMF of Geometric Distribution (p = 0.2)", xlab = "Number of Failures", ylab = "Probability")


# 6: Samples from discrete distribution
folder_path <- "./samples_from_discrete_distributions/data"
files <- list.files(path = folder_path, full.names = TRUE)
prob_list <- c("geo 0.05","geo 0.1","geo 0.2","geo 0.4","geo 0.8","zeta 1.5","zeta 2.5",
          "zeta 2","zeta 3.5","zeta 3")

# Visualization
for (i in 1:length(files)){
  degree_sequence = read.table(files[i], header = FALSE)
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum, main = paste("Distribution :",prob_list[i]), xlab = "degree", ylab = "number of vertices", log = "y")
}

#
# WORK IN PROGRESS
#
prob.df <- data.frame()

for (x in 1:length(files)) {
  prob.df <- write_table(prob_list[x], files[x],prob.df)
}

colnames(prob.df) <- c("Distribution", "N", "Maximum degree", "M/N", "N/M")
prob.df

mle_calc <- function(i,prob,file,param.df,AIC.df){
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
  
  x <- read.table(file, header = FALSE)$V1
  mle_geo <- mle(minus_log_like_geo,
                 start = list(p = prob.df$"N/M"[i]),
                 method = "L-BFGS-B",
                 lower = c(0.01), #problems with zeta: gamma = 1.5, with smaller values it crashes
                 upper = c(0.99))
  mle_pois <- mle(minus_log_like_pois,
                  start = list(lambda = prob.df$"M/N"[i]),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  mle_zeta <- mle(minus_log_like_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  mle_zeta_trunc <- mle(minus_log_like_zeta_trunc,
                        start = list(gamma = 2),
                        method = "L-BFGS-B",
                        lower = c(1.0000001))
  
  # 4: Finding the best models
  best_p_geo <- attributes(summary(mle_geo))$coef[1]
  best_lambda_pois <- attributes(summary(mle_pois))$coef[1]
  best_gamma_zeta <- attributes(summary(mle_zeta))$coef[1]
  best_gamma_zeta_trunc <- attributes(summary(mle_zeta_trunc))$coef[1]
  best_h_max_zeta_trunc <- max(x)
  param.df <- rbind(param.df,data.frame(prob, best_lambda_pois, best_p_geo, 
                                        best_gamma_zeta, best_gamma_zeta_trunc,
                                        best_h_max_zeta_trunc))
  
  # 5: Best model selection
  get_AIC <- function(m2logL,K,N){
    m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
  }
  
  AIC_list <- c()
  AIC_list[1] <- get_AIC(attributes(summary(mle_geo))$m2logL,1,length(x))
  AIC_list[2] <- get_AIC(attributes(summary(mle_pois))$m2logL,1,length(x))
  AIC_list[3] <- get_AIC(attributes(summary(mle_zeta))$m2logL,1,length(x))
  mle_zeta2 <- 2*minus_log_like_zeta2()
  AIC_list[4] <- get_AIC(mle_zeta2,0,length(x))
  AIC_list[5] <- get_AIC(attributes(summary(mle_zeta_trunc))$m2logL,1,length(x))
  
  best_AIC <- min(AIC_list)
  AIC_list <- AIC_list-best_AIC
  
  AIC.df <- rbind(AIC.df,data.frame(prob,t(AIC_list)))
  return(list(param.df, AIC.df))
}

param.df <- data.frame() #table with most likely parameters
AIC.df <- data.frame()

for (i in 1:length(files)){
  tables <- mle_calc(i, prob_list[i], files[i],param.df,AIC.df)
  param.df <- tables[[1]]
  AIC.df <- tables[[2]]
}

colnames(param.df) <- c("Distribution", "lambda", "p", "gamma 1","gamma 2","k max")
colnames(AIC.df) <- c("Distribution", "1", "2", "3", "4","5")

param.df
AIC.df
