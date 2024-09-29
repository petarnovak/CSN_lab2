# 1: Introduction
write_table <- function(language,file,lang.df) {
   degree_sequence = read.table(file, header = FALSE)
   lang.df <- rbind(lang.df,data.frame(language, length(degree_sequence$V1), max(degree_sequence$V1),
                                       sum(degree_sequence$V1)/length(degree_sequence$V1),
                                       length(degree_sequence$V1)/sum(degree_sequence$V1)))
   return(lang.df)
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
  barplot(degree_spectrum, main = language, xlab = "degree", ylab = "number of vertices", log = "y")
}

for (x in 1:nrow(source)){
  degree_plot(source$language[x], source$file[x])
}

# 3: Log-likelihood function
require(stats4) # for MLE
require(VGAM) # for the Riemann-zeta function

# 4: Finding the best models
mle_calc <- function(i,language,file,param.df){
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
    length(x) * log(pi^2/6) + 2 * sum(log(x))
  }
  
  # Right-truncated zeta distribution
  minus_log_like_zeta_trunc <- function(gamma,h_max){
    length(x) * log(sum((1:h_max)^(-gamma))) + gamma * sum(log(x)) 
  }
  
  x <- read.table(file, header = FALSE)$V1
  mle_geo <- mle(minus_log_like_geo,
                 start = list(p = lang.df$"N/M"[i]),
                 method = "L-BFGS-B",
                 lower = c(0.0000001),
                 upper = c(0.9999999))
  mle_pois <- mle(minus_log_like_pois,
                  start = list(lambda = lang.df$"M/N"[i]),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  mle_zeta <- mle(minus_log_like_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  #This gives error for some languages
  #mle_zeta_trunc <- mle(minus_log_like_zeta_trunc,
  #                      start = list(gamma = 2, h_max = max(x)),
  #                      method = "L-BFGS-B",
  #                      lower = c(1.0000001, 1),
  #                      upper = c(Inf, max(x)+1))
  best_p_geo <- attributes(summary(mle_geo))$coef[1]
  best_lambda_pois <- attributes(summary(mle_pois))$coef[1]
  best_gamma_zeta <- attributes(summary(mle_zeta))$coef[1]
  #best_gamma_zeta_trunc <- attributes(summary(mle_zeta_trunc))$coef[1]
  #best_h_max_zeta_trunc <- attributes(summary(mle_zeta_trunc))$coef[2]
  param.df <- rbind(param.df,data.frame(language, best_lambda_pois, best_p_geo, best_gamma_zeta))
  return(param.df)
}

param.df <- data.frame() #table with most likely parameters

for (i in 1:nrow(source)){
  param.df <- mle_calc(i, source$language[i], source$file[i],param.df)
}

colnames(param.df) <- c("Language", "lambda", "p", "gamma_1")
#Table with the best parameters for every language
param.df

# 5: Best model selection
get_AIC <- function(m2logL,K,N){
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

# 6: Samples from discrete distribution
folder_path <- "./samples_from_discrete_distributions/data"
files <- list.files(path = folder_path, full.names = TRUE)
# Visualization (TO ADD FILE NAMES)
for (file in files){
  degree_sequence = read.table(file, header = FALSE)
  degree_spectrum = table(degree_sequence)
  barplot(degree_spectrum, xlab = "degree", ylab = "number of vertices", log = "y")
}
geom_005 <- files[1]
