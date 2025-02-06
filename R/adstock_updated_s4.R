install.packages(readLines("requirements.R"))

setClass("AdstockModel", 
         slots = list(data = "data.frame", results = "data.frame", dates = "Date", saturation_models = "list", distribution_params = "list"))

setGeneric("runAdstock", function(object) standardGeneric("runAdstock"))
setGeneric("computeCorrelations", function(object) standardGeneric("computeCorrelations"))

setMethod("runAdstock", "AdstockModel", function(object) {
  library(data.table)
  library(dplyr)
  
  # Extract the input data and dates
  data <- object@data
  dates <- object@dates
  dist_params <- object@distribution_params
  
  # Ensure all variables are numeric except for date
  data <- data %>% mutate(across(where(is.character), as.numeric))
  
  # Define Lag Function
  lagpad <- function(x, k) {
    if (k > 0) return(c(rep(NA, k), x)[1:length(x)])
    else return(c(x[(-k + 1):length(x)], rep(NA, -k)))
  }
  
  # Define Normalization Function
  normalize <- function(x) x / sum(x, na.rm = TRUE)
  
  # Define Saturation Models
  saturationModels <- list(
    linear = function(x) x,
    diminishing = function(x) x / (1 + x),
    log = function(x) log(1 + x),
    exponential = function(x) 1 - exp(-x)
  )
  
  # Define Standardization Function
  ztran <- function(x) {
    mns <- colMeans(x, na.rm = TRUE)
    sds <- apply(x, 2, sd, na.rm = TRUE)
    sweep(sweep(x, 2, mns, "-"), 2, sds, "/")
  }
  
  # Generate Distributions
  generateDistributions <- function(dist, params, func) {
    t <- data.frame(X = seq(1, 104))
    for (i in seq(params$param1[1], params$param1[2], params$param1[3])) {
      for (j in seq(params$param2[1], params$param2[2], params$param2[3])) {
        t[[paste0(dist, "_", i, j)]] <- sapply(t$X, function(x) func(x, i, j))
      }
    }
    return(t)
  }
  
  distributions <- list(
    Weibull = generateDistributions("Weibull", dist_params$Weibull, dweibull),
    Lognormal = generateDistributions("lnorm", dist_params$Lognormal, dlnorm),
    Logistic = generateDistributions("logis", dist_params$Logistic, dlogis),
    NegExp = generateDistributions("exp", dist_params$NegExp, dexp),
    ChiSquare = generateDistributions("chisq", dist_params$ChiSquare, dchisq)
  )
  
  final_results <- data.frame()
  
  for (dist_name in names(distributions)) {
    dist_data <- distributions[[dist_name]]
    dist_data <- dist_data %>% mutate(across(-X, normalize))
    dist_data <- rbind(dist_data, matrix(0, 104, ncol(dist_data)))
    
    for (sat_name in names(saturationModels)) {
      # Lagging and Multiplication with Media Variables
      result_data <- data.frame(X = dist_data$X)
      for (m in 2:ncol(dist_data)) {
        mult <- as.data.frame(dist_data[[m]])
        for (i in 1:208) {
          mult[[paste0("Week_", i)]] <- lagpad(mult[[1]], i - 1)
        }
        mult[[1]] <- NULL
        mult[is.na(mult)] <- 0
        transformed_data <- rowSums(t(t(mult) * data[[m]]))
        result_data[[paste0(colnames(dist_data)[m], "_", sat_name)]] <- saturationModels[[sat_name]](transformed_data)
      }
      result_data <- ztran(result_data)
      result_data$Week <- dates
      
      object@data <- result_data
      result_data <- computeCorrelations(object)
      
      result_data$Vehicle <- paste0(dist_name, "_", sat_name)
      final_results <- rbind(final_results, result_data)
    }
  }
  
  object@results <- final_results
  return(object)
})

setMethod("computeCorrelations", "AdstockModel", function(object) {
  df <- data.frame(Variable = character(), Correl = numeric(), Pval = numeric())
  
  # Compute correlation between media variables and EQ/Volume
  for (i in 2:ncol(object@data)) {
    cor_result <- cor.test(object@data[[i]], object@data[[1]])  # First column assumed as EQ/Volume
    df <- rbind(df, data.frame(Variable = colnames(object@data)[i], 
                               Correl = cor_result$estimate, 
                               Pval = cor_result$p.value))
  }
  
  df <- df[complete.cases(df) & df$Correl > 0, ]
  df <- df[order(df$Correl, decreasing = TRUE), ]
  return(df)
})


