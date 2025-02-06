setClass("AdstockModel", 
         slots = list(data = "data.frame", eq_data = "data.frame", results = "data.frame", dates = "Date"))

setGeneric("runAdstock", function(object) standardGeneric("runAdstock"))
setGeneric("computeCorrelations", function(object) standardGeneric("computeCorrelations"))

setMethod("runAdstock", "AdstockModel", function(object) {
  library(data.table)
  library(dplyr)
  library(readxl)
  
  # Data Processing
  EQ <- object@eq_data
  EQ$ZIP <- NULL
  EQ <- EQ[EQ$Date > object@dates, ]
  EQ <- EQ %>% group_by(TDLINXID) %>% mutate(EQ_z = scale(EQ_proj))
  EQ$EQ_proj <- NULL
  EQ[is.na(EQ)] <- 0
  
  med_data <- object@data
  dates <- object@dates
  med_data$Week <- as.Date(med_data$Week)
  
  # Define Lag Function
  lagpad <- function(x, k) {
    if (k > 0) return(c(rep(NA, k), x)[1:length(x)])
    else return(c(x[(-k + 1):length(x)], rep(NA, -k)))
  }
  
  # Define Normalization Function
  normalize <- function(x) x / sum(x, na.rm = TRUE)
  
  # Define Standardization Function
  ztran <- function(x) {
    mns <- colMeans(x, na.rm = TRUE)
    sds <- apply(x, 2, sd, na.rm = TRUE)
    sweep(sweep(x, 2, mns, "-"), 2, sds, "/")
  }
  
  # Generate Distributions
  generateDistributions <- function(dist, param1, param2, func) {
    t <- data.frame(X = seq(1, 104))
    for (i in seq(param1[1], param1[2], param1[3])) {
      for (j in seq(param2[1], param2[2], param2[3])) {
        t[[paste0(dist, "_", i, j)]] <- sapply(t$X, function(x) func(x, i, j))
      }
    }
    return(t)
  }
  
  distributions <- list(
    Weibull = generateDistributions("Weibull", c(0.1, 5, 0.1), c(1, 70, 1), dweibull),
    Lognormal = generateDistributions("lnorm", c(1, 20, 1), c(10, 20, 1), dlnorm),
    Logistic = generateDistributions("logis", c(1, 16, 1), c(1, 16, 1), dlogis),
    NegExp = generateDistributions("exp", c(0.1, 1, 0.1), c(NA, NA, NA), dexp),
    ChiSquare = generateDistributions("chisq", c(1, 70, 1), c(NA, NA, NA), dchisq)
  )
  
  final_results <- data.frame()
  
  for (dist_name in names(distributions)) {
    dist_data <- distributions[[dist_name]]
    dist_data <- dist_data %>% mutate(across(-X, normalize))
    dist_data <- rbind(dist_data, matrix(0, 104, ncol(dist_data)))
    
    # Lagging and Multiplication with Med Data
    result_data <- data.frame(X = dist_data$X)
    for (m in 2:ncol(dist_data)) {
      mult <- as.data.frame(dist_data[[m]])
      for (i in 1:208) {
        mult[[paste0("Week_", i)]] <- lagpad(mult[[1]], i - 1)
      }
      mult[[1]] <- NULL
      mult[is.na(mult)] <- 0
      result_data[[colnames(dist_data)[m]]] <- rowSums(t(t(mult) * med_data[[m]]))
    }
    result_data <- ztran(result_data)
    result_data$Week <- dates
    
    object@data <- result_data
    object@eq_data <- EQ
    result_data <- computeCorrelations(object)
    
    result_data$Vehicle <- dist_name
    final_results <- rbind(final_results, result_data)
  }
  
  object@results <- final_results
  return(object)
})

setMethod("computeCorrelations", "AdstockModel", function(object) {
  df <- data.frame(Variable = character(), Correl = numeric(), Pval = numeric())
  merged_data <- merge(object@eq_data, object@data, by.x = "Date", by.y = "Week")
  
  for (i in 5:ncol(merged_data)) {
    merged_data <- merged_data %>% group_by(TDLINXID) %>% mutate_at(vars(i), scale)
    merged_data[is.na(merged_data)] <- 0
    cor_result <- cor.test(merged_data[[i]], merged_data$EQ_z)
    df <- rbind(df, data.frame(Variable = colnames(merged_data)[i], 
                               Correl = cor_result$estimate, 
                               Pval = cor_result$p.value))
  }
  
  df <- df[complete.cases(df) & df$Correl > 0, ]
  df <- df[order(df$Correl, decreasing = TRUE), ]
  return(df)
})

# Run the Adstock Model with Input Data and Dates
adstock_model <- new("AdstockModel", data = input_data, eq_data = eq_input_data, dates = input_dates)
adstock_model <- runAdstock(adstock_model)
write.csv(adstock_model@results, "adstock_new_iterations.csv")
