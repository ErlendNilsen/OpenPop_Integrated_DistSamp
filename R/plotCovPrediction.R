#' Plot vital rates as functions of fitted covariates
#'
#' Note that per now, this only supports covariates included in hierarchical
#' models for recruitment rate. However, the function can easily be extended
#' to also work for covariate effects on survival.
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param effectParam character string naming the effect parameter to use for plotting.
#' @param covName character string that defines the covariate.
#' @param minCov numeric. Minimum covariate value to predict for.
#' @param maxCov numeric. Maximum covariate value to predict for.
#' @param meanCov numeric. Mean value of the original covariate (used for z-
#' standardization).
#' @param sdCov numeric. Standard deviation of the original covariate (used for 
#' z-standardization).
#' @param covData numeric. Matrix of covariate values used in analysis (columns = years, rows = areas).
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param fitRodentCov logical. If TRUE, plots are made for predictions with 
#' rodent covariate. If FALSE, no plots are made. 
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/CovPredictions.
#' @export
#'
#' @examples

plotCovPrediction <- function(mcmc.out, 
                              effectParam, covName,
                              minCov, maxCov,
                              meanCov, sdCov,
                              covData,
                              N_areas, area_names,
                              fitRodentCov){
  
  if(fitRodentCov){
    ## Make sequence of absolute and z-standardized covariate values to predict for
    cov_abs <- seq(minCov, maxCov, length.out = 100)
    cov <- (cov_abs - meanCov) / sdCov
  
    ## Assemble dataframe for storing posterior summaries of predictions
    cov.pred.data <- data.frame()
    
    for(i in 1:N_areas){
      
      ## Extract posterior samples of relevant effect slope and intercept
      beta <- as.matrix(mcmc.out)[, paste0(effectParam, "[", i, "]")]
      Mu.R <- as.matrix(mcmc.out)[, paste0("Mu.R[", i, "]")]
      
      ## Make, summarise, and store predictions for each covariate value
      for(x in 1:100){
        R.pred <- exp(log(Mu.R) + beta*cov[x])
        cov.pred.temp <- data.frame(Area = area_names[i], 
                                    covValue_z = cov[x], 
                                    covValue = cov_abs[x], 
                                    pred_Median = median(R.pred),
                                    pred_lCI = unname(quantile(R.pred, probs = 0.025)),
                                    pred_uCI = unname(quantile(R.pred, probs = 0.975)))
        cov.pred.data <- rbind(cov.pred.data, cov.pred.temp)
      }
    }
    
    ## Write covariate data into a data frame
    cov.raw.data <- data.frame()
    for(i in 1:N_areas){
      
      # Assemble data
      raw.data <- data.frame(Area = area_names[i],
                             covValue_z = covData[i,],
                             covValue = (covData[i,] * sdCov) + meanCov,
                             yearIdx = 1:length(covData[i,]))
      # Count number of times same value appears (for spacing in plot)
      raw.data <- raw.data %>%
        dplyr::add_count(covValue) %>%
        dplyr::mutate(occur_n = 1)

      for(j in 1:nrow(raw.data)){
        if(raw.data$n[j] > 1){
          raw.data$occur_n[j] <- nrow(subset(raw.data[1:j,], covValue == covValue[j]))
        }
      }
      cov.raw.data <- rbind(cov.raw.data, raw.data)
    }
    
    
    ## Plot predictions
    ifelse(!dir.exists("Plots/CovPredictions"), dir.create("Plots/CovPredictions"), FALSE) ## Check if folder exists, if not create folder
    
    pdf(paste0("Plots/CovPredictions/Rep_", effectParam, ".pdf"), width = 6, height = 4) 
    for(i in 1:N_areas){
      pred.max <- max(subset(cov.pred.data, Area == area_names[i])$pred_uCI)
      sub.cov.data <- subset(cov.raw.data, Area == area_names[i]) %>%
        dplyr::mutate(y_pos = pred.max -(0.1*(occur_n-1)))
      
      print(
        ggplot(subset(cov.pred.data, Area == area_names[i]), aes(x = covValue, y = pred_Median)) +
          geom_line(color = "forestgreen") + 
          geom_ribbon(aes(ymin = pred_lCI, ymax = pred_uCI), alpha = 0.5, fill = "forestgreen") + 
          geom_point(data = sub.cov.data, aes(x = covValue, y = y_pos)) + 
          xlab(covName) +
          ylab("Recruitment rate") + 
          ggtitle(area_names[i]) + 
          theme_classic()
      )
    }
    dev.off()
    
    ## Return plot paths
    plot.paths <- paste0("Plots/CovPredictions/Rep_", effectParam, ".pdf")
    return(plot.paths)
    
  }else{
    message("No plots produced since fitRodentCov = FALSE (no covariate effect fitted")
  }
}