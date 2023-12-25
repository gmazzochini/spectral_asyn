# Load rainfall data from 1980 to 2018 into a matrix.
rainfall_1980_2018 <-  matrix(c(540.7, 542.4, 308.9, 415.6, 562, 956.4, 664.3, 425.5, 651.7, 786.3, 292.8, 360.9, 509.7, 
                                215.8, 486.9, 455.3, 526.5, 458.9, 387.2, 442.5, 609.4, 474.7, 560.8, 488.6, 525.4, 506.8, 
                                588.9, 442.5, 715.4, 581.2, 482.5, 615.6, 267.6, 426, 609.8, 399.6, 376.7, 382.9, 448.5),
                              ncol = 39, nrow = 1, dimnames = list(1, 1980:2018))

# Load EVI data from 2001 to 2018 into a matrix.
ievi_2001_2018 <- matrix(c(3.228700, 3.564800, 3.537150, 3.708675, 3.845050, 4.362100, 3.120375, 4.270150, 4.386000, 
                           3.881350, 4.093050, 2.269400, 2.486450, 3.179300, 2.827050, 3.254400, 2.524025, 3.034800),
                         ncol = 18, nrow = 1, dimnames = list(1, 2001:2018))

# Define a function to calculate the resistance and recovery of EVI during and after a drought year.
drought_response <- function(rain_data, evi_data){
  
  # Initialize a data frame to hold the resistance and recovery values.
  response_df <- data.frame(resistance = NA,recovery = NA)
  
  for(i in 1:nrow(rain_data)){
    # Scale the logarithm of the rainfall data.
    rainfall_scaled <- scale(log(rain_data[i,]))
    
    # Identify drought years as those with scaled rainfall values less than -1.5.
    dry_years <- which(rainfall_scaled[colnames(evi_data),] < -1.5)
    
    # Scale the EVI data and convert it to a vector.
    y <- as.vector(scale(evi_data[i,]))
    
    # Scale the logarithm of the rainfall data for EVI data years and convert it to a vector.
    x <- as.vector(scale(log(rain_data[i,colnames(evi_data)])))
    
    # Create a dataframe for modeling.
    d1 <- data.frame(x,y)
    
    # Fit a generalized least squares (GLS) model with an autoregressive correlation structure.
    m <- gls(y ~ x, correlation = corAR1(form = ~1), data = d1[-dry_years,])
    
    # Calculate and return resistance and recovery.
    # Resistance: Mean difference between observed EVI in drought years and the model's prediction.
    # Recovery: Mean difference between observed EVI in the year following drought and the model's prediction.
    response_df[i, "resistance"] <- mean(y[dry_years] - predict(m,d1[dry_years,]),na.rm = T)
    response_df[i, "recovery"] <- mean(y[dry_years+1] - predict(m,d1[dry_years+1,]),na.rm = T)
  }
  return(response_df)
}

# Execute the drought response function with rainfall and EVI data.
drought_response(rain_data = rainfall_1980_2018,
                 evi_data = ievi_2001_2018)
