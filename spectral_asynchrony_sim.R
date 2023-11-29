# Set parameters for the simulation
n_years <- 5
n_rows <- 9
n_bands <- 10

# Define a function to add correlated noise to a base vector.
# This function creates a blend of the base values and random noise.
add_correlated_noise <- function(base, noise_scale, correlation_factor) {
  base * correlation_factor + runif(length(base), min=-noise_scale, max=noise_scale) * (1 - correlation_factor)
}

# Define a function to create a single matrix, either correlated or independent.
# For correlated matrices, it applies 'add_correlated_noise' to each row.
# For independent matrices, it fills the matrix with uniform random values.
create_matrix <- function(correlated = TRUE, base_values = NULL, noise_scale = 0.1, correlation_factor = 0.8) {
  matrix_ <- matrix(nrow = n_rows, ncol = n_years)
  if (correlated) {
    for (i in 1:n_rows) {
      matrix_[i, ] <- add_correlated_noise(base_values, noise_scale, correlation_factor)
    }
  } else {
    matrix_ <- matrix(runif(n_years * n_rows), nrow = n_rows)
  }
  return(matrix_)
}

# Initialize 3D arrays for storing the correlated and independent matrices.
correlated_array <- array(dim = c(n_rows, n_years, n_bands))
independent_array <- array(dim = c(n_rows, n_years, n_bands))

# Populate the arrays with either correlated or independent matrices.
for (band in 1:n_bands) {
  base_values <- runif(n_years)  # Generate new base values for each band
  correlated_array[,,band] <- create_matrix(correlated = TRUE, base_values, noise_scale = 0.5, correlation_factor = 0.7)
  independent_array[,,band] <- create_matrix(correlated = FALSE)
}

# Function to calculate and store correlation values for each band in a 3D array.
calculate_correlations <- function(array_3d) {
  n_bands <- dim(array_3d)[3]
  n_pixels <- dim(array_3d)[1]
  
  # Initialize a matrix to store the correlation values.
  correlation_values_matrix <- matrix(nrow = n_pixels * (n_pixels - 1) / 2, ncol = n_bands)
  
  # Loop through each band to calculate correlations.
  for (band in 1:n_bands) {
    # Calculate the correlation matrix for the transpose of the band matrix.
    corr_matrix <- cor(t(array_3d[,,band]))
    
    # Store the upper triangle of the correlation matrix in the matrix.
    correlation_values_matrix[,band] <- corr_matrix[upper.tri(corr_matrix)]
  }
  
  return(correlation_values_matrix)
}


# Define a function to calculate spectral asynchrony from correlation values.
# This function uses Fisher z-transformation to normalize correlation values,
# calculates the mean z-score for each band and an overall mean,
# and then back-transforms the mean to a correlation value.
spectral_asynchrony <- function(x){
  cor_matrix <- calculate_correlations(x)
  rz <- 0.5 * log((1 + cor_matrix)/(1 - cor_matrix))
  rz_mean <- colMeans(rz)
  spec_asyn <- 1 - (exp(2 * mean(rz_mean)) - 1)/(exp(2 * mean(rz_mean)) + 1)
  return(spec_asyn)
}

# Calculate spectral asynchrony for the correlated matrix correlations
spectral_asynchrony(correlated_array)

# Calculate spectral asynchrony for the independent matrix correlations
spectral_asynchrony(independent_array)
