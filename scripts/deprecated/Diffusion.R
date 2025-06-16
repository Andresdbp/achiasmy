#------------------------------------------------------------------
# Finite Difference Scheme for the Backward Kolmogorov Equation
# and a Parameter Sweep over s_driver and c_decline
#------------------------------------------------------------------

# Define a function that runs the simulation and returns the fixation probability vector u.
simulate_fixation <- function(s_driver, c_decline, M = 100, N_time = 1000, T_final = 1, Ne = 1000) {
  dx <- 1 / M
  dt <- T_final / N_time
  
  # Initial condition: linear profile u(x) = x for x in [0, 1]
  u <- seq(0, 1, length.out = M + 1)
  
  # Backward Euler time-stepping (integrate backward from T_final to 0)
  for (j in 1:N_time) {
    t <- T_final - j * dt
    # Compute the time-dependent selection coefficient and bound it between -10 and 10
    s_t <- max(-10, min(s_driver - c_decline * t, 10))
    
    n_int <- M - 1   # Number of interior grid points (indices 2 through M)
    b_vec <- u[2:M]  # Right-hand side from previous time step for interior points
    A <- matrix(0, n_int, n_int)  # Tridiagonal matrix for the implicit system
    
    for (i in 2:M) {
      k <- i - 1          # Interior index (1 to n_int)
      x_i <- (i - 1) * dx # Spatial coordinate for u[i]
      
      # Compute the drift and diffusion coefficients at x_i
      drift_coef <- s_t * x_i * (1 - x_i)
      diff_coef  <- (1 / (2 * Ne)) * x_i * (1 - x_i)
      
      # Finite-difference coefficients:
      # Left neighbor coefficient (for u[i-1])
      L_coef <- dt * (drift_coef / (2 * dx) - diff_coef / (dx^2))
      # Center coefficient (for u[i])
      D_coef <- 1 + 2 * dt * (diff_coef / (dx^2))
      # Right neighbor coefficient (for u[i+1])
      U_coef <- - dt * (drift_coef / (2 * dx) + diff_coef / (dx^2))
      
      # Fill the tridiagonal system
      A[k, k] <- D_coef
      if (k > 1) {
        A[k, k - 1] <- L_coef
      }
      if (k < n_int) {
        A[k, k + 1] <- U_coef
      } else {
        # For the last interior point (i = M), account for the boundary condition at x=1 (u(1)=1)
        b_vec[k] <- b_vec[k] - U_coef * 1
      }
    }
    
    # Solve the linear system to obtain the updated interior values
    u[2:M] <- solve(A, b_vec)
  }
  
  return(u)
}

#------------------------------------------------------------------
# Parameter Sweep: Define ranges for s_driver and c_decline
#------------------------------------------------------------------

c_max <- (1e-08)*(5e6)*(0.03)

# Increase the resolution of the parameter grid:
s_values <- seq(0, 0.5, length.out = 51)    # s_driver from 0 to 0.5
c_values <- seq(0, c_max, length.out = 51)      # c_decline from 0 to 0.2

# We'll record the fixation probability at a chosen allele frequency.
# With M=100 (grid from 0 to 1 in steps of 0.01), the allele frequency x = 0.1
# corresponds to index: x_index = 0.1/0.01 + 1 = 11.
fix_prob_matrix <- matrix(NA, nrow = length(s_values), ncol = length(c_values))
x_index <- 2  # index corresponding to the allele frequency of interest

# Loop over the parameter grid
for (i in seq_along(s_values)) {
  cat("Processing s_driver index:", i, "\n")
  for (j in seq_along(c_values)) {
    u_final <- simulate_fixation(s_driver = s_values[i], c_decline = c_values[j])
    fix_prob_matrix[i, j] <- u_final[x_index]
  }
}

#------------------------------------------------------------------
# Plotting: Visualize the effect of s_driver and c_decline on the fixation probability
# at the chosen allele frequency.
# We want c_decline on the x-axis and s_driver on the y-axis.
#------------------------------------------------------------------
library(RColorBrewer)

# Define the number of color levels and create level breakpoints
nlevels <- 100
levels_vec <- seq(min(t(fix_prob_matrix)), max(t(fix_prob_matrix)), length.out = nlevels)

# Image plot with the same color gradient as the filled.contour() example:
image(c_values, s_values, t(fix_prob_matrix),
      xlab = "c_decline", ylab = "s_driver",
      main = "Fixation Probability at Allele Frequency (index 2)",
      col = colorRampPalette(brewer.pal(11, "Blues"))(nlevels - 1),
      #col = viridis::mako(nlevels - 1),
      breaks = levels_vec)

