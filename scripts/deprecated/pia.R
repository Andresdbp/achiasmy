# Pairwise Invasibility Analysis

library(RColorBrewer)


compute_lambda <- function(S_B, L, S_d = 0.03, N = 101) {
  U <- 1e-8 * L  # Genomic mutation rate.
  A <- matrix(0, nrow = N, ncol = N)
  
  # Loop over mutation classes k = 0 to N-2.
  for (k in 0:(N - 2)) {
    f_k <- (1 + S_B) * ((1 - S_d)^k)
    # Offspring remaining in class k.
    A[k + 1, k + 1] <- f_k * (1 - U)
    # Offspring acquiring one extra deleterious mutation (moving to class k+1).
    A[k + 2, k + 1] <- f_k * U
  }
  
  # For the highest mutation class (k = N - 1), assume an absorbing state.
  k <- N - 1
  f_k <- (1 + S_B) * ((1 - S_d)^k)
  A[N, N] <- f_k
  
  lambda <- max(Re(eigen(A)$values))
  return(lambda)
}

# Define parameter sequences.
S_b_values <- seq(0, 0.05, length.out = 100)
L_values   <- seq(0, 5e6, length.out = 100)

# Build the lambda matrix:
# Rows correspond to different L values and columns to S_b values.
lambda_matrix <- matrix(NA, nrow = length(L_values), ncol = length(S_b_values))
for (i in seq_along(L_values)) {
  print(i)
  for (j in seq_along(S_b_values)) {
    lambda_matrix[i, j] <- compute_lambda(S_b_values[j], L_values[i], S_d = 0.03, N = 101)
  }
}

diverge_pal <- colorRampPalette(brewer.pal(11, "RdBu"))(1000)

breaks_vec <- seq(0.95, 1.05, length.out = 1001)

# plotting.
image(y = S_b_values, x = L_values, z = lambda_matrix,
      breaks = breaks_vec, col = diverge_pal,
      ylab = expression(S[b]),
      xlab = "Recombining region size, L (Mb)",
      main = expression("PIA of " * lambda * " vs. " * S[b] * " and L"),
      yaxt='n', xaxt='n',cex.lab = 1.6, cex.main = 1.8)
axis(2, at = seq(from = 0, to = 0.05, length.out=6),
     labels = seq(from = 0, to = 0.05, length.out=6),
     cex.axis = 1.5)
axis(1, at = seq(from = 0, to = 5e6, length.out=6),
     labels = seq(from = 0, to = 5, length.out=6),
     cex.axis = 1.5)
contour(y = S_b_values, x = L_values, z = lambda_matrix, 
        add = TRUE, drawlabels = TRUE, col = "black")

