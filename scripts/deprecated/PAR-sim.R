#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Parallelized simulation of mutation accumulation in the PAR
# -----------------------------------------------------------------------------

# Load required package
library(parallel)

# -----------------------------------------------------------------------------
# 0) Helper: Determine if an individual is male or female.
# -----------------------------------------------------------------------------
IsMale <- function(individual) {
  any(individual[, 1] == 1)
}

# -----------------------------------------------------------------------------
# 1) GetMaleGametes()
#    - Generate N gametes from the male subpopulation.
#    - Each male is chosen with probability proportional to its fitness.
#    - Recombination is restricted to columns mincol..L.
# -----------------------------------------------------------------------------
GetMaleGametes <- function(malePop, wMale, mincol, N) {
  gampool <- vector("list", N)
  parents <- sample(seq_along(malePop), size = N, replace = TRUE, prob = wMale)
  for (i in seq_len(N)) {
    mat <- malePop[[parents[i]]]
    L <- ncol(mat)
    start_row  <- sample(1:2, 1)
    random_col <- sample(mincol:L, 1)
    result     <- numeric(L)
    if (random_col > 1) {
      result[1:(random_col - 1)] <- mat[start_row, 1:(random_col - 1)]
    }
    result[random_col:L] <- mat[3 - start_row, random_col:L]
    gampool[[i]] <- result
  }
  gampool
}

# -----------------------------------------------------------------------------
# 2) GetFemaleGametes()
#    - Generate N gametes from the female subpopulation.
#    - Each female is chosen with probability proportional to its fitness.
#    - Recombination can occur anywhere from 1..L.
# -----------------------------------------------------------------------------
GetFemaleGametes <- function(femalePop, wFemale, N) {
  gampool <- vector("list", N)
  parents <- sample(seq_along(femalePop), size = N, replace = TRUE, prob = wFemale)
  for (i in seq_len(N)) {
    mat <- femalePop[[parents[i]]]
    L <- ncol(mat)
    start_row  <- sample(1:2, 1)
    random_col <- sample(1:L, 1)
    result     <- numeric(L)
    if (random_col > 1) {
      result[1:(random_col - 1)] <- mat[start_row, 1:(random_col - 1)]
    }
    result[random_col:L] <- mat[3 - start_row, random_col:L]
    gampool[[i]] <- result
  }
  gampool
}

# -----------------------------------------------------------------------------
# 3) CombineGametes()
# -----------------------------------------------------------------------------
CombineGametes <- function(spermPool, eggPool) {
  N <- length(spermPool)
  newpop <- vector("list", N)
  spermOrder <- sample(seq_len(N))
  eggOrder   <- sample(seq_len(N))
  for (i in seq_len(N)) {
    newpop[[i]] <- rbind(
      spermPool[[spermOrder[i]]],
      eggPool[[eggOrder[i]]]
    )
  }
  newpop
}

# -----------------------------------------------------------------------------
# 4) GetFit()
#    - Excludes the first locus (sex-determining) from fitness calculation.
# -----------------------------------------------------------------------------
GetFit <- function(pop) {
  w <- numeric(length(pop))
  for (i in seq_along(pop)) {
    mat <- pop[[i]]
    L   <- ncol(mat)
    if (L >= 2) {
      w[i] <- prod(colMeans(mat[, 2:L, drop = FALSE]))
    } else {
      w[i] <- NA_real_
    }
  }
  w
}

# -----------------------------------------------------------------------------
# 5) MutPop()
#    - Prevents mutations at column 1 (sex-determining locus).
# -----------------------------------------------------------------------------
MutPop <- function(pop, num.mutes, s) {
  mutes <- sample(seq_along(pop), num.mutes, replace = FALSE)
  for (i in mutes) {
    mat <- pop[[i]]
    L   <- ncol(mat)
    r   <- sample(1:2, 1)
    c   <- sample(2:L, 1)
    mat[r, c] <- s
    pop[[i]]  <- mat
  }
  pop
}

# -----------------------------------------------------------------------------
# 6) Single‐replicate runner: one “j” from original RunSim()
# -----------------------------------------------------------------------------
RunOneRep <- function(seed, s, N, gen, mincol, num.mutes, num_loci) {
  set.seed(seed)
  
  # initialize population
  pop <- vector("list", N)
  for (i in seq_len(N)) {
    individual     <- matrix(1, nrow = 2, ncol = num_loci)
    individual[,1] <- sample(c(0,1), 2, replace = TRUE)
    pop[[i]]       <- individual
  }
  
  # evolve for gen generations
  mfit <- numeric(gen)
  for (g in seq_len(gen)) {
    w          <- GetFit(pop)
    male_idx   <- which(sapply(pop, IsMale))
    female_idx <- setdiff(seq_along(pop), male_idx)
    
    if (length(male_idx)==0 || length(female_idx)==0) {
      mfit[g:gen] <- NA
      break
    }
    
    w_male    <- w[male_idx]
    w_female  <- w[female_idx]
    malePop   <- pop[male_idx]
    femalePop <- pop[female_idx]
    
    spermPool <- GetMaleGametes(malePop,   w_male,   mincol, N)
    eggPool   <- GetFemaleGametes(femalePop, w_female, N)
    newpop    <- CombineGametes(spermPool, eggPool)
    pop       <- MutPop(newpop, num.mutes, s)
    
    w_new     <- GetFit(pop)
    mfit[g]   <- mean(w_new, na.rm = TRUE)
  }
  
  list(mfit = mfit, pop = pop)
}

# -----------------------------------------------------------------------------
# 7) Parallel “RunSim”: one replicate per core
# -----------------------------------------------------------------------------
RunSimParallel <- function(s         = 0.03,
                           iter      = 10,
                           N         = 500,
                           gen       = 1000,
                           mincol    = 100,
                           num.mutes = 500,
                           num_loci  = 10000,
                           ncores    = detectCores() - 1) {
  
  # spawn one seed per replicate
  reps <- mclapply(
    X        = seq_len(iter),
    FUN      = RunOneRep,
    s        = s,
    N        = N,
    gen      = gen,
    mincol   = mincol,
    num.mutes= num.mutes,
    num_loci = num_loci,
    mc.cores = min(ncores, iter)
  )
  
  # stitch together fitness matrix
  mfitMat  <- do.call(rbind, lapply(reps, `[[`, "mfit"))
  # keep the final-pop from the last replicate
  finalPop <- reps[[iter]][["pop"]]
  
  list(result = mfitMat, pop = finalPop)
}

# -----------------------------------------------------------------------------
# 8) Sweep over a vector of mincol values & collect results
# -----------------------------------------------------------------------------
# Set your simulation parameters here:
s         <- 0.03
iter      <- 10
N         <- 500
gen       <- 1000
num.mutes <- 500
num_loci  <- 10000
ncores    <- detectCores() - 1

# The different mincol breakpoints you want to try:
mincolVec <- c(10, 100, 1000, 5000, 7500, 9990)

# Run them all and store in a named list:
results <- lapply(mincolVec, function(m) {
  RunSimParallel(
    s         = s,
    iter      = iter,
    N         = N,
    gen       = gen,
    mincol    = m,
    num.mutes = num.mutes,
    num_loci  = num_loci,
    ncores    = ncores
  )
})

names(results) <- paste0("mincol.", mincolVec)

# -----------------------------------------------------------------------------
# Now `results$mincol.10`, `results$mincol.5000`, and `results$mincol.9990`
# each contain:
#   • $result: an iter×gen matrix of average fitness trajectories
#   • $pop:    the final population list from that mincol setting
# -----------------------------------------------------------------------------

# (Optional) save to disk for later analysis:
# save(results, file = "simulation_results.RData")
