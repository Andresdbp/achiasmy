##############################################################################
# Helper: Determine if an individual is male or female.
#   - We define "male" if there's ANY '1' in the first column (row 1 or row 2).
#   - Otherwise, female.
##############################################################################
IsMale <- function(individual) {
  # individual is a 2x100 matrix
  return(any(individual[, 1] == 1))
}
##############################################################################
# 1) GetMaleGametes()
#    - Generate N gametes from the male subpopulation.
#    - Each male is chosen with probability proportional to its fitness.
#    - Recombination is restricted to columns mincol..100.
##############################################################################
GetMaleGametes <- function(malePop, wMale, mincol, N) {
  # malePop: list of individuals that are male (each 2x100)
  # wMale: fitness values for these male individuals
  # mincol: boundary for restricted recombination
  # N: number of gametes to produce (i.e., population size)
  
  gampool <- vector("list", N)  # store N male gametes
  
  # Sample parents (with replacement) according to fitness
  parents <- sample(
    x = seq_along(malePop),
    size = N,
    replace = TRUE,
    prob = wMale
  )
  
  for (i in seq_len(N)) {
    # Choose one of the two rows at random (row1 or row2)
    start_row <- sample(1:2, 1)
    mat       <- malePop[[parents[i]]]
    
    # Restricted recombination breakpoint
    random_col <- sample(mincol:ncol(mat), 1)
    
    # Build the gamete
    result <- numeric(100)
    if (random_col > 1) {
      result[1:(random_col - 1)] <- mat[start_row, 1:(random_col - 1)]
    }
    result[random_col:100] <- mat[3 - start_row, random_col:100]
    
    gampool[[i]] <- result
  }
  
  return(gampool)
}
##############################################################################
# 2) GetFemaleGametes()
#    - Generate N gametes from the female subpopulation.
#    - Each female is chosen with probability proportional to its fitness.
#    - Recombination can occur anywhere from 1..100.
##############################################################################
GetFemaleGametes <- function(femalePop, wFemale, N) {
  # femalePop: list of individuals that are female (each 2x100)
  # wFemale: fitness values for these female individuals
  # N: number of gametes to produce (i.e., population size)
  
  gampool <- vector("list", N)  # store N female gametes
  
  # Sample parents (with replacement) according to fitness
  parents <- sample(
    x = seq_along(femalePop),
    size = N,
    replace = TRUE,
    prob = wFemale
  )
  
  for (i in seq_len(N)) {
    # Choose one of the two rows at random
    start_row <- sample(1:2, 1)
    mat       <- femalePop[[parents[i]]]
    
    # Unrestricted recombination breakpoint
    random_col <- sample(1:ncol(mat), 1)
    
    # Build the gamete
    result <- numeric(100)
    if (random_col > 1) {
      result[1:(random_col - 1)] <- mat[start_row, 1:(random_col - 1)]
    }
    result[random_col:100] <- mat[3 - start_row, random_col:100]
    
    gampool[[i]] <- result
  }
  
  return(gampool)
}
##############################################################################
# 3) CombineGametes()
#    - Randomly pair N sperm with N eggs to form N new diploid offspring.
#    - Each offspring is a 2x100 matrix (row1 = sperm, row2 = egg).
##############################################################################
CombineGametes <- function(spermPool, eggPool) {
  N <- length(spermPool)
  newpop <- vector("list", N)
  
  # Shuffle the order of sperm and eggs to randomize pairings
  spermOrder <- sample(seq_len(N))
  eggOrder   <- sample(seq_len(N))
  
  for (i in seq_len(N)) {
    newpop[[i]] <- rbind(spermPool[[spermOrder[i]]],
                         eggPool[[eggOrder[i]]])
  }
  
  return(newpop)
}
##############################################################################
# 4) GetFit()
#    - Excludes the first locus (sex-determining) from fitness calculation.
#    - Fitness = product of colMeans for columns 2..100.
##############################################################################
GetFit <- function(pop) {
  w <- numeric(length(pop))
  
  for (i in seq_along(pop)) {
    # Exclude col 1 from fitness
    w[i] <- prod(colMeans(pop[[i]][, 2:100]))
  }
  
  return(w)
}
##############################################################################
# 5) MutPop()
#    - Prevents mutations at column 1 (the sex-determining locus).
#    - Each mutation sets one random locus in columns 2..100 to 's'.
##############################################################################
MutPop <- function(pop, num.mutes, s) {
  # pop: list of length N, each element is 2x100
  # num.mutes: number of individuals to mutate
  # s: new (mutant) value, e.g. 0.95
  
  # Randomly select which individuals to mutate
  mutes <- sample(seq_along(pop), num.mutes, replace = FALSE)
  
  for (i in mutes) {
    # Randomly pick which row (1 or 2)
    r <- sample(1:2, 1)
    # Randomly pick which column (2..100)
    c <- sample(2:100, 1)
    
    # Apply mutation
    pop[[i]][r, c] <- s
  }
  
  return(pop)
}
##############################################################################
# 6) RunSim()
#    - Creates an initial population of size N, each 2x100 matrix of all 1's.
#    - Randomly assigns sex to each row at locus1 (0 => X, 1 => Y).
#    - Splits population into malePop and femalePop each generation.
#    - Produces N sperm, N eggs, then combines them into N offspring.
#    - Mutates each new population and repeats for 'gen' generations.
#    - Repeats for 'iter' runs, returns a matrix (iter x gen) of avg fitness.
##############################################################################
RunSim <- function(s = 0.95, iter, N, gen, 
                   mincol, num.mutes) {
  
  result <- vector("list", iter)
  
  for (j in seq_len(iter)) {
    cat("Run", j, "\n")
    
    # Initialize population
    pop <- vector("list", N)
    for (i in seq_len(N)) {
      individual <- matrix(1, nrow = 2, ncol = 100)
      # Randomly assign X/Y (0/1) to each row at locus 1
      individual[, 1] <- sample(c(0, 1), 2, replace = TRUE)
      pop[[i]] <- individual
    }
    
    # Track mean fitness each generation
    mfit <- numeric(gen)
    
    for (g in seq_len(gen)) {
      cat("Generation", g, "\n")
      
      # Calculate fitness of the current population
      w <- GetFit(pop)
      
      # Separate males and females (by entire individual)
      # An individual is male if it has any '1' in its first column
      # across the two rows.
      male_idx   <- which(sapply(pop, IsMale))
      female_idx <- setdiff(seq_along(pop), male_idx)
      
      # Fitness sub-vectors
      w_male   <- w[male_idx]
      w_female <- w[female_idx]
      
      malePop   <- pop[male_idx]
      femalePop <- pop[female_idx]
      
      # If no males or no females exist, the population can't reproduce
      # (in a purely sexual scenario). This code will produce an empty
      # new population. You might handle that edge case separately.
      if (length(malePop) == 0 || length(femalePop) == 0) {
        # Optional: break or fill with NA
        cat("Warning: All individuals are one sex. Stopping.\n")
        mfit[g:gen] <- NA
        break
      }
      
      # Produce sperm (N total) from the male subpop
      spermPool <- GetMaleGametes(
        malePop = malePop,
        wMale   = w_male,
        mincol  = mincol,
        N       = N
      )
      
      # Produce eggs (N total) from the female subpop
      eggPool <- GetFemaleGametes(
        femalePop = femalePop,
        wFemale   = w_female,
        N         = N
      )
      
      # Combine sperm and eggs
      newpop <- CombineGametes(spermPool, eggPool)
      
      # Mutate new population
      pop <- MutPop(newpop, num.mutes, s)
      
      # Record mean fitness
      w_new    <- GetFit(pop)
      mfit[g]  <- mean(w_new)
    }
    
    # Store the run's fitness trace
    result[[j]] <- mfit
  }
  
  # Combine all runs into a matrix: iter rows × gen columns
  result <- do.call(rbind, result)
  foo <- list(result, pop)
  return(foo)
}
##############################################################################
# 7) GetAllYChromosomes() and GetAllXChromosomes()
#    - Same as before: gather rows from each 2x100 whose first column is 1 (Y)
#      or 0 (X) into one matrix.
##############################################################################
GetAllYChromosomes <- function(pop) {
  Ylist <- lapply(pop, function(ind) {
    idx <- which(ind[, 1] == 1)
    if (length(idx) > 0) {
      ind[idx, , drop = FALSE]
    } else {
      NULL
    }
  })
  Ylist <- Filter(Negate(is.null), Ylist)
  if (length(Ylist) == 0) {
    return(matrix(nrow = 0, ncol = 100))
  } else {
    return(do.call(rbind, Ylist))
  }
}
GetAllXChromosomes <- function(pop) {
  Xlist <- lapply(pop, function(ind) {
    idx <- which(ind[, 1] == 0)
    if (length(idx) > 0) {
      ind[idx, , drop = FALSE]
    } else {
      NULL
    }
  })
  Xlist <- Filter(Negate(is.null), Xlist)
  if (length(Xlist) == 0) {
    return(matrix(nrow = 0, ncol = 100))
  } else {
    return(do.call(rbind, Xlist))
  }
}


s = 0.95
iter <- 1
N <- 100
gen <- 100
mincol <- 25
num.mutes <- 50

popmat <- do.call(rbind,pop)
plot(colMeans(popmat))


foo <- RunSim(s = 0.95, iter = 1, N = 500, 
              gen = 5000, mincol = 25, 
              num.mutes = 20)
plot(foo[[1]][1,], type="l")
plot(colMeans(do.call(rbind,foo[[2]])[,2:100]),type="l")


bary <- GetAllYChromosomes(foo[[2]])
plot(colMeans(bary[,2:100]))
barx <- GetAllXChromosomes(foo[[2]])
plot(colMeans(barx[,2:100]))

plot(colMeans(bary[,2:100])~colMeans(barx[,2:100]))
