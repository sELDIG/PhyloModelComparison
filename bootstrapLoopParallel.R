# IMPORTANT: Calculations.qmd needs to be run indepdently of this to generate the non-bootstrapped results

library(parallel)

# Set up cluster
n_cores <- 10  # note random forest is also parallelized, so make sure to adjust cores here together with the RF in Calculations.qmd
cl <- makeCluster(n_cores)

# Export any needed variables/packages to workers
clusterEvalQ(cl, library(rmarkdown))

# Run in parallel
parLapply(cl, 1:250, function(seed) {
  rmarkdown::render('Calculations.qmd',
                    params = list(seed = seed),
                    output_file = paste("bootstrap/bootstrap", seed, ".html", sep = ""))
})

# Stop cluster
stopCluster(cl)



