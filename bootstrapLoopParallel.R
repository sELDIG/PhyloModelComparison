library(parallel)

# Set up cluster
n_cores <- 30  # Leave one core free
cl <- makeCluster(n_cores)

# Export any needed variables/packages to workers
clusterEvalQ(cl, library(rmarkdown))

# Run in parallel
parLapply(cl, 3:4, function(seed) {
  rmarkdown::render('Analysis.qmd',
                    params = list(seed = seed),
                    output_file = paste("bootstrap/bootstrap", seed, ".html", sep = ""))
})

# Stop cluster
stopCluster(cl)