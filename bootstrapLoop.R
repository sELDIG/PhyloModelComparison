# This is just for internal testing - use the parallelized script for full calculations

for (seed in 3:3) {
  rmarkdown::render('Calculations.qmd',
                    params = list(seed = seed),
                    output_file = paste("bootstrap/bootstrap", seed, ".html", sep = ""))
}


