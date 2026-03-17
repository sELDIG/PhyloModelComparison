for (seed in 3:3) {
  rmarkdown::render('Analysis.qmd',
                    params = list(seed = seed),
                    output_file = paste("bootstrap/bootstrap", seed, ".html", sep = ""))
}



