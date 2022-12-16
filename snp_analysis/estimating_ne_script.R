########################################################
#     Estimating Ne on D. pseudoobscura time series    #
########################################################

##### To estimate Ne for a given data file #####
library(data.table)

### Get data from read.sync data
ne_params = function(x, y, ne_data){
  one_rep_data <- data.frame(p = eval(parse(text = paste0("ne_data@alleles$F", x, ".R", y, ".freq"))), # extract allele frequencies from data.table
                           cov = eval(parse(text = paste0("ne_data@alleles$F", x, ".R", y, ".cov"))), # extract coverages from data.table
                           pos = ne_data@alleles$pos) # extract chr positions
  colnames(one_rep_data) <- c("p", "cov", "pos")
  rownames(one_rep_data) <- NULL
  one_rep_data <- one_rep_data[order(one_rep_data$pos),]
  return(one_rep_data)
}

estimate_ne = function(dataFileName, gen, repl, N, output, pool_size, window_size, window_par){
  library(poolSeq)
  to_estNe = function(tp, N, pool_size, window_par){
  # To compute an Ne estimate per time point averaged across windows
	the_snps = sapply(tp, function(x) checkSNP(p0 = params_0$p,
                                               pt = params_t$p,
                                               cov0 = params_0$cov,
                                               covt = params_t$cov))
	params_0 = lapply(params_0, function(y) {y[the_snps[,1]]})
	params_t = lapply(params_t, function(y) {y[the_snps[,1]]})
  if (window_par == "snps"){
  # Run if fixed no. of SNPs per window
  window_max <- length(params_0$p) - (length(params_0$p) %% window_size)
  windows = sapply(seq(window_size*0.1 + 1, window_max, by = window_size), function(x) {paste(x - window_size*0.1 , x + window_size, sep=':')})
  the_estimates = sapply(windows, function(x) estimateNe(p0 = params_0$p[eval(parse(text = x))],
                                                         pt = params_t$p[eval(parse(text = x))],
                                                         cov0 = params_0$cov[eval(parse(text = x))],
                                                         covt = params_t$cov[eval(parse(text = x))],
                                                         t = tp[2] - tp[1],
                                                         method = "P.planII",
                                                         Ncensus = N,
                                                         poolSize = pool_size))

  } else if (window_par == "bp"){
  # Run if fixed kbp-long window
  window_max <- max(params_0$pos) - ((max(params_0$pos) - min(params_0$pos)) %% window_size)
  windows = sapply(seq(min(params_0$pos + window_size*0.1), window_max, by = window_size),
                   function(x) {paste(x - window_size*0.1 , x + window_size*1.1, sep=',')})
  # Check which SNP belongs to which interval
  checkSNP_interval = function(snp, windows){
    last_interval <- unlist(strsplit(windows[length(windows)], ","))
    if (snp > last_interval[2]){
      return(NA)
    } else {
    for (interval in windows){
      interval_boundaries <- unlist(strsplit(interval, ","))
      this_snp_interval <- (snp <= as.numeric(interval_boundaries[2]) & snp >= as.numeric(interval_boundaries[1]))
      if (this_snp_interval){
        snp_interval <- interval
        break
      }
    }
      return(snp_interval)
    }
    }
  # Create vector of intervals
  the_snp_windows <- sapply(params_0$pos, function(x) checkSNP_interval(snp = x,
                                                                       windows = windows))

  snps_per_window <- melt(table(the_snp_windows))
  print(paste0("Replicate ", j, ": min SNPs per window is ", min(snps_per_window$value)))
  print(paste0("Replicate ", j, ": max SNPs per window is ", max(snps_per_window$value)))
  print(paste0("Replicate ", j, ": median SNPs per window is ", median(snps_per_window$value)))

  # Add vector to params_0 and params_t
  params_0$interval <- the_snp_windows
  params_t$interval <- the_snp_windows
  # Produce estimates per interval
  the_estimates = sapply(windows, function(x) estimateNe(p0 = params_0$p[which(params_0$interval == x)],
                                                         pt = params_t$p[which(params_t$interval == x)],
                                                         cov0 = params_0$cov[which(params_0$interval == x)],
                                                         covt = params_t$cov[which(params_t$interval == x)],
                                                         t = tp[2] - tp[1],
                                                         method = "P.planII",
                                                         Ncensus = N,
                                                         poolSize = pool_size))
  } else if (window_par == "all snps") {
    the_estimates = estimateNe(p0 = params_0$p,
                             pt = params_t$p,
                             cov0 = params_0$cov,
                             covt = params_t$cov,
                             t = tp[2] - tp[1],
                             method = "P.planII",
                             Ncensus = N,
                             poolSize = pool_size)
  }

	if (output == "windows"){
      return(list(estimates = the_estimates))
    } else if (output == "mean") {
      return(mean(the_estimates, na.rm = TRUE))
    } else if (output == "median") {
      return(median(the_estimates, na.rm = TRUE))
    } else {
      return(the_estimates)
    }
  }

  ne_data = read.sync(file = dataFileName,
                      gen = gen,
                      repl = 1:repl)

  estimates = list()

  for (j in 1:repl){
    params_0 = ne_params(x = gen[j], y = j, ne_data = ne_data)
    params_t = ne_params(x = gen[j + repl], y = j, ne_data = ne_data)
    est_ne = to_estNe(tp = c(gen[j], gen[j + repl]), N = N, pool_size = pool_size, window_par = window_par)
    estimates[j] = est_ne
  }

  rm(ne_data)
  rm(params_0)
  rm(params_t)
  gc()
  return(estimates)
}

##### Arguments #####
# Run if M or P
time_series <-  c(21, 61, 114, 162, 200)
N <- 160
repl <- 4
# Run if BL
time_series <- c(26, 69, 136, 170, 202)
N <- 500
repl <- 1

pool_size <- rep(40, times = 2)
window_size <- 2000
window_par <- "snps" # "snps" or "bp" or "all snps"
output <- "median" # windows, mean, median or other

file_directory <- "/media/barata/Elements/dpseudo_proj/snp_data/intersected_sync_files"
# Run if M or P
filenames <- list.files(path = file_directory, pattern = "^R4T2_(.*)_M_snp_only_nohead.sync", full.names = TRUE)

# Run if BL
filenames <- list.files(path = file_directory, pattern = "^R1T2_(.*)_BL_snp_only_nohead.sync", full.names = TRUE)

for (i in 1:20){
  if (i <= 4){
    gen <-  rep(c(time_series[1], time_series[2]), each = repl)
  } else if (i > 4 & i <= 8){
    gen <-  rep(c(time_series[1], time_series[5]), each = repl)
  } else if (i > 8 & i <= 12){
    gen <-  rep(c(time_series[2], time_series[3]), each = repl)
  } else if (i > 12 & i <= 16){
    gen <-  rep(c(time_series[3], time_series[4]), each = repl)
  } else if (i > 16){
    gen <-  rep(c(time_series[4], time_series[5]), each = repl)
  }
  print(filenames[i])
  chr_data <-  estimate_ne(dataFileName = filenames[i],
                    gen = gen,
                    repl = repl,
                    N = N,
                    pool_size = pool_size,
                    output = output,
                    window_size = window_size,
                    window_par = window_par)

  #if (window_par == "bp" && output == "windows"){
  if (output == "windows"){
    if (repl > 1) {
      # Check for total number of windows
      number_windows <- min(c(length(unlist(chr_data[1])),
                              length(unlist(chr_data[2])),
                              length(unlist(chr_data[3])),
                              length(unlist(chr_data[4]))))

      all_data <- data.frame(rep1 = unlist(chr_data[1])[1:number_windows],
                             rep2 = unlist(chr_data[2])[1:number_windows],
                             rep3 = unlist(chr_data[3])[1:number_windows],
                             rep4 = unlist(chr_data[4])[1:number_windows])

      new_filename <- paste0(unlist(strsplit(filenames[i], "_"))[6], "_",
                             unlist(strsplit(filenames[i], "_"))[7], "_ne_estimates_P.txt") # update if M/P
    } else {
      all_data <- data.frame(rep1 = unlist(chr_data))
      new_filename <- paste0(unlist(strsplit(filenames[i], "_"))[6], "_",
                           unlist(strsplit(filenames[i], "_"))[7], "_ne_estimates_BL.txt")
    }

    write.table(all_data, file = new_filename, sep = "\t", col.names = TRUE)
  } else if (window_par == "all snps") {
    print(chr_data)
  } else if (output == "median") {
    print(chr_data)
  } else {
    #Figure this out later
  }
}
