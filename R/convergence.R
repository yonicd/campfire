#' @importFrom tibble as_tibble rownames_to_column tibble
#' @importFrom rstan monitor
compute_window_convergence <- function(usamples, window_size, target_rhat, target_ess,verbose = TRUE) {

  end    <- dim(usamples)[1]
  starts <- seq(1, end, by = window_size)
  rhats  <- c()
  ess    <- c()

  for(start in starts) {

    res <- rstan::monitor(usamples[start:end,,], warmup = 0, print = FALSE)
    mdf <- as.data.frame(mdf)
    mdf <- tibble::rownames_to_column(mdf)
    mdf <- tibble::as_tibble(mdf)
    mdf <- tibble::rowid_to_column(mdf,var = 'r')
    mdf <- mdf[,order(mdf$Rhat,decreasing = TRUE)]

    rhats <- c(rhats, as.matrix(mdf[1, "Rhat"])[1])

    mdf <- mdf[,order(mdf$Bulk_ESS)]

    ess <- c(ess, as.matrix(mdf[1, "Bulk_ESS"])[1])
  }

  meets_target <- (rhats < target_rhat) & (ess > target_ess)

  converged <- any(meets_target)

  if(converged) {

    options <- which(meets_target)
    choice <- options[order(ess[options], decreasing = TRUE)[1]]

  } else {

    choice <- order(ess, decreasing = TRUE)[1]

  }

  if(verbose){

    cat("Adaptive warmup debug info (figure out how many initial windows to drop):\n\n")

    res <- tibble::tibble(max_Rhat = rhats,
                          min_Bulk_ESS = ess,
                          drop_first_n_windows = 0:(length(rhats) - 1),
                          picked = ifelse(1:length(rhats) == choice, "picked", ""))

    cols <- c('drop_first_n_windows','picked')

    res <- res[,c(cols,setdiff(names(res),cols))]

    print(res)
  }


  ret <- list(rhat      = rhats[choice],
              ess       = ess[choice],
              start     = starts[choice],
              converged = converged)

  invisible(ret)

}
