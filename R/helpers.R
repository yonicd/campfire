extras_ <- c('lp__', 'accept_stat__', 'stepsize__', 'treedepth__', 'n_leapfrog__', 'divergent__', 'energy__')

env <- new.env()

#' @importFrom readr read_csv cols col_double
read_files <- function(fit,cache = TRUE){

  if(cache|isFALSE(exists('diag_files',envir = env))){

    res <- lapply(fit$diagnostic_files(), function(file) {
      readr::read_csv(file, comment = "#", col_types = readr::cols(.default = readr::col_double()))
    })

    assign(x = 'diag_files',value = res,envir = env)
  }

  get(x = 'diag_files',envir = env)
}


getUnconstrainedSamples <- function(fit) {

  usamples_list <- lapply(read_files(fit), function(res) {

    extras_grep <- grep('^(p|g)_',names(res),value = TRUE)

    as.matrix(res[,-c(extras_,extras_grep)])

  })

  usamples <- array(0, dim = c(nrow(usamples_list[[1]]),
                               length(usamples_list),
                               ncol(usamples_list[[1]])))

  for(i in 1:length(usamples_list)) {
    usamples[, i,] <- usamples_list[[i]]
  }

  return(usamples)
}

getExtras <- function(fit) {
  lapply(read_files(fit), function(res) {
    res[,extras_]
  })
}

getStepsizes <- function(fit) {
  lapply(read_files(fit), function(res) {

    res[['stepsize__']][nrow(res)]

  })
}

#' @importFrom rstan constrain_pars stan_rdump
getInitFile <- function(stan_fit, ldraw) {

  init      <-  rstan::constrain_pars(stan_fit, as.matrix(ldraw))
  init_file <- tempfile("init", fileext = ".dat")

  rstan::stan_rdump(names(init), init_file, env = list2env(init))

  return(init_file)
}


