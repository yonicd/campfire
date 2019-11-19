compute_window_convergence = function(usamples, window_size, target_rhat, target_ess) {
  end = dim(usamples)[1]
  starts = seq(1, end, by = window_size)
  rhats = c()
  ess = c()

  for(start in starts) {
    mdf = as_tibble(rownames_to_column(as.data.frame(monitor(usamples[start:end,,], warmup = 0, print = FALSE)))) %>%
      mutate(r = row_number()) %>%
      select(rowname, r, everything()) %>%
      arrange(-Rhat)

    rhats = c(rhats, as.matrix(mdf[1, "Rhat"])[1])
    ess = c(ess, as.matrix((mdf %>% arrange(n_eff))[1, "Bulk_ESS"])[1])
  }

  meets_target = (rhats < target_rhat) & (ess > target_ess)

  converged = any(meets_target)
  if(converged) {
    choice = which(meets_target)[1]
  } else {
    choice = order(rhats)[1]
  }

  print("Adaptive warmup debug info (figure out how many initial windows to drop):")
  print(tibble(max_Rhat = rhats,
               min_Bulk_ESS = ess,
               drop_first_n_windows = 0:(length(rhats) - 1)) %>%
          select(drop_first_n_windows, everything()))

  return(list(rhat = rhats[choice],
              ess = ess[choice],
              start = starts[choice],
              converged = converged))
}