ConleySEs <- function(reg,
                      unit,
                      time,
                      lat,
                      lon,
                      kernel = "bartlett",
                      dist_fn = "Haversine",
                      dist_cutoff = 500,
                      lag_cutoff = 5,
                      lat_scale = 111,
                      verbose = FALSE,
                      cores = 1,
                      balanced_pnl = FALSE) {
  Fac2Num <- function(x) {
    as.numeric(as.character(x))
  }

  if (cores > 1 & !requireNamespace("parallel", quietly = TRUE)) {
      stop("parallel pkg needed for this function to work with multiple cores. Please install it.",
           call. = FALSE)
  }

  if (class(reg) == "felm") {
    Xvars <- rownames(reg$coefficients)
    dt = data.table::data.table(
      reg$cY,
      reg$cX,
      fe1 = Fac2Num(reg$fe[[1]]),
      fe2 = Fac2Num(reg$fe[[2]]),
      coord1 = Fac2Num(reg$clustervar[[1]]),
      coord2 = Fac2Num(reg$clustervar[[2]])
    )
    data.table::setnames(dt,
                         c("fe1", "fe2", "coord1", "coord2"),
                         c(names(reg$fe)[1:2], names(reg$clustervar)))
    dt = dt[, e := as.numeric(reg$residuals)]

  } else {
    message("Model class not recognized.")
    break
  }

  n <- nrow(dt)
  k <- length(Xvars)

  # Renaming variables:
  orig_names <- c(unit, time, lat, lon)
  new_names <- c("unit", "time", "lat", "lon")
  data.table::setnames(dt, orig_names, new_names)

  # Empty Matrix:
  XeeX <- matrix(nrow = k, ncol = k, 0)

  #================================================================
  # Correct for spatial correlation:
  timeUnique <- unique(dt[, time])
  Ntime <- length(timeUnique)
  data.table::setkey(dt, time)

  if (verbose) {
    message("Starting to loop over time periods...")
  }

  if (balanced_pnl) {
    sub_dt <- dt[time == timeUnique[1]]
    lat <- sub_dt[, lat]
    lon <- sub_dt[, lon]
    rm(sub_dt)

    if (balanced_pnl &
        verbose) {
      message("Computing Distance Matrix...")
    }

    d <- DistMat(cbind(lat, lon), cutoff = dist_cutoff, kernel, dist_fn)
    rm(list = c("lat", "lon"))
  }


  if (requireNamespace("parallel", quietly = TRUE) & cores > 1) {
    XeeXhs <- parallel::mclapply(timeUnique,function(t) {
      iterateObs(dt, Xvars,
                 sub_index = t,
                 type = "spatial",
                 cutoff = dist_cutoff,
                 balanced_pnl = balanced_pnl, d=d,
                 kernel = kernel,
                 dist_fn = dist_fn,
                 verbose = verbose)},
      mc.cores = cores)
  } else {
    XeeXhs <- lapply(timeUnique, function(t) {
      iterateObs(dt, Xvars,
                 sub_index = t,
                 type = "spatial",
                 cutoff = dist_cutoff,
                 balanced_pnl = balanced_pnl, d=d,
                 kernel = kernel,
                 dist_fn = dist_fn,
                 verbose = verbose)
    })
  }


  if(balanced_pnl){ rm(d) }
  # First Reduce:
  XeeX <- Reduce("+",  XeeXhs)

  X <- as.matrix(dt[, eval(Xvars), with = FALSE])
  invXX <- solve(t(X) %*% X) * n

  #================================================================
  # Correct for serial correlation:
  panelUnique <- unique(dt[, unit])
  Npanel <- length(panelUnique)
  data.table::setkey(dt, unit)

  if (verbose) {
    message("Starting to loop over units...")
  }

  if (requireNamespace("parallel", quietly = TRUE) & cores > 1) {
    XeeXhs <- parallel::mclapply(panelUnique,function(t) {
      iterateObs(dt, Xvars,
                 sub_index = t,
                 type = "serial",
                 cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl,
                 kernel = kernel,
                 dist_fn = dist_fn,
                 verbose = verbose)},
      mc.cores = cores)
  } else {
    XeeXhs <- lapply(panelUnique, function(t) {
      iterateObs(dt, Xvars,
                 sub_index = t,
                 type = "serial",
                 cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl,
                 kernel = kernel,
                 dist_fn = dist_fn,
                 verbose = verbose)
    })
  }

  XeeX_serial <- Reduce("+",  XeeXhs)

  XeeX <- XeeX + XeeX_serial

  V_spatial_HAC <- invXX %*% (XeeX / n) %*% invXX / n
  V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2
  return(V_spatial_HAC)
}

iterateObs <- function(dt, Xvars, sub_index, type, cutoff, balanced_pnl, d,
                       verbose, kernel, dist_fn) {
  k <- length(Xvars)
  if (type == "spatial" & balanced_pnl) {
    sub_dt <- dt[time == sub_index]
    n1 <- nrow(sub_dt)
    if (n1 > 1000 &
        verbose) {
      message(paste("Starting on sub index:", sub_index))
    }

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]

    XeeXhs <- Bal_XeeXhC(d, X, e, n1, k)

  } else if (type == "spatial" & !balanced_pnl) {
    sub_dt <- dt[time == sub_index]
    n1 <- nrow(sub_dt)
    if (n1 > 1000 &
        verbose) {
      message(paste("Starting on sub index:", sub_index))
    }

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]
    lat <- sub_dt[, lat]
    lon <- sub_dt[, lon]

    # If n1 >= 50k obs, then avoiding construction of distance matrix.
    # This requires more operations, but is less memory intensive.
    if (n1 < 5 * 10 ^ 4) {
      XeeXhs <- XeeXhC(cbind(lat, lon), cutoff, X, e, n1, k,
                       kernel, dist_fn)
    } else {
      XeeXhs <- XeeXhC_Lg(cbind(lat, lon), cutoff, X, e, n1, k,
                          kernel, dist_fn)
    }

  } else if (type == "serial") {
    sub_dt <- dt[unit == sub_index]
    n1 <- nrow(sub_dt)
    if (n1 > 1000 &
        verbose) {
      message(paste("Starting on sub index:", sub_index))
    }

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]
    times <- sub_dt[, time]

    XeeXhs <- TimeDist(times, cutoff, X, e, n1, k)
  }

  XeeXhs
}
