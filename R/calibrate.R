#' Calibrate trial specification to yield desired value of metric
#'
#'
#'
#' @inheritParams run_trials
#' @param target double, the target value
#' @param search_range two-element vector, the lower and upper values of the range to explore. The default is `c(0.9, 1)` corresponding to explore superiority thresholds between 90% and 100%.
#' @param acq_fun string, name of acqusition function to use. One of `"ei"`, `"poi"` and `"lcb"` (default) for expected improvement, probability of improvement and lower confidence bound.
#' @param superior_only logical, should the trial end with superiority only (`TRUE`, default) or is it sufficient for the trial to be conclusive (`FALSE`)?
#' @param inferiority logical, optimise also the inferiority threshold (`TRUE`, default) as `1 - p(superiority)` or nor (`FALSE`)?
#' @param verbose logical, should the optimiser write updates to the console and visualise its progress with plots? (default: `FALSE`)
#' @param base_seed integer used for reproducible results. For optimisation, this really should be provided, to ensure that the optimiser is using equivalent results; if not, random fluctuations will hamper its exploration.
#' @param optims_controls settings for the optimiser with sensible defaults for calibrating w.r.t the type 1 error rate. See Details.
#'
#' @details
#' - `n0` the number of points where to evaluate the GP
#' - `max_iter` the maximum number of calibration iterations
#' - `grid_res` the resolution of the grid on which the acquisition function is evaluated
#' - `tol` tolerace, if |target - metric| <= tol, the optimiser halts
#' - `pad_sd` if the optimiser lands on an already explored value, a
#'     Gaussian noise of the form N(0, pad_sd) is added to move the optimiser to
#'     a new value. Most relevant with a very coarse grid and low tolerance.
#' - `kappa` tunable hyperparamter for the LCB acquisition function
#'
#' @return A list with four elements
#' - `trial_spec`: a calibrated version of the input `trial_spec`
#' - `pred_plots`: a list of plots, one for each optimisation iteration showing the predictions with standard deviations (lines + ribbon) and actually computed values (points, red = last, black = first)
#' - `acq_plots`: a list of plots, one for each optimisation iterations showing the values of the acquisition functions and the value at which to evaluate the trial (vertical dashed line)
#' - `evaluations`: a data frame with all evaluations
#'
#' @export
#'
#' @examples [Coming]
#'
calibrate <- function(
    trial_spec,
    target = 0.05,
    search_range = c(0.9, 1),
    acq_fun = "ei",
    superior_only = TRUE,
    inferiority = TRUE,
    cores = 1,
    n_rep = 100,
    progress = NULL,
    base_seed = NULL,
    verbose = FALSE,
    controls = list(n_initial = 4, n_max = 25, grid_res = 10000, tol = 0.005, pad_sd = 0.01, kappa = 2)
) {

  # Function to optimise
  f <- function(thres) {
    thres <- thres * diff(search_range) + search_range[1] # scale to natural scale
    trial_spec$superiority <- thres
    if (isTRUE(inferiority)) trial_spec$inferiority <- 1 - thres

    try(
      {sims <- run_trials(trial_spec, n_rep = n_rep, cores = cores, base_seed = base_seed, progress = progress)},
      silent = TRUE
    )
    if (!exists("sims", inherits = FALSE)) return(1 - target)

    estimates <- check_performance(sims)
    return(abs(target - estimates$est[estimates$metric == metric]))
  }

  rescale <- function(x) {
    # Put back to original scale
    x * diff(search_range) + search_range[1]
  }

  # Acquisition functions
  ## Probability of improvement -- higher is better
  poi <- function(mu, sigma) {
    ifelse(sigma == 0, 0, pnorm((y_best - mu) / sigma)) # = poi
    # 1 - poi (if maximising instead of minimising)
  }

  ## Expected improvement -- higher is better
  ei <- function(mu, sigma) {
    gamma <- ifelse(sigma == 0, 0, (y_best - mu)/sigma)
    phi <- pnorm(gamma)
    sigma * (gamma * phi + dnorm(gamma))
  }

  ## Lower confidence bound -- lower is better
  lcb <- function(mu, sigma) {
    mu - controls$kappa * sigma
  }

  # Enforce defaults in control argument if unspecified
  controls$n_initial <- controls$n_initial %||% 4
  controls$n_max <- controls$n_max %||% 25
  controls$grid_res <- controls$grid_res %||% 10000
  controls$tol <- controls$tol %||% 0.005
  controls$pad_sd <- controls$pad_sd %||% 0.01
  controls$kappa <- controls$kappa %||% 2

  # Initial values
  if (isTRUE(verbose)) message("Setting up initial grid evaluation")

  pred_df <- data.frame(x = seq(0, 1, length.out = controls$grid_res))
  metric <- if (superior_only) "prob_superior" else "prob_conclusive"

  evaluations <- data.frame(y = numeric(0), x = numeric(0))
  for (x in seq(0, 1, length.out = controls$n_initial)) {
    if (isTRUE(verbose)) message("Adding evaluation at initial grid point")
    evaluations <- rbind(evaluations, data.frame(y = f(x), x = x))
  }

  pred_plots <- list()
  acq_plots <- list()

  for (n in seq_len(controls$n_max)) {
    y_best <- min(evaluations[, "y"]) # best value so far

    if (isTRUE(verbose)) {
      message(
        "Iteration no. ", n,
        ". Best threshold undtil now: ",
        rescale(evaluations[which.min(evaluations[, "y"]), "x"])
      )
    }

    fit <- GP_fit(
      X = evaluations$x,
      Y = evaluations$y,
      corr = list(type = "exponential", power = 1.95)
    )

    df <- as.data.frame(predict.GP(fit, xnew = pred_df)$complete_data)
    colnames(df) <- c("x_grid", "mu", "mse")
    df$x_original_scale <- rescale(df$x_grid)
    df$sigma <- sqrt(df$mse)
    df$acq_value <- do.call(acq_fun, list(mu = df$mu, sigma = df$sigma))

    dir <- ifelse(acq_fun == "lcb", -1, 1)
      # invert ranking for lower conf. bound (to harmonise with POI and EI)
    next_x_idx <- which.max(rank(dir * df$acq_value, ties.method = "random"))
    next_x <- df[next_x_idx, "x_grid"]

    while (next_x %in% evaluations[, "x"]) {
      if (isTRUE(verbose)) message("next_x was ", rescale(next_x), " -- nudging")
      next_x <- max(0, min(1, next_x + rnorm(1, 0, controls$pad_sd)))
    }

    points_df <- as.data.frame(evaluations)
    points_df$colour <- seq_len(nrow(points_df))
    x_lab <- ifelse(isTRUE(superior_only), "Prob. of superiority", "Prob. of conclusiveness")

    pred_plots[[n]] <- ggplot2::ggplot(df, ggplot2::aes(x = x_original_scale)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = mu - sigma, ymax = mu + sigma), alpha = 0.2) +
      ggplot2::geom_line(aes(y = mu)) +
      ggplot2::geom_point(ggplot2::aes(x * diff(search_range) + search_range[1], y, colour = colour),
                          points_df, show.legend = FALSE) +
      ggplot2::scale_colour_gradient(low = "black", high = "red") +
      ggplot2::labs(x = x_lab, y = NULL, title = "Gaussian process-predicted values")

    acq_plots[[n]] <- ggplot2::ggplot(df, ggplot2::aes(x_original_scale, acq_value)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = rescale(next_x), linetype = 2) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Value of acquisition function")

    if (isTRUE(verbose)) {
      print(fit_plots[[n]] / acq_plots[[n]])
    }

    if (y_best <= controls$tol) break

    evaluations <- rbind(evaluations, c(f(next_x), next_x))
  }

  evaluations$x <- rescale(evaluations$x)
  best_x <- evaluations[which.min(evaluations[, "y"]), "x"]
  trial_spec$superiority <- best_x
  if (isTRUE(inferiority)) trial_spec$inferiority <- 1 - best_x

  list(
    trial_spec = trial_spec,
    pred_plots = pred_plots,
    acq_plots = acq_plots,
    evaluations = evaluations
  )
}
