#' Calibrate trial specification to yield desired value of metric
#'
#'
#'
#' @inheritParams run_trials
#' @param target double, the target value
#' @param search_range two-element vector, the lower and upper values of the
#'   range to explore. The default is `c(0.9, 1)` corresponding to explore
#'   superiority thresholds between 90% and 100%.
#' @param acq_fun string, name of acqusition function to use. One of `"ei"`,
#'   `"poi"` and `"lcb"` (default) for expected improvement, probability of
#'   improvement and lower confidence bound. See Details for more information.
#' @param superior_only logical, should the trial end with superiority only
#'   (`TRUE`, default) or is it sufficient for the trial to be conclusive
#'   (`FALSE`)?
#' @param inferiority logical, optimise also the inferiority threshold (`TRUE`,
#'   default) as `1 - p(superiority)` or nor (`FALSE`)?
#' @param verbose logical, should the optimiser write updates to the console and
#'   visualise its progress with plots? (default: `FALSE`)
#' @param base_seed integer used for reproducible results. For optimisation,
#'   this really should be provided, to ensure that the optimiser is using
#'   equivalent results; if not, random fluctuations will hamper its
#'   exploration.
#' @param controls settings for the optimiser with sensible defaults for
#'   calibrating w.r.t the type 1 error rate. See Details.
#'
#' @details
#' The acquisition function uses the surrogate model to identify the next place
#' to evaluate the trial. Probability of improvement estimates the next value
#' most likely to yield an improvement but tends to get stuck in local minima,
#' while expected improvement strikes a good balance between exploring uncertain
#' regions and exploiting promising regions already found (the same applies to
#' the lower confidence bound method).
#'
#' The `controls` argument must be a list specifying any of the following values
#' (those unspecified will use the default values):
#' - `n_initial` the number of initial points where to evaluate the surrogate
#' function
#' - `max_iter` the maximum number of calibration iterations
#' - `grid_res` the resolution of the grid on which the acquisition function is
#' evaluated
#' - `tol` tolerace, if `abs(target - metric <= tol`, the optimiser halts
#' - `nudge_sd` if the optimiser lands on an already seen value, a Gaussian
#' noise of the form N(0, nudge_sd) is added to move the optimiser to a new
#' value. Most relevant with a very coarse grid and low tolerance.
#' - `kappa` tunable hyperparamter for the LCB acquisition function
#'
#' @return An object of class `calibrated_trial_spec` (inheriting from
#'   `trial_spec`) with updated superiority (and, if chosen, inferiority)
#'   thresholds. The object also has an additional elements named `calibration`,
#'   itself a three-element list:
#' - `surrogate_plots`: a list of plots, one for each optimisation iteration,
#' with the surrogate function and its standard deviations (lines + ribbon),
#' actually computed values (points, red = last, black = first) and the location
#' of the next evaluation (vertical dashed line)
#' - `acquisition_plots`: a list of plots, one for each optimisation iteration,
#' with the acquisition function and the locations of the next evaluation
#' evaluate the trial (vertical dashed line)
#' - `evaluations`: a data frame with all evaluations
#'
#' @export
#'
#' @examples [Coming]
#'
calibrate <- function(
    trial_spec,
    target = 0.05,
    search_range = c(0.5, 1),
    acq_fun = "lcb",
    superior_only = TRUE,
    inferiority = TRUE,
    cores = 1,
    n_rep = 1000,
    progress = NULL,
    base_seed = NULL,
    verbose = FALSE,
    controls = list(
      n_initial = 4,
      n_max = 25,
      grid_res = 10000,
      tol = 0.001,
      nudge_sd = 0.01,
      kappa = 2
    )
) {

  assert_pkgs("GPfit")

  if (!is.list(controls)) {
    warning0("The controls argument must be a list, falling back to default")
    controls <- formals()$controls
  }
  # Enforce defaults in control argument if unspecified
  controls$n_initial <- controls$n_initial %||% 4
  controls$n_max <- controls$n_max %||% 25
  controls$grid_res <- controls$grid_res %||% 10000
  controls$tol <- controls$tol %||% 0.001
  controls$nudge_sd <- controls$nudge_sd %||% 0.01
  controls$kappa <- controls$kappa %||% 2

  # Expensive-to-evaluate function to approximate with GP
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

  # Helper functions
  to_original_scale <- function(x, n_digits = ceiling(log10(controls$grid_res))) {
    # Put back to original scale, round according to grid resolution
    round(x * diff(search_range) + search_range[1], n_digits)
  }

  # Probability of improvement -- higher is better
  poi <- function(mu, sigma) {
    ifelse(sigma == 0, 0, pnorm((y_best - mu) / sigma)) # = poi
    # 1 - poi (if maximising instead of minimising)
  }

  # Expected improvement -- higher is better
  ei <- function(mu, sigma) {
    gamma <- ifelse(sigma == 0, 0, (y_best - mu) / sigma)
    phi <- pnorm(gamma)
    sigma * (gamma * phi + dnorm(gamma))
  }

  # Lower confidence bound -- lower is better
  lcb <- function(mu, sigma) {
    mu - controls$kappa * sigma
  }

  # Setup
  if (isTRUE(verbose)) message("Setting up initial grid evaluation")

  pred_df <- data.frame(x = seq(0, 1, length.out = controls$grid_res))
  metric <- if (superior_only) "prob_superior" else "prob_conclusive"

  evaluations <- data.frame(y = numeric(0), x = numeric(0))
  for (x in seq(0, 1, length.out = controls$n_initial)) {
    if (isTRUE(verbose)) message("Evaluating at initial grid point")
    evaluations <- rbind(evaluations, data.frame(y = f(x), x = x))
  }

  surr_plots <- list()
  acq_plots <- list()

  # Optimise w.r.t. chosen metric and target
  for (n in seq_len(controls$n_max)) {
    y_best <- min(evaluations$y) # best value so far

    if (isTRUE(verbose)) {
      message(
        "Iteration no. ", n, ". Best threshold undtil now: ",
        to_original_scale(evaluations[which.min(evaluations[, "y"]), "x"])
      )
    }

    fit <- GPfit::GP_fit(
      X = evaluations$x,
      Y = evaluations$y,
      corr = list(type = "exponential", power = 1.95)
    )

    df <- as.data.frame(GPfit::predict.GP(fit, xnew = pred_df)$complete_data)
    colnames(df) <- c("x_grid", "mu", "mse")
    df$x_original_scale <- to_original_scale(df$x_grid)
    df$sigma <- sqrt(df$mse)
    df$acq_value <- do.call(acq_fun, list(mu = df$mu, sigma = df$sigma))

    dir <- ifelse(acq_fun == "lcb", -1, 1)
      # invert ranking for lower conf. bound (to harmonise with POI and EI)
    next_x_idx <- which.max(rank(dir * df$acq_value, ties.method = "random"))
    next_x <- df[next_x_idx, "x_grid"]

    while (next_x %in% evaluations$x) {
      next_x <- rbeta(1, next_x * nrow(df), (1 - next_x) * nrow(df))
        # beta distribution centered at next_x with mean/sample-size parameterisation
    }

    points_df <- as.data.frame(evaluations)
    points_df$colour <- seq_len(nrow(points_df))
    x_lab <- ifelse(isTRUE(superior_only), "Prob. of superiority", "Prob. of conclusiveness")

    surr_plots[[n]] <- ggplot2::ggplot(df, ggplot2::aes(x = x_original_scale)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = mu - sigma, ymax = mu + sigma), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = mu), size = 0.2) +
      ggplot2::geom_point(ggplot2::aes(to_original_scale(x), y, colour = colour),
                          points_df, show.legend = FALSE) +
      ggplot2::scale_colour_gradient(low = "black", high = "red") +
      ggplot2::geom_vline(xintercept = to_original_scale(next_x), linetype = 2, size = 0.2) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Gaussian process-predicted values")

    acq_plots[[n]] <- ggplot2::ggplot(df, ggplot2::aes(x_original_scale, acq_value)) +
      ggplot2::geom_line(size = 0.2) +
      ggplot2::geom_vline(xintercept = to_original_scale(next_x), linetype = 2, size = 0.2) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Value of acquisition function")

    if (isTRUE(verbose)) {
      print(surr_plots[[n]])
    }

    if (y_best <= controls$tol) break

    evaluations <- rbind(evaluations, c(f(next_x), next_x))
  }

  evaluations$x <- to_original_scale(evaluations$x)
  best_x <- evaluations[which.min(evaluations[, "y"]), "x"]

  # Update and extend trial_spec object
  trial_spec$superiority <- best_x
  if (isTRUE(inferiority)) trial_spec$inferiority <- 1 - best_x

  trial_spec$calibration <- list(
    surrogate_plots = surr_plots,
    acquisition_plots = acq_plots,
    evaluations = evaluations
  )

  class(trial_spec) <- c("calibrated_trial_spec", class(trial_spec))

  return(trial_spec)
}
