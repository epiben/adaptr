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
    superior_only = TRUE,
    inferiority = TRUE,
    cores = 1,
    n_initial = 4,
    n_iter_max = 25,
    n_rep = 1000,
    grid_res = 10000,
    grid_res_incr = 10,
    grid_res_max = grid_res * 1000,
    tol = 0.0001,
    base_seed = as.numeric(Sys.Date()),
    progress = NULL,
    verbose = FALSE
) {

  # Housekeeping
  assert_pkgs(c("GPfit", "ggplot2"))

  # TODO: not sure this functionality is useful
  if (n_iter_max %% length(n_rep) != 0) {
    stop0("Ensure that n_iter_max is a multiple of n_rep if n_rep is not a scalar.")
  } else {
    n_rep <- rep(n_rep, each = n_iter_max / length(n_rep))
  }

  # Helpers
  to_original_scale <- function(x) {
    x * diff(search_range) + search_range[1]
  }

  # Expensive-to-evaluate function to approximate with GP
  f <- function(thres, n_rep = NULL) {
    tmp_trial_spec <- trial_spec
    tmp_trial_spec$superiority <- to_original_scale(thres)
    if (isTRUE(inferiority)) {
      tmp_trial_spec$inferiority <- 1 - tmp_trial_spec$superiority
    }

    try(
      sims <- run_trials(
        tmp_trial_spec,
        n_rep = n_rep,
        cores = cores,
        base_seed = base_seed,
        progress = progress
      ),
      silent = TRUE
    )
    if (!exists("sims", inherits = FALSE)) return(1)

    estimates <- check_performance(sims)
    return(estimates$est[estimates$metric == metric])
  }

  # Setup
  if (isTRUE(verbose)) {
    message("Setting up initial grid evaluation")
  }

  metric <- if (isTRUE(superior_only)) "prob_superior" else "prob_conclusive"

  evaluations <- data.frame(y = numeric(0), x = numeric(0))
  for (x in seq(0, 1, length.out = n_initial)) {
    if (isTRUE(verbose)) {
      message("Evaluating at initial grid point")
    }

    evaluations <- rbind(evaluations, data.frame(y = f(x, n_rep[1]), x = x))
  }

  surr_plots <- list() # plots of the surrogat function (with evaluated points)

  # Optimise w.r.t. chosen metric and target
  for (n in seq_len(n_iter_max)) {
    best_idx <- which.min(abs(evaluations$y - target))

    if (isTRUE(verbose)) {
      message(sprintf(
        "Iteration no. %s. Best threshold so far: %s. Best metric value so far: %s",
        n,
        to_original_scale(evaluations[best_idx, "x"]),
        round(evaluations[best_idx, "y"], 1 + ceiling(-log10(tol)))
      ))
    }

    fit <- GPfit::GP_fit(
      X = evaluations$x,
      Y = evaluations$y,
      corr = list(type = "exponential", power = 1.95)
    )

    pred_df <- data.frame(x = seq(0, 1, length.out = grid_res))
    df <- as.data.frame(GPfit::predict.GP(fit, xnew = pred_df)$complete_data)
    colnames(df) <- c("x_grid", "mu", "mse")
    df$x_original_scale <- to_original_scale(df$x_grid)
    df$sigma <- sqrt(df$mse)

    next_x <- df[which.min(abs(df$mu - target)), "x_grid"]
    next_x_natural <- to_original_scale(next_x)

    x_lab <- if (isTRUE(superior_only) "Prob. of superiority" else "Prob. of conclusiveness"

    surr_plots[[n]] <- ggplot2::ggplot(df, ggplot2::aes(x = x_original_scale)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = mu - sigma, ymax = mu + sigma), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = mu), size = 0.2) +
      ggplot2::geom_point(ggplot2::aes(to_original_scale(x), y, colour = seq_along(x)),
                          evaluations, show.legend = FALSE) +
      ggplot2::scale_colour_gradient(low = "black", high = "red") +
      ggplot2::geom_vline(xintercept = next_x_natural, linetype = 2, size = 0.2) +
      ggplot2::geom_hline(yintercept = target, linetype = 2, size = 0.2) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Approximate surrogate function") +
      ggplot2::theme_minimal()

    if (isTRUE(verbose)) {
      print(surr_plots[[n]])
    }

    if (abs(evaluations[best_idx, "y"] - target) <= tol) {
      if (isTRUE(verbose)) {
        message(sprintf(
          "=== Stopping because %s falls within %s +/- %s ===",
          round(evaluations[best_idx, "y"], ceiling(-log10(tol))),
          target,
          tol
        ))
      }
      break
    }

    if (next_x %in% evaluations$x) {
      if (grid_res >= grid_res_max){
        if (isTRUE(verbose)) {
          message(sprintf(
            "Reached maximum grid resolution without finding a metric value ",
            "within %s +/- %s (target +/- tol). Consider increasing grid_res_max",
            target,
            tol
          ))
        }
        break
      }

      grid_res <- round(min(grid_res_max, grid_res * grid_res_incr))
      if (isTRUE(verbose)) {
        message(sprintf(
          "Increasing grid resolution by a factor %s to %i",
          grid_res_incr,
          grid_res
        ))
      }
    } else {
      evaluations <- rbind(evaluations, c(f(next_x, n_rep[n]), next_x))
    }
  }

  # Update and extend trial_spec object
  evaluations$x <- to_original_scale(evaluations$x)

  trial_spec$superiority <- evaluations[best_idx, "x"]
  if (isTRUE(inferiority)) trial_spec$inferiority <- 1 - trial_spec$superiority

  trial_spec$calibration <- list(
    surrogate_plots = surr_plots,
    evaluations = evaluations
  )

  class(trial_spec) <- c("calibrated_trial_spec", class(trial_spec))

  return(trial_spec)
}
