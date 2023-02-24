plot_figure_a1 <- function(data_path) {
  data_path |>
    open_dataset() |>
    filter(a == 1, b == 1, p == .5, id == 1, n == 10000) |>
    collect() |>
    ggplot(aes(x = x, group = c, fill = factor(c))) +
    geom_mirror_histogram(bins = 30) +
    labs(fill = "confounder") +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = c("cornflower blue", "orange")) +
    theme_classic()
}

plot_figure_a2 <- function(data_path) {
  data_path |>
    open_dataset() |>
    filter(a == 10, b == 1, p == .5, id == 1, n == 10000) |>
    collect() |>
    ggplot(aes(x = x, group = c, fill = factor(c))) +
    geom_mirror_histogram(bins = 30) +
    labs(fill = "confounder") +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = c("cornflower blue", "orange")) +
    theme_classic()
}

plot_figure_weight <- function(data_path) {
  .df <- data_path |>
    open_dataset() |>
    filter(a == 10, b == 1, id == 1, p == .5, n == 10000) |>
    collect()

  den <- lm(x ~ c, data = .df)
  .df |>
    mutate(
      ps = dnorm(x, fitted(den), sigma(den)),
      weight = dnorm(x, mean(x), sd(x)) / ps
    ) |>
    ggplot(aes(x = x, group = c, fill = factor(c))) +
    geom_histogram(bins = 30, fill = "grey") +
    geom_histogram(bins = 30, aes(weight = weight), alpha = 0.5) +
    labs(fill = "confounder") +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = c("cornflower blue", "orange")) +
    theme_classic()
}

plot_figure_skew <- function(.df) {
  .df <- .df |>
    filter(a == 10, b == 1, p == .5, n == 10000)

  h_wt <- .df |>
    filter(fit == "propensity score weighted (ATE)")

  ggplot(h_wt, aes(x = bias)) +
    geom_histogram(
      binwidth = 0.005,
      aes(y = ggplot2::after_stat(density))
    ) +
    geom_density() +
    stat_function(
      fun = dnorm,
      args = list(
        mean = mean(h_wt$bias),
        sd = sd(h_wt$bias)
      ),
      lty = 2,
      color = "orange"
    ) +
    labs(x = "coefficient")
}

plot_figure_sims <- function(.df) {
  results <- .df |>
    filter(a == 1, b == 1, p == .5)

  results_ex <- .df |>
    filter(a == 10, b == 1, p == .5)

  results_ou <- .df |>
    filter(a == 10, b == 10, p == .5)

  ggplot(results, aes(n, bias, color = fit, group = fit)) +
    geom_point(alpha = 0.1) +
    geom_smooth() +
    # geom_smooth(color = "white") +
    geom_hline(yintercept = 0, lty = 2) +
    theme_classic() +
    scale_color_manual(values = c("orange", "cornflower blue", "pink")) +
    theme(legend.position = "none") -> p1

  ggplot(results_ex, aes(n, bias, color = fit, group = fit)) +
    geom_point(alpha = 0.1) +
    geom_smooth() +
    # geom_smooth(color = "white") +
    geom_hline(yintercept = 0, lty = 2) +
    theme_classic() +
    scale_color_manual(values = c("orange", "cornflower blue", "pink")) +
    theme(legend.position = "none") -> p2

  ggplot(results_ou, aes(n, bias, color = fit, group = fit)) +
    geom_point(alpha = 0.1) +
    geom_smooth() +
    # geom_smooth(color = "white") +
    geom_hline(yintercept = 0, lty = 2) +
    theme_classic() +
    scale_color_manual(values = c("orange", "cornflower blue", "pink")) +
    theme(legend.position = "bottom") -> p3

  p1 / p2 / p3 + plot_annotation(tag_levels = "A")
}


plot_figure_bias <- function(.df) {
  .df |>
    filter(
      fit == "propensity score weighted (ATE)",
      a %in% c(1, 2, 3, 4, 10), b == 1,
      p %in% seq(0.1, 0.5, by = 0.1),
      n == 10000
    ) |>
    mutate(a = factor(a)) |>
    ggplot(aes(x = p, y = bias, color = a, group = a)) +
    geom_jitter(height = 0, alpha = 0.5, color = "grey80") +
    geom_hline(color = "grey70", yintercept = 0) +
    geom_smooth() +
    scale_color_okabe_ito() +
    theme_classic()
}

plot_figure_variance <- function(.df) {
  .df |>
    filter(
      fit == "propensity score weighted (ATE)",
      a %in% seq(2, 10, by = 2), b == 1,
      p %in% seq(0.1, 0.5, by = 0.1),
      n == 10000
    ) |>
    group_by(a, b, p) |>
    summarise(sd = sd(bias), .groups = "drop") |>
    mutate(a = factor(a)) |>
    ggplot(aes(x = p, y = sd, color = a, group = a)) +
    geom_point() +
    geom_line() +
    scale_color_okabe_ito() +
    theme_classic()
}

