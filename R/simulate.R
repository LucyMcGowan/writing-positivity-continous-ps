params <- function() {
  bind_rows(
    expand_grid(
      .n = c(
        100,
        500,
        seq(1000, 10000, by = 1000),
        seq(10000, 100000, by = 10000)
      ),
      .a = c(1, 10),
      .b = c(1, 10),
      .p = 0.5,
      .id = 1:100
    ),
    expand_grid(
      .n = 10000,
      .a = seq(2, 10, by = 2),
      .b = 1,
      .p = seq(0.1, 0.5, by = 0.1),
      .id = 1:1000
    )
  ) |>
    distinct()
}

sum_ps <- function(x, denominator_model) {
  sum(dnorm(x, fitted(denominator_model), sigma(denominator_model)))
}

simulate_data <- function(.n, .a = 1, .b = 1, .p = 0.5, .id = 1) {
  .c <- rbinom(.n, 1, .p)
  .x <- .a * .c + rnorm(.n)
  .y <- .b * .c + rnorm(.n)

  arrow_table(
    y = .y,
    x = .x,
    c = .c,
    n = .n,
    a = .a,
    b = .b,
    p = .p,
    id = .id
  )
}

fit_models <- function(data_path, .n, .a, .b, .p, .id) {
  .df <- open_dataset(data_path) |>
    filter(n == .n, a == .a, b == .b, p == .p, id == .id) |>
    collect()

  denominator_model <- lm(
    x ~ c,
    data = .df
  )

  ps <- dnorm(
    .df$x,
    fitted(denominator_model),
    sigma(denominator_model)
  )

  weight <- dnorm(.df$x, mean(.df$x), sd(.df$x)) / ps

  # overlap <- map_dbl(x, f, denominator_model)^(-1) / ps

  unadjusted <- lm(y ~ x, data = .df)
  ate_w <- lm(y ~ x, weights = weight, data = .df)
  # ato_w <- lm(y ~ x, weights = overlap_w)
  gform <- lm(y ~ x + c, data = .df)
  tibble(
    bias = c(
      coef(unadjusted)[2],
      coef(ate_w)[2],
      # coef(ato_w)[2],
      coef(gform)[2]
    ),
    variance = c(
      summary(unadjusted)$coefficients[2, 2],
      summary(ate_w)$coefficients[2, 2],
      # summary(ato_w)$coefficients[2, 2],
      summary(gform)$coefficients[2, 2]
    ),
    coverage = c(
      1 * (confint(unadjusted)[2, 1] < 0 &
        confint(unadjusted)[2, 2] > 0),
      1 * (confint(ate_w)[2, 1] < 0 &
        confint(ate_w)[2, 2] > 0),
      1 * (confint(gform)[2, 1] < 0 &
        confint(gform)[2, 2] > 0)
    ),
    fit = c(
      "unadjusted",
      "propensity score weighted (ATE)",
      # "propensity score weighted (Overlap)",
      "covariate adjustment"
    ),
    n = unique(.df$n),
    a = unique(.df$a),
    b = unique(.df$b),
    p = unique(.df$p),
    id = unique(.df$id)
  )
}

fit_simulated_models <- function(data_path, parameters) {
  parameters |>
    future_pmap_dfr(
      fit_models,
      data_path = data_path
    )
}

simulate_arrow_dataset <- function(parameters) {
  tbls <- parameters |>
    pmap(simulate_data)

  rlang::inject(concat_tables(!!!tbls)) |>
    write_dataset("data/", partitioning = c("n", "a", "b", "p"))

  "data/"
}

open_dataset <- function(data_path) {
  arrow::open_dataset(
    data_path,
    schema = arrow::schema(
      x = double(),
      y = double(),
      c = arrow::int32(),
      id = arrow::int32(),
      n = arrow::int32(),
      a = arrow::int32(),
      b = arrow::int32(),
      p = double()
    )
  )
}
