---
title: The bad news about non-structural lack of positivity for
  continuous exposures
toc-title: Table of contents
---

The *assignment mechanism* is key to any causal analysis. There are
three assumptions: (1) the assignment is *individualistic*, that is an
individual's assignment probability cannot depend on the values of
covariates or potential outcomes of other individuals, (2) the
assignment is *unconfounded*, the assignment mechanism cannot depend on
the potential outcomes, and *probabilistic*, there is a non-zero
probability for each treatment value for every individual.
[@imbens2015causal] This paper focuses on the final assumption,
sometimes referred to as *positivity*. Formally, we denote the
continuous exposure as $T$, which takes on values in a set
$\mathscr{T}$. For individual $i$ and each value of the exposure, $t$,
there is a potential outcome $Y_i(t)$. [@imbens2000role] Each individual
also has a vector of pre-exposure covariates, $X_i$. In the most general
case, the assignment mechanism is defined as
$P(T=t|\mathbf{X}, \mathbf{Y}(t))$ (with the unconfounded assumption
this reduced to $P(T=t|\mathbf{X})$. An assignment mechanism is
*probabilisitic* if the probability of assignment to each exposure level
is strictly between 0 and 1:

$$0 < P_i(X, Y(t)) < 1, \textrm{ for each possible } X, Y(t)$$

for all $i=1, \dots, N$ In other words, the probabilistic (positivity)
assumption states that for all individuals within strata of $X$ there
exists a non-zero probability of receiving every exposure level.

## GPS

A generalized propensity score is defined as:

$$e(t, x) = P(T = t | X = x)$$ Often, for continuous exposures this is
estimated by first fitting a linear regression model prediction the
exposure from a set of pre-exposure covariates. We then use the fitted
values and the model variance in a Gaussian probability density
function. [@austin2019assessing]

$$e(t, x) = f_{T|X}(t|x) = \frac{1}{\sqrt{2\pi\hat\sigma^2}}\exp\left\{-\frac{(t-X\hat\beta)^2}{2\pi\hat\sigma^2}\right\}$$
The following stabilized weight is then used, where the numerator is the
marginal density of the exposure.

$$w = \frac{f_T(t)}{f_{T|X}(t|x)} = \frac{\hat{\sigma}_{t|x}}{\hat\sigma_t}\exp\left\{\frac{(t-X\hat\beta)^2}{2\hat\sigma_{t|x}^2}-\frac{(t-\mu_t)^2}{2\hat\sigma_t^2}\right\}$$
\## Positivity

Violations (or near violations) of the probabilistic assumption can
increase both the bias and variance of causal effect estimates.
[@petersen2012diagnosing]. There are two ways this violation can arise:
structural non-positivity and random non-positivity. The former suggests
that there is a structural mechanism that makes certain levels of the
exposure impossible for a subset of individuals. The latter are random
violations that can occur in finite samples due to chance. This paper
will focus on random violations of positivity, demonstrating that
"finite" does not mean small, but rather in the case of continuous
exposures can be any sample size that is not infinite.

## Simulations

We examine 3 scenarios:

1.  A single binary confounder with a moderate effect on the exposure
    and outcome
2.  A single binary confounder with a large effect on the exposure and
    moderate effect on the outcome
3.  A single binary confounder with a large effect on the exposure and
    large effect on the outcome.

We generate a binary confounder, $X$, from a Bernoulli distribution with
probability 0.1. The continuous exposure, $T$ is generated as follows:

$$T = aX + \varepsilon_{t|x}$$ Where $\varepsilon_{t|x}\sim N(0,1)$ and
$Y$ such that the "true" effect of $T$ is 0:

$$Y = bX + \varepsilon_{y|x}$$

where $\varepsilon_{y|x}\sim N(0,1)$. We examine $a = 1, 10$ and
$b=1,10$ for moderate and large effects, respectively. We vary the
sample size from 100 to 1,000,000, examining 100 replicates of each.

::: cell
``` {.r .cell-code}
library(tidyverse)

s <- function(n, a = 1, b = 1) {
  c <- rbinom(n, 1, 0.1)
  x <- a * c + rnorm(n)
  y <- b * c + rnorm(n)
  
  denominator_model <- lm(
    x ~ c
  )
  
  weight <- dnorm(x, mean(x), sd(x)) /
    dnorm(x, fitted(denominator_model), sigma(denominator_model))
  tibble(
    n = c(n, n, n),
    bias = c(coef(lm(y ~ x))[2],
             coef(lm(y ~ x, weights = weight))[2],
             coef(lm(y ~ x + c))[2]),
    coverage = c(1*(confint(lm(y ~ x))[2, 1] < 0 &
                      confint(lm(y ~ x))[2, 2] > 0),
                 1*(confint(lm(y ~ x, weights = weight))[2, 1] < 0 &
                      confint(lm(y ~ x, weights = weight))[2, 2] > 0),
                 1*(confint(lm(y ~ x + c))[2, 1] < 0 &
                      confint(lm(y ~ x + c))[2, 2] > 0)),
    fit = c("unadjusted", "propensity score weighted", "covariate adjustment")
  )
}
```
:::

When $a = 1$, the positivity assumption is not violated, for example
[Figure 1](#fig-a1) shows a mirrored histogram for a simulation as
specified above with a sample size of 10,000. When $a = 10$, however, we
do see a near positivity violation, for all sample sizes. Here, because
$T$ is continuous, the individuals at the population level *do* have a
non-zero probability of receiving all other values of the exposure,
however in practice it happens very infrequently. See
[Figure 2](#fig-a2), a mirrored histogram for a simulation where $a=10$
and a sample size of 10,000.

::: cell
``` {.r .cell-code}
a <- 1
b <- 1
n <- 10000
set.seed(1)

df <- tibble(
  x = rbinom(n, 1, 0.1),
  t = a * x + rnorm(n),
  y = b * x + rnorm(n)
)

library(halfmoon)
ggplot(df, aes(x = t, group = x, fill = factor(x))) + 
  geom_mirror_histogram(bins = 30) + 
  labs(fill = "X") + 
  scale_y_continuous(label = abs) +
  scale_fill_manual(values = c("cornflower blue", "orange")) + 
  theme_classic()
```

::: cell-output-display
![Figure 1: Mirrored Histogram showing overlap. a = 1, b = 1, n =
10,000](README_files/figure-markdown/fig-a1-1.png){#fig-a1}
:::
:::

::: cell
``` {.r .cell-code}
a <- 10
set.seed(1)
df <- tibble(
  x = rbinom(n, 1, 0.1),
  t = a * x + rnorm(n),
  y = b * x + rnorm(n)
)

ggplot(df, aes(x = t, group = x, fill = factor(x))) + 
  geom_mirror_histogram(bins = 30) + 
  labs(fill = "X") + 
  scale_y_continuous(label = abs) +
  scale_fill_manual(values = c("cornflower blue", "orange")) + 
  theme_classic()
```

::: cell-output-display
![Figure 2: Mirrored histogram with positivity
near-violation](README_files/figure-markdown/fig-a2-1.png){#fig-a2}
:::
:::

::: cell
``` {.r .cell-code}
n <- c(100, 500,
       seq(1000, 10000, by = 1000),
       seq(10000, 100000, by = 10000))
B <- 100
```
:::

::: cell
``` {.r .cell-code}
set.seed(1)
results <- map_df(1:B, ~map_df(n, s))
results_ex <- map_df(1:B, ~map_df(n, s, a = 10))
results_ou <- map_df(1:B, ~map_df(n, s, a = 10, b = 10))
```
:::

[Figure 3](#fig-sims) demonstrates that in the case of a continuous
exposure, when there is near violation of the positivity assumption, as
seen in panels B and C, the causal estimate using the stabilized
generalized propensity score is both biased and has large variance,
which does not appreciably decrease despite the large sample size.

::: cell
``` {.r .cell-code}
ggplot(results, aes(n, bias, color = fit, group = fit)) +
  geom_point(alpha = 0.1) +
  # geom_line(data = bias %>% 
  #             group_by(n, fit) %>% 
  #             summarise(bias = mean(bias), .groups = "drop"), color = "white") +
  geom_smooth(color = "white") + 
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  scale_color_manual(values = c("orange", "cornflower blue", "pink")) + 
  theme(legend.position = "none") -> p1

ggplot(results_ex, aes(n, bias, color = fit, group = fit)) +
  geom_point(alpha = 0.1) +
  geom_smooth(color = "white") + 
  # geom_line(data = bias_ex %>% 
  #             group_by(n, fit) %>% 
  #             summarise(bias = mean(bias), .groups = "drop"), color = "white") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  scale_color_manual(values = c("orange", "cornflower blue", "pink")) + 
  theme(legend.position = "none") -> p2

ggplot(results_ou, aes(n, bias, color = fit, group = fit)) +
  geom_point(alpha = 0.1) +
  geom_smooth(color = "white") +
  # geom_line(data = bias_ou %>% 
  #             group_by(n, fit) %>% 
  #             summarise(bias = mean(bias), .groups = "drop"), color = "white") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  scale_color_manual(values = c("orange", "cornflower blue", "pink")) + 
  theme(legend.position = "bottom") -> p3

library(patchwork)

p1 / p2 / p3 + plot_annotation(tag_levels = "A")
```

::: {.cell-output .cell-output-stderr}
    `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
    `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
    `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
:::

``` {.r .cell-code}
print(head(results))
```

::: {.cell-output .cell-output-stdout}
    # A tibble: 6 × 4
          n    bias coverage fit                      
      <dbl>   <dbl>    <dbl> <chr>                    
    1   100  0.129         1 unadjusted               
    2   100  0.106         1 propensity score weighted
    3   100  0.0665        1 covariate adjustment     
    4   500  0.0192        1 unadjusted               
    5   500 -0.129         0 propensity score weighted
    6   500 -0.0677        1 covariate adjustment     
:::

::: cell-output-display
![Figure 3: Boo. Look at all that bias and
variance](README_files/figure-markdown/fig-sims-1.png){#fig-sims}
:::
:::

[Figure 4](#fig-sims2) shows similar results for a coverage metric, with
some undercoverage even observable in Panel A.

::: cell
``` {.r .cell-code}
results2 <- results %>%
  group_by(n, fit) %>%
  summarise(mean_coverage = mean(coverage),
            mcse = sqrt(mean_coverage*(1-mean_coverage) / B))
```

::: {.cell-output .cell-output-stderr}
    `summarise()` has grouped output by 'n'. You can override using the `.groups`
    argument.
:::

``` {.r .cell-code}
results_ex2 <- results_ex %>%
  group_by(n, fit) %>%
  summarise(mean_coverage = mean(coverage),
            mcse = sqrt(mean_coverage*(1-mean_coverage) / B))
```

::: {.cell-output .cell-output-stderr}
    `summarise()` has grouped output by 'n'. You can override using the `.groups`
    argument.
:::

``` {.r .cell-code}
results_ou2 <- results_ou %>%
  group_by(n, fit) %>%
  summarise(mean_coverage = mean(coverage),
            mcse = sqrt(mean_coverage*(1-mean_coverage) / B))
```

::: {.cell-output .cell-output-stderr}
    `summarise()` has grouped output by 'n'. You can override using the `.groups`
    argument.
:::

``` {.r .cell-code}
ggplot(results2, aes(n, mean_coverage, color = fit, group = fit)) +
  geom_point(position = position_dodge(1000)) +
  geom_line(position = position_dodge(1000)) +
  geom_errorbar(aes(ymin = mean_coverage - 1.96*mcse,
                    ymax = mean_coverage + 1.96*mcse), width = .1,
                position = position_dodge(1000)) +
  geom_hline(yintercept = 0.95, lty = 2) +
  theme_classic() +
  scale_color_manual(values = c("orange", "cornflower blue", "pink")) + 
  theme(legend.position = "none") -> p4

ggplot(results_ex2, aes(n, mean_coverage, color = fit, group = fit)) +
  geom_point(position = position_dodge(1000)) +
  geom_line(position = position_dodge(1000)) +
  geom_errorbar(aes(ymin = mean_coverage - 1.96*mcse,
                    ymax = mean_coverage + 1.96*mcse), width = .1,
                position = position_dodge(1000)) +
  geom_hline(yintercept = 0.95, lty = 2) +
  theme_classic() +
  scale_color_manual(values = c("orange", "cornflower blue", "pink")) + 
  theme(legend.position = "none") -> p5

ggplot(results_ou2, aes(n, mean_coverage, color = fit, group = fit)) +
  geom_point(position = position_dodge(1000)) +
  geom_line(position = position_dodge(1000)) +
  geom_errorbar(aes(ymin = mean_coverage - 1.96*mcse,
                    ymax = mean_coverage + 1.96*mcse), width = .1,
                position = position_dodge(1000)) +
  geom_hline(yintercept = 0.95, lty = 2) +
  theme_classic() +
  scale_color_manual(values = c("orange", "cornflower blue", "pink")) + 
  theme(legend.position = "none") -> p6

p4 / p5 / p6 + plot_annotation(tag_levels = "A")
```

::: {.cell-output .cell-output-stderr}
    Warning: `position_dodge()` requires non-overlapping x intervals
    `position_dodge()` requires non-overlapping x intervals
    `position_dodge()` requires non-overlapping x intervals
    `position_dodge()` requires non-overlapping x intervals
    `position_dodge()` requires non-overlapping x intervals
    `position_dodge()` requires non-overlapping x intervals
:::

::: cell-output-display
![Figure 4: The coverage ain't
coveraging](README_files/figure-markdown/fig-sims2-1.png){#fig-sims2}
:::
:::
