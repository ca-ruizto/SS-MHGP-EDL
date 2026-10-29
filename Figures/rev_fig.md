---
title: "rev_fig"
output: 
  html_document: 
    fig_caption: true
    keep_md: true
date: "2026-09-23"
---






``` r
revision_root <- "C:/Users/AbsoluteArm/OneDrive - University of Oklahoma (1)/Metal_proj_c/SS-MHGP-EDL-revision/SS-MHGP-EDL-revision"
figure_dir <- file.path(revision_root, "Figures", "Revision Output")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metals <- c("Mo", "Nb", "Ta", "V", "W")
properties <- c(
  "log_max_strain",
  "bc_1", "bc_2", "bc_3", "bc_4", "bc_5", "bc_6"
)

required_files <- c(
  file.path(revision_root, "Results", "ss_all_ort_coefs_c_uns_revision.csv"),
  file.path(revision_root, "Results", "ss_all_ort_coefs_t_uns_revision.csv"),
  file.path(
    revision_root, "Results", "Revision Evidence", "report",
    "active_learning_integrated_variance_summary.csv"
  ),
  file.path(
    revision_root, "Results", "Revision Tables",
    "table_03_whole_curve_metrics_with_intervals.csv"
  ),
  file.path(
    revision_root, "Results", "Revision Tables",
    "table_03_regional_curve_metrics_with_intervals.csv"
  )
)

composition_data <- bind_rows(
  readr::read_csv(
    file.path(revision_root, "Results", "ss_all_ort_coefs_c_uns_revision.csv"),
    show_col_types = FALSE
  ) %>% mutate(type = "c"),
  readr::read_csv(
    file.path(revision_root, "Results", "ss_all_ort_coefs_t_uns_revision.csv"),
    show_col_types = FALSE
  ) %>% mutate(type = "t")
) %>%
  filter(!(alloy %in% metals)) %>%
  mutate(
    elements = (Mo > 0) + (Nb > 0) + (Ta > 0) + (V > 0) + (W > 0),
    binary_pair = alloy,
    element_percentage = 100 * case_when(
      W > 0  ~ W,
      V > 0  ~ V,
      Ta > 0 ~ Ta,
      Nb > 0 ~ Nb,
      TRUE   ~ NA_real_
    )
  ) %>%
  filter(
    elements == 2,
    round(pmax(Mo, Nb, Ta, V, W), 1) %in% c(0.1, 0.3, 0.5, 0.7, 0.9)
  ) %>%
  pivot_longer(
    cols = all_of(properties),
    names_to = "property",
    values_to = "val_true"
  ) %>%
  group_by(type, alloy_cnt, property) %>%
  mutate(residual = val_true - mean(val_true)) %>%
  ungroup() %>%
  group_by(type, property) %>%
  mutate(
    pooled_residual_sd =
      sqrt(sum(residual^2) / (n() - n_distinct(alloy_cnt))),
    scale_location = sqrt(abs(residual / pooled_residual_sd))
  ) %>%
  ungroup() %>%
  mutate(
    property_math = factor(
      property,
      levels = properties,
      labels = c("ln(epsilon[f])", paste0("beta[", 1:6, "]"))
    )
  )

scale_location_ymax <- 1.05 * max(composition_data$scale_location)
```





``` r
composition_compression <- composition_data %>%
  filter(type == "c") %>%
  ggplot(aes(element_percentage, scale_location)) +
  geom_point(size = 1.15, alpha = 0.58, colour = "black") +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    method.args = list(method = "REML"),
    colour = "blue",
    fill = "blue",
    alpha = 0.16,
    linewidth = 0.75
  ) +
  facet_grid(
    property_math ~ binary_pair,
    labeller = labeller(property_math = label_parsed)
  ) +
  scale_x_continuous(breaks = c(10, 50, 90), limits = c(7, 93)) +
  coord_cartesian(ylim = c(0, scale_location_ymax)) +
  labs(
    x = "Atomic percentage of the second named element",
    y = "Square root of absolute residual"
  ) +
  theme(
    strip.text = element_text(size = 11),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "res_composition_compression_paper.png"),
  composition_compression,
  width = 18, height = 13, dpi = 400, bg = "white"
)

composition_compression
```

![](rev_fig_files/figure-html/paper_composition_compression-1.png)<!-- -->






``` r
composition_tension <- composition_data %>%
  filter(type == "t") %>%
  ggplot(aes(element_percentage, scale_location)) +
  geom_point(size = 1.15, alpha = 0.58, colour = "black") +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    method.args = list(method = "REML"),
    colour = "blue",
    fill = "blue",
    alpha = 0.16,
    linewidth = 0.75
  ) +
  facet_grid(
    property_math ~ binary_pair,
    labeller = labeller(property_math = label_parsed)
  ) +
  scale_x_continuous(breaks = c(10, 50, 90), limits = c(7, 93)) +
  coord_cartesian(ylim = c(0, scale_location_ymax)) +
  labs(
    x = "Atomic percentage of the second named element",
    y = "Square root of absolute residual"
  ) +
  theme(
    strip.text = element_text(size = 11),
    legend.position = "none"
  )

ggsave(
  file.path(figure_dir, "res_composition_tension_paper.png"),
  composition_tension,
  width = 18, height = 13, dpi = 400, bg = "white"
)

composition_tension
```

![](rev_fig_files/figure-html/paper_composition_tension-1.png)<!-- -->





``` r
imspe_data <- readr::read_csv(
  file.path(
    revision_root, "Results", "Revision Evidence", "report",
    "active_learning_integrated_variance_summary.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    Strategy = factor(
      dplyr::recode(
        strategy,
        imspe = "IMSPE",
        maximin = "Maximin",
        max_variance = "Maximum variance",
        random = "Random"
      ),
      levels = c("IMSPE", "Maximin", "Maximum variance", "Random")
    ),
    Loading = factor(
      dplyr::recode(
        loading,
        combined = "Combined",
        c = "Compression",
        t = "Tension"
      ),
      levels = c("Combined", "Compression", "Tension")
    )
  )

imspe_plot <- ggplot(
  imspe_data,
  aes(
    x = step,
    y = mean_integrated_variance,
    colour = Strategy,
    group = Strategy
  )
) +
  geom_ribbon(
    data = filter(imspe_data, Strategy == "Random"),
    aes(
      x = step,
      ymin = integrated_variance_low,
      ymax = integrated_variance_high
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.30
  ) +
  geom_line(linewidth = 1.25) +
  geom_point(size = 2) +
  facet_wrap(~Loading, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 0:5) +
  scale_colour_brewer("Selection strategy", palette = "Dark2") +
  labs(
    x = "Added higher-order compositions",
    y = "Integrated standardized latent variance"
  ) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))

ggsave(
  file.path(figure_dir, "imspe_others_paper.png"),
  imspe_plot,
  width = 11, height = 4.8, dpi = 400, bg = "white"
)

imspe_plot
```

![](rev_fig_files/figure-html/paper_imspe_comparison-1.png)<!-- -->





``` r
whole_metrics <- readr::read_csv(
  file.path(
    revision_root, "Results", "Revision Tables",
    "table_03_whole_curve_metrics_with_intervals.csv"
  ),
  show_col_types = FALSE
) %>%
  transmute(
    Cohort, Loading,
    Region = "Whole modeled curve",
    Metric, Mean,
    lower = `CI low`,
    upper = `CI high`
  )

regional_metrics <- readr::read_csv(
  file.path(
    revision_root, "Results", "Revision Tables",
    "table_03_regional_curve_metrics_with_intervals.csv"
  ),
  show_col_types = FALSE
) %>%
  transmute(
    Cohort, Loading, Region, Metric, Mean,
    lower = `CI low`,
    upper = `CI high`
  )

regional_plot_data <- bind_rows(whole_metrics, regional_metrics) %>%
  filter(
    Metric %in% c(
      "ISE (GPa^2 strain)", "IAE (GPa strain)",
      "RMSE (GPa)", "MAE (GPa)", "MaxAE (GPa)", "CRPS (GPa)"
    )
  ) %>%
  mutate(
    Region = factor(
      Region,
      levels = c(
        "Whole modeled curve", "Elastic window",
        "Nonlinear window", "Near-failure window"
      )
    ),
    Metric = factor(
      Metric,
      levels = c(
        "ISE (GPa^2 strain)", "IAE (GPa strain)",
        "RMSE (GPa)", "MAE (GPa)", "MaxAE (GPa)", "CRPS (GPa)"
      )
    ),
    Loading = factor(Loading, levels = c("Compression", "Tension")),
    Cohort = factor(
      Cohort,
      levels = c("Binary interpolation test", "Higher-order test")
    )
  )

regional_plot <- ggplot(
  regional_plot_data,
  aes(x = Region, y = Mean, colour = Cohort)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.16,
    linewidth = 0.35,
    position = position_dodge(width = 0.46)
  ) +
  geom_point(
    size = 2.1,
    position = position_dodge(width = 0.46)
  ) +
  facet_grid(Metric ~ Loading, scales = "free_y") +
  scale_colour_brewer("Held-out set", palette = "Dark2") +
  labs(x = NULL, y = "Error estimate") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "Curve_preds_regions_paper.png"),
  regional_plot,
  width = 10.5, height = 12.5, dpi = 400, bg = "white"
)

regional_plot
```

![](rev_fig_files/figure-html/paper_regional_curve_errors-1.png)<!-- -->




