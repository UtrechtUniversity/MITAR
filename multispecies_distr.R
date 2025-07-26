#### Visualise density of multiple truncated distributions ####
library(TruncatedNormal)
library(ggplot2)
library(ggrepel)


##### Settings  #####
distr_means_intra <- c(-1, 0)
distr_sd_intra <- 0.9
ub_intra <- 0
distr_means_inter <- c(-1, 0, 0.5)
distr_sd_inter <- 0.5
ub_inter <- Inf
n_species <- 3
n_vals <- 2^23 # is rounded up to power of 2, see ?density

n_sample_intra <- n_species
n_sample_inter <- n_species * (n_species - 1)
DateTimeStamp <- format(Sys.time(), format = "%Y_%m_%d_%H_%M")
set.seed(seed = 314, kind = "default", normal.kind = "default", sample.kind = "default")


#### Intraspecies ####

##### Calculations #####
# Notes:
# - In the main script the following code is used to obtain the intraspecies
#   interaction coefficients:
#   rtnorm(n = nspecies, mu = selfintmean, sd = selfintsd, lb = -Inf, ub = 0,
#          method = "invtransfo")
# To do:
# - When calculating normalized y normalize by part not truncated? E.g., look
#   for an y for which TruncatedNormal::norminvp(y, -Inf, Inf) gives ...
vals_dens_ext <- vector(mode = "list", length = length(distr_means_intra))
vals_sample_ext <- vector(mode = "list", length = length(distr_means_intra))
dens_obj_ext <- vector(mode = "list", length = length(distr_means_intra))
for(i in seq_along(distr_means_intra)) {
  vals_distr <- rtnorm(n = n_vals, mu = distr_means_intra[i], sd = distr_sd_intra,
                       lb = -Inf, ub = ub_intra)
  if(is.finite(ub_intra)) {
    dens_obj <- density(vals_distr, to = ub_intra)
  } else {
    dens_obj <- density(vals_distr)
  }
  dens_obj_ext[[i]] <- dens_obj
  vals_dens <- data.frame(x = dens_obj$x, y = dens_obj$y, drawn = FALSE,
                          order = "", set = i)
  vals_dens$ynorm <- vals_dens$y / sum(vals_dens$y)
  sample_x <- sample(x = vals_distr, size = n_sample_intra, replace = FALSE)
  sample_y <- numeric(length = length(sample_x))
  for(j in seq_along(sample_x)) {
    index_draw <- which.min(abs(sample_x[j] - vals_dens$x))
    sample_y[j] <- vals_dens$y[index_draw]
    vals_dens[index_draw, "order"] <- as.character(j)
    vals_dens[index_draw, "drawn"] <- TRUE
  }
  vals_dens_ext[[i]] <- vals_dens
  vals_sample_ext[[i]] <- data.frame(x = sample_x, y = sample_y, set = i)
}
vals_dens$set <- as.factor(vals_dens$set)
vals_dens_ext <- do.call(rbind, vals_dens_ext)
vals_dens_ext$set <- as.factor(vals_dens_ext$set)
vals_sample_ext <- do.call(rbind, vals_sample_ext)
vals_sample_ext$set <- as.factor(vals_sample_ext$set)

##### Plots #####
# This part can also be plotted in base-R with
# plot(dens_obj_ext[[1]], col = 1, ylim = c(0, max(vals_dens_ext$y)));
# for(index_plot in 2L:length(dens_obj_ext)) {
#   lines(dens_obj_ext[[index_plot]], col = index_plot)
# }; grid() or with
# plot(x = vals_dens_ext$x, y = vals_dens_ext$y, col = vals_dens_ext$set,
#      pch = 16, cex = 0.5, type = "p", main = "Density plot", xlab = "Value",
#      ylab = "Density"); grid(); abline(h = 0, col = "darkgrey")
# The vertical lines can be added using:
# abline(v = vals_sample_ext$x, col = "darkgrey", lwd = 2, lty = 2)
labset <- paste0("mean: ", signif(distr_means_intra, 3),
                 ", sd: ", signif(distr_sd_intra, 3))
names(labset) <- seq_along(distr_means_intra)
my_labeller_dens <- labeller(set = labset, .default = label_value)

# To have same colours for the same mean values of the intraspecies and
# interspecies interaction coefficients
n_col_max <- max(length(distr_means_intra), length(distr_means_inter))
if(n_col_max != 3L) {
  stop("Adjust ordering of palette colours to the changed number of interactions.")
}
# my_cols <- ((scales::hue_pal()(n_col_max))[c(1, 3, 2)])
my_cols <- c("#e30613", "#82368c", "#009640")
my_cols_intra <- my_cols[seq_along(distr_means_intra)]
names(my_cols_intra) <- levels(vals_dens_ext$set)

# Notes:
# - Does geom_text_repel() change the seed of the random number generator?
# To do:
# - Put all data in geom_text_repel, and add condition on which values should be
#   labelled, so that text also evades the curve:
#   geom_label_repel(mapping = aes(label = ifelse(order != "", signif(x, 2), ""),
#                                  x = x, y = 0.5 * y, col = set),
#                    data = vals_dens_ext, size = 3, label.size = NA,
#                    min.segment.length = 0, max.overlaps = 100,
#                    show.legend = FALSE, inherit.aes = FALSE)
p_intra <- ggplot(data = vals_dens_ext,
                  aes(x = x, y = y, col = set)) +
  theme_bw(base_size = 21) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        plot.tag.position = c(0.0125, 0.9875)) +
  scale_colour_manual(values = my_cols_intra) +
  geom_line(size = 0.75) +
  coord_cartesian(xlim = c(-3, 0), expand = FALSE) +
  labs(x = "Intraspecies interaction coefficient",
       y = "Probability density (arbitrary units)",
       tag = "A")
p_intra
ggsave(paste0("distr_intra", DateTimeStamp, ".png"),
       width = 1650, height = 2675, units = "px", dpi = 300)
ggsave(paste0(DateTimeStamp, "FigS01A.png"),
       width = 1650, height = 2675, units = "px", dpi = 300)

p_intra_lines <- p_intra +
  facet_wrap(facets = vars(set), ncol = 1, labeller = my_labeller_dens) +
  geom_linerange(mapping = aes(x = x, ymin = -Inf, ymax = Inf, col = set),
                 data = vals_sample_ext, size = 0.75, linetype = 2,
                 show.legend = FALSE, inherit.aes = FALSE)

p_intra_lines +
  geom_label_repel(mapping = aes(label = ifelse(order != "", order, ""),
                                 x = x, y = y, col = set),
                   data = vals_dens_ext, size = 4, box.padding = 0.05,
                   label.padding = 0.15, label.size = 0.15,
                   min.segment.length = 100, max.overlaps = 100,
                   show.legend = FALSE, direction = "y", inherit.aes = FALSE)
ggsave(paste0("distr_intra_numbered", DateTimeStamp, ".png"),
       width = 1650, height = 2675, units = "px", dpi = 300)

p_intra_lines +
  geom_label_repel(mapping = aes(label = signif(x, 3), x = x, y = 0.75 * y,
                                 col = set),
                   data = vals_sample_ext, size = 4, box.padding = 0.05,
                   label.padding = 0.15, label.size = 0.15,
                   min.segment.length = 1e100, show.legend = FALSE,
                   direction = "y", inherit.aes = FALSE)
ggsave(paste0("distr_intra_labelled", DateTimeStamp, ".png"),
       width = 1650, height = 2675, units = "px", dpi = 300)


#### Interspecies ####

##### Calculations #####
# Notes:
# - A faster way to obtain the first columns of vals_dens_ext would be
#   n_out <- 301
#   p_x <- rep(seq(from = -3, to = 2, length.out = n_out),
#              each = length(distr_means_inter))
#   p_y <- dnorm(x = p_x, mean = distr_means_inter, sd = distr_sd_inter)
#   p_mean <- rep(distr_means_inter, n_out)
#   p_data <- data.frame(x = p_x, y = p_y, set = p_mean)
vals_dens_ext <- vector(mode = "list", length = length(distr_means_inter))
vals_sample_ext <- vector(mode = "list", length = length(distr_means_inter))
dens_obj_ext <- vector(mode = "list", length = length(distr_means_inter))
for(i in seq_along(distr_means_inter)) {
  vals_distr <- rtnorm(n = n_vals, mu = distr_means_inter[i], sd = distr_sd_inter,
                       lb = -Inf, ub = ub_inter)
  if(is.finite(ub_inter)) {
    dens_obj <- density(vals_distr, to = ub_inter)
  } else {
    dens_obj <- density(vals_distr)
  }
  dens_obj_ext[[i]] <- dens_obj
  vals_dens <- data.frame(x = dens_obj$x, y = dens_obj$y, drawn = FALSE, set = i)
  vals_dens$ynorm <- vals_dens$y / sum(vals_dens$y)
  sample_x <- sample(x = vals_distr, size = n_sample_inter, replace = FALSE)
  sample_y <- numeric(length = length(sample_x))
  for(j in seq_along(sample_x)) {
    index_draw <- which.min(abs(sample_x[j] - vals_dens$x))
    sample_y[j] <- vals_dens$y[index_draw]
    vals_dens[index_draw, "order"] <- as.character(j)
    vals_dens[index_draw, "drawn"] <- TRUE
  }
  vals_dens_ext[[i]] <- vals_dens
  vals_sample_ext[[i]] <- data.frame(x = sample_x, y = sample_y, set = i)
}
vals_dens$set <- as.factor(vals_dens$set)
vals_dens_ext <- do.call(rbind, vals_dens_ext)
vals_dens_ext$set <- as.factor(vals_dens_ext$set)
vals_sample_ext <- do.call(rbind, vals_sample_ext)
vals_sample_ext$set <- as.factor(vals_sample_ext$set)

##### Plots #####
# This part can also be plotted in base-R with:
# plot(dens_obj_ext[[1]], col = 1, ylim = c(0, max(vals_dens_ext$y)));
# for(index_plot in 2L:length(dens_obj_ext)) {
#   lines(dens_obj_ext[[index_plot]], col = index_plot)
# }; grid() or with
# plot(x = vals_dens_ext$x, y = vals_dens_ext$y, col = vals_dens_ext$set,
#      pch = 16, cex = 0.5, type = "p", main = "Density plot", xlab = "Value",
#      ylab = "Density"); grid(); abline(h = 0, col = "darkgrey")
# The vertical lines can be added using:
# abline(v = vals_sample_ext$x, col = "darkgrey", lwd = 2, lty = 2)
labset <- paste0("mean: ", signif(distr_means_inter, 3))
names(labset) <- seq_along(distr_means_inter)
my_labeller_dens <- labeller(set = labset, .default = label_value)

my_cols_inter <- my_cols[seq_along(distr_means_inter)]
names(my_cols_inter) <- levels(vals_dens_ext$set)


# Notes:
# - See 'Notes' above in the section 'Intraspecies'
p_inter <- ggplot(data = vals_dens_ext,
                  aes(x = x, y = y, col = set)) +
  theme_bw(base_size = 21) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        plot.tag.position = c(0.0125, 0.9875)) +
  scale_colour_manual(values = my_cols_inter) +
  geom_line(size = 0.75) +
  coord_cartesian(xlim = c(-3, 2), expand = FALSE) +
  labs(x = "Interspecies interaction coefficient",
       y = "Probability density (arbitrary units)",
       tag = "B")
p_inter
ggsave(paste0("distr_inter", DateTimeStamp, ".png"),
       width = 2 * 1650 * 5/6, height = 2675, units = "px", dpi = 300)
ggsave(paste0(DateTimeStamp, "FigS01B.png"),
       width = 2 * 1650 * 5/6, height = 2675, units = "px", dpi = 300)

p_inter_lines <- p_inter +
  facet_wrap(facets = vars(set), ncol = 1, labeller = my_labeller_dens) +
  geom_linerange(mapping = aes(x = x, ymin = -Inf, ymax = Inf, col = set),
                 data = vals_sample_ext, size = 0.75, linetype = 2,
                 show.legend = FALSE, inherit.aes = FALSE)

p_inter_lines +
  geom_label_repel(mapping = aes(label = ifelse(order != "", order, ""),
                                 x = x, y = y, col = set),
                   data = vals_dens_ext, size = 4, box.padding = 0.05,
                   label.padding = 0.15, label.size = 0.15,
                   min.segment.length = 100, max.overlaps = 100,
                   show.legend = FALSE, direction = "y", inherit.aes = FALSE)
ggsave(paste0("plot_distr_inter_numbered", DateTimeStamp, ".png"),
       width = 2675, height = 1650, units = "px", dpi = 300)

p_inter_lines +
  geom_label_repel(mapping = aes(label = signif(x, 3), x = x, y = 0.75 * y,
                                 col = set),
                   data = vals_sample_ext, size = 4, box.padding = 0.05,
                   label.padding = 0.15, label.size = 0.15,
                   min.segment.length = 1e100, show.legend = FALSE,
                   direction = "y", inherit.aes = FALSE)
ggsave(paste0("plot_distr_inter_labelled", DateTimeStamp, ".png"),
       width = 2675, height = 1650, units = "px", dpi = 300)
