## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

library(ggplot2)
theme_set(theme_gray(base_size = 22))
library(dplyr)
source("src/poststats/R_rainclouds.R")
library(reshape2) # acast()
library(plotly)
library(parallel) # mclapply

INPATH <- 'out/pymvpa/'
OUTPATH <- 'out/poststats/'
SAMPLE_SIZE <- 16
TIME_STEP <- 200 # ms

plot_timeseries <- function(df) {
    best <- user_maxima(df)
    xbreaks <- seq(0, 19800, TIME_STEP)
    xlabels <- sapply(xbreaks,
                      function(t) {if (t %% 1000 == 0) as.character(t) else ""})
    ybreaks = c(seq(0, 1, .1), mean(df$mean_accuracy), 1/3)
    ylabels = c(seq(0, 1, .1), "media", "azar")
    ggplot(df, aes(x = ms,
                   y = mean_accuracy,
                   group = subject,
                   color = subject)) +
        geom_line(aes(alpha=.01), show.legend = FALSE) +
        geom_point(aes(x = ms, y = mean_accuracy, color = subject), best) +
        scale_x_continuous(breaks = xbreaks, minor_breaks = NULL,
                           labels = xlabels) +
        scale_y_continuous(breaks = ybreaks, minor_breaks = NULL,
                           labels = ylabels) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x="Latencia estímulo-respuesta (ms)",
             y="Exactitud de clasificación") +
        guides(colour = FALSE)
}

plot_timeseries_2 <- function(df2, order) {
    xbreaks <- seq(0, 19800, TIME_STEP)
    xlabels <- sapply(xbreaks,
                      function(t) {if (t %% 1000 == 0) as.character(t) else ""})
    ggplot(df2, aes(x = ms, y = factor(subject, levels = order),
                    z = mean_accuracy)) +
        geom_tile(aes(fill = mean_accuracy)) +
        ## scale_x_discrete(breaks = as.character(xbreaks), labels = xlabels) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "Latencia estímulo-respuesta (ms)", y = "Sujeto") +
        scale_fill_gradient2(name = "Exactitud de clasificación",
                             low = "red", mid = "white",
                             high = "blue", midpoint = mean(df2$mean_accuracy))
}

plot_mean_timeseries_denoise <- function(df, radii = 0) {
    df2 <- do.call(rbind, lapply(radii, function(r) {
        times <- unique(df$ms)
        data.frame("mean_accuracy" = sapply(times, function(t) {
                      mean(df[t - r <= df$ms & t + r >= df$ms, 'mean_accuracy'])
                   }),
                   "ms" = times,
                   "radius" = rep(r, length(times))
        )
    }))
    xbreaks <- seq(0, 19800, TIME_STEP)
    xlabels <- sapply(xbreaks,
                      function(t) {if (t %% 1000 == 0) as.character(t) else ""})
    ybreaks = c(seq(0, 1, .02), mean(df$mean_accuracy), 1/3)
    ylabels = c(seq(0, 1, .02), "media", "azar")
    ggplot(df2, aes(x = ms, y = mean_accuracy)) +
        geom_line(aes(group = radius, color = radius, size = as.factor(radius),
                      alpha = as.factor(radius))) +
        scale_x_continuous(breaks = xbreaks, minor_breaks = NULL,
                           labels = xlabels) +
        scale_y_continuous(breaks = ybreaks, minor_breaks = NULL,
                           labels = ylabels) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "Latencia estímulo-respuesta (ms)",
             y = "Exactitud de clasificación (media móvil)") +
        scale_color_gradient(trans = "log2") +
        scale_size_discrete(name = "Radio de suavización (ms)",
                            range = c(.5, 2.5)) +
        scale_alpha_discrete(name = "Radio de suavización (ms)",
                             range = c(1, 1/3)) +
        guides(color = FALSE)
}

user_maxima <- function(df) {
    best <- do.call(rbind, lapply(unique(df$subject), function(subject) {
        subdf <- df[df$subject == subject, ]
        best <- subdf[subdf$mean_accuracy == max(subdf$mean_accuracy), ]
        cbind(best[1, ], ocurrences = nrow(best))
    }))
    best <- best[order(best$mean_accuracy, decreasing = TRUE), ]
    best$subject <- factor(best$subject,
                                  levels = unique(best$subject))
    return(best)
}

plot_maxima_rank <- function(best) {
    breaks = seq(0, 1, .025)
    labels = sapply(breaks,
                    function(a) {if (a %% .05 == 0) as.character(a) else ""})
    ggplot(best, aes(x = subject, y = mean_accuracy)) +
        geom_line(aes(color = subject), show.legend = FALSE) +
        geom_point(#show.legend = FALSE,
                   shape = 21, stroke = 2.5,
                   aes(color = subject,
                       ## fill = factor(sample_size),
                       size = ocurrences,
                       alpha = best$ms)) +
        scale_alpha(name = "Latencia (ms)", range = c(1, .2)) +
        ## scale_fill_grey(name = "Muestras por clase", start = .8, end = .2) +
        scale_size(name = "Ocurrencias") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "Sujeto", y = "Exactitud de clasificación (maxima)") +
        scale_y_continuous(breaks = c(breaks, 1/3, mean(best$mean_accuracy)),
                           labels = c(labels, "azar", "media"),
                           minor_breaks = NULL) +
        guides(colour = FALSE)
}

plot_statistical_test <- function(nulls, best) {
    ggplot() +
        geom_flat_violin(aes(x = rep(-.25, nrow(nulls)), y = mean_accuracy,
                             group=subject, fill = subject), nulls,
                         adjust = .5, trim = FALSE, color = NA, alpha = .5) +
        geom_flat_violin(aes(x = rep(0, nrow(nulls)), y = mean_accuracy), nulls,
                         adjust = 1, trim = FALSE, color = NA) +
        geom_boxplot(aes( x = rep(0, nrow(nulls)), y = mean_accuracy),
                     nulls, outlier.shape = NA, alpha = .3, width = .01) +
        geom_flat_violin(aes(x = rep(0, nrow(best)), y = mean_accuracy), best,
                         fill = 1, adjust = .2, trim = FALSE, color = NA,
                         alpha = .3) +
        geom_boxplot(aes(x = rep(0, nrow(best)), y = mean_accuracy),
                     best, outlier.shape = NA, alpha = .3, width = .01) +
        geom_point(aes(x = .0129 * as.numeric(subject) - .48,
                       y = mean_accuracy, color = subject), best, size = 2) +
        scale_x_continuous(breaks = NULL) +
        scale_fill_discrete(name = "Sujeto") +
        scale_color_discrete(name = "Sujeto") +
        labs(x = "", y = "Exactitud de clasificación (máxima)") +
        scale_y_continuous(breaks = c(seq(0, 1, .1), 1/3,
                                      mean(best$mean_accuracy)),
                           labels = c(seq(0, 1, .1), "azar", "media"),
                           minor_breaks = NULL) +
        # guides(colour = FALSE, fill = FALSE) +
        coord_flip()
}

p_values <- function(nulls, best, bonferroni) {
    sapply(unique(best$subject), function(s) {
        h0 <- nulls[nulls$subject == s, 'mean_accuracy']
        p <- length(h0[h0 >= best[best$subject == s, 'mean_accuracy']]) /
            length(h0)
        p <- min(1 / bonferroni, p)
        names(p) <- s
        return(p)
    }) * bonferroni
}

corrected_best <- function(nulls, best, bonferroni) {
    sapply(unique(best$subject), function(s) {
        h0 <- sort(nulls[nulls$subject == s, 'mean_accuracy'])
        p <- bonferroni * 
             (length(h0[h0 >= best[best$subject == s, 'mean_accuracy']]) /
              length(h0))
        p <- min(1, p)
        new_acc <- if (p == 1) 0 else h0[length(h0) * (1 - p)]
        names(new_acc) <- s
        return(new_acc)
    })
}

best_sampling <- function(df2,
                          sampling_periods,
                          time_limits,
                          p_instead_of_acc = FALSE) {
    best_sampling <- matrix(ncol = length(time_limits),
                            nrow = length(sampling_periods))
    best_sampling <- do.call(rbind, mclapply(sampling_periods, function(i) {
        sapply(time_limits, function(j) {
            df3 <- df2[df2$ms %% i == 0 & df2$ms <= j, ]
            if (p_instead_of_acc) {
                mean(p_values(nulls,
                              user_maxima(df3),
                              length(unique(df3$ms))))
            } else {
                mean((user_maxima(df3))$mean_accuracy)
            }
        })
    }, mc.cores = detectCores() - 1))
    colnames(best_sampling) <- time_limits
    rownames(best_sampling) <- 1000 / sampling_periods
    return(best_sampling)
}

plot_best_sampling <- function(best_sampling, ztitle) {
    plot_ly(type = "surface",
            x = as.numeric(colnames(best_sampling)),
            y = as.numeric(rownames(best_sampling)),
            z = best_sampling) %>%
        layout(scene = list(xaxis = list(title = "Límite de tiempo (ms)"),
                            yaxis = list(title = "Tasa de muestreo (Hz)"),
                            zaxis = list(title = ztitle)))
}

# dataset loading and preparation ##############################################

time_series_files <- list.files(path = INPATH,
                                pattern = "result-time-series.txt",
                                full.names = TRUE,
                                recursive = TRUE)
null_dist_files <- list.files(path = INPATH,
                              pattern = "null-dist.txt",
                              full.names = TRUE,
                              recursive = TRUE)

df <- do.call(rbind, lapply(time_series_files, function(file) {
    df <- read.csv(file, header = FALSE, sep = " ")
    subject <- as.factor(regmatches(file, regexpr("\\d{3}", file)))
    cbind(df, rep(subject, nrow(df)), 0:(nrow(df) - 1) * TIME_STEP)
}))
names(df) <- c("sample_size", "mean_accuracy", "voxel_prop", "subject", "ms")

nulls <- do.call(rbind, lapply(null_dist_files, function(file) {
    df <- read.csv(file, header = FALSE, sep = " ")
    subject <- as.factor(regmatches(file, regexpr("\\d{3}", file)))
    cbind(df, rep(subject, nrow(df)))
}))
names(nulls) <- c("mean_accuracy", "subject")

# subject ids <= 526 have missing events on eprime files. discard them
df2 <- df[as.numeric(as.character(df$subject)) > 526 &
          df$sample_size == SAMPLE_SIZE, ]
## df2 <- df
nulls <- nulls[nulls$subject %in% df2$subject, ]
df2 <- df2[df2$subject %in% nulls$subject, ]

best <- user_maxima(df2)
nulls$subject <- factor(nulls$subject, levels = levels(best$subject))

# plots ########################################################################

svg(paste0(OUTPATH, '/timeseries.svg'), width = 20, height = 7)
plot(plot_timeseries(df2))
dev.off()

order <- list()
order$best <- rev(levels(best$subject))
order$first_max <-
    as.character(best[order(best$ms, best$mean_accuracy, decreasing = TRUE),
                      'subject'])
df2_matrix <- acast(df2[, c("subject", "ms", "mean_accuracy")],
                    subject ~ ms,
                    value.var = "mean_accuracy")
df2_matrix <- as.data.frame(df2_matrix)
subject_cluster <- hclust(dist(df2_matrix, method = "euclidean"),
                          method = "ward.D")
order$cluster <- sort(as.character(best$subject))[rev(subject_cluster$order)]
for (i in c('best', 'first_max', 'cluster')) {
    svg(paste0(OUTPATH, '/timeseries2-', i, '.svg'), width = 20, height = 7)
    plot(plot_timeseries_2(df2, order[i][[1]]))
    dev.off()
}
svg(paste0(OUTPATH, '/timeseries2-cluster-dendo.svg'))
plot(subject_cluster, xlab = "Sujeto", ylab = "Distancia (Ward)")
dev.off()

svg(paste0(OUTPATH, '/timeseries-mean.svg'), width = 20, height = 7)
plot(plot_mean_timeseries_denoise(df2, c(1, 500, 2000, 20000)))
dev.off()

svg(paste0(OUTPATH, '/user_maxima.svg'), width = 12, height = 10)
plot(plot_maxima_rank(best))
dev.off()

svg(paste0(OUTPATH, '/test.svg'), width = 10, height = 10)
plot(plot_statistical_test(nulls, best))
dev.off()

sampling_periods <- seq(max(df2$ms) / TIME_STEP, 1) * TIME_STEP
time_limits <- seq(200, 19800, TIME_STEP)
acc <- best_sampling(df2, sampling_periods, time_limits,
                     p_instead_of_acc = FALSE)
pval <- best_sampling(df2, sampling_periods, time_limits,
                      p_instead_of_acc = TRUE)
plot_best_sampling(acc, "Exactitud de clasificación máxima media")
plot_best_sampling(pval, "Valor p medio (Bonferroni)")
plot_best_sampling(acc / pval, "Exactitud máxima media / valor p (Bonferroni)")

## select samples with best corrected mean p-value
df3 <- df2[df2$ms %% 1400 == 0 & df2$ms <= 18200, ]
best2 <- user_maxima(df3)

svg(paste0(OUTPATH, '/user_maxima2.svg'), width = 12, height = 10)
plot(plot_maxima_rank(best2))
dev.off()

svg(paste0(OUTPATH, '/test2.svg'), width = 10, height = 10)
plot(plot_statistical_test(nulls, best2))
dev.off()

best3 <- best2
best3$mean_accuracy <- corrected_best(nulls, best2, length(unique(df3$ms)))
best3 <- best3[rev(order(best3$mean_accuracy)), ]
best3$subject <- as.factor(best3$subject)
svg(paste0(OUTPATH, '/test3.svg'), width = 10, height = 10)
plot(plot_statistical_test(nulls, best3))
dev.off()

svg(paste0(OUTPATH, '/pvalues.svg'))
plot(rev(sort(p_values(nulls, best2, length(unique(df3$ms))))), xlab = "Rango", ylab = "Valor p")
abline(.05, 0)
dev.off()
