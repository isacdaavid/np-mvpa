## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

library(ggplot2)
library(dplyr)
source("src/poststats/R_rainclouds.R")

INPATH <- 'out/pymvpa/'
OUTPATH <- 'out/poststats/'
SAMPLE_SIZE <- 16
TIME_STEP <- 200 # ms

plot_timeseries <- function(df) {
    best <- user_maxima(df)
    xbreaks <- seq(0, 19800, TIME_STEP)
    xlabels <- sapply(xbreaks,
                      function(t) {if (t %% 1000 == 0) as.character(t) else ""})
    ybreaks = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1,
                mean(df$mean_accuracy), 1/3)
    ylabels = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1,
                "media", "azar")
    ggplot(df, aes(x = ms,
                   y = mean_accuracy,
                   group = subject,
                   color = subject)) +
        geom_line(aes(alpha=.01), show.legend = FALSE) +
	geom_point(aes(x = ms, y = mean_accuracy, color = subject), best) +
        scale_x_continuous(breaks = xbreaks, minor_breaks = NULL, labels = xlabels) +
	scale_y_continuous(breaks = ybreaks, minor_breaks = NULL, labels = ylabels) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
	labs(x="Retraso estímulo-respuesta (ms)", y="Exactitud de clasificación") +
	guides(colour = FALSE)
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
        geom_line(aes(group = -radius, color = radius, size = as.factor(radius))) +
        scale_x_continuous(breaks = xbreaks, minor_breaks = NULL, labels = xlabels) +
	scale_y_continuous(breaks = ybreaks, minor_breaks = NULL, labels = ylabels) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
	labs(x="Retraso estímulo-respuesta (ms)", y="Exactitud media de clasificación") +
        scale_color_gradient(name = "Radio de suavización (ms)", trans = "log2") +
        scale_size_manual(values = c(.5, 1, 1.5, 2, 2.5)) +
        guides(size = FALSE)
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
    colors = c()
    col = TRUE
    prev = best[1, "subject"]
    for (s in best$subject) {
        if (s != prev) { col = !col }
        colors = c(colors, col)
        prev = s
    }

    ggplot(best, aes(x = subject, y = mean_accuracy)) +
        geom_line(aes(color = subject), show.legend = FALSE) +
        geom_point(#show.legend = FALSE,
                   shape = 21, stroke = 1.5,
                   aes(color = subject,
                       fill = factor(sample_size),
                       size = ocurrences,
                       alpha = abs(mean(best$ms) - best$ms))) +
        scale_alpha(name = "Separación a tiempo medio (ms)", range = c(1, .2)) +
        scale_fill_grey(name = "Muestras por clase", start = .8, end = .2) +
	scale_size(name = "Ocurrencias") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
	labs(x = "Sujeto", y = "Exactitud máxima de clasificación") +
        guides(colour = FALSE)
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
          df$sample_size >= SAMPLE_SIZE, ]
# df2 <- df
nulls <- nulls[nulls$subject %in% df2$subject, ]
df2 <- df2[df2$subject %in% nulls$subject, ]

best <- user_maxima(df2)
nulls$subject <- factor(nulls$subject, levels = levels(best$subject))

p_values <- sapply(unique(best$subject), function(s) {
    h0 <- nulls[nulls$subject == s, 'mean_accuracy']
    p <- length(h0[h0 >= best[best$subject == s, 'mean_accuracy']]) / length(h0)
    names(p) <- s
    return(p)
})

# plots ########################################################################

svg(paste0(OUTPATH, '/timeseries.svg'), width = 20, height = 3)
plot(plot_timeseries(df2))
dev.off()

svg(paste0(OUTPATH, '/timeseries-mean.svg'), width = 20, height = 3)
plot(plot_mean_timeseries_denoise(df2, c(1, 500, 1000, 2000, 20000)))
dev.off()

svg(paste0(OUTPATH, '/user_maxima.svg'))
plot(plot_maxima_rank(best))
dev.off()

svg(paste0(OUTPATH, '/test.svg'))
# revert order because fill=subject is an idiot who also does it
levels(best$subject) <- sort(levels(best$subject))
ggplot() +
    geom_abline(slope = 0, intercept = 1/3, color='white') +
    geom_flat_violin(aes(x=rep(-.25, nrow(nulls)), y=mean_accuracy, fill=subject),
                     nulls, adjust = .5, trim = FALSE, color = NA, alpha = .5) +
    geom_flat_violin(aes(x=rep(0, nrow(nulls)), y=mean_accuracy),
                     nulls, adjust = 1, trim = FALSE, color = NA, alpha = 1) +
    geom_boxplot(aes(x=rep(0, nrow(nulls)), y=mean_accuracy),
                 nulls, outlier.shape=NA, alpha = .3, width = .01) +
    geom_flat_violin(aes(x=rep(0, nrow(best)), y=mean_accuracy, fill=as.factor(680)),
                     best, adjust = .2, trim = FALSE, color = NA, alpha = .5) +
    geom_boxplot(aes(x=rep(0, nrow(best)), y=mean_accuracy),
                 best, outlier.shape=NA, alpha = .3, width = .01) +
    geom_jitter(aes(x=rep(-.25, nrow(best)),
                    y=mean_accuracy, color=subject),
                best, size = 2, alpha = .5, height = 0, width = .2) +
    guides(colour = FALSE, fill = FALSE) +
    scale_x_continuous(breaks = NULL) +
    labs(x="", y="Classification accuracy") +
    coord_flip()
dev.off()
