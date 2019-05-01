## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

library(ggplot2)
library(dplyr)
source("src/poststats/R_rainclouds.R")

INPATH <- 'out/pymvpa/'
OUTPATH <- 'out/poststats/'
SAMPLE_SIZE <- 16

plot_timeseries <- function(df) {
    ggplot(df, aes(x = ms,
                       y = mean_accuracy,
                       group = subject,
                       color = subject)) +
    geom_line(aes(alpha=.01), show.legend = FALSE) +
    scale_x_continuous(breaks = seq(0, 19800, 200)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
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
        scale_alpha(range = c(1, .2)) +
        scale_fill_grey(start = .8, end = .2) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        guides(colour = FALSE)
}

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
    cbind(df, rep(subject, nrow(df)), 0:(nrow(df) - 1) * 200)
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

ggplot() +
geom_line(aes(x=unique(df2$ms), y=sapply(unique(df2$ms), function(t) {mean(df2[df2$ms == t, 'mean_accuracy'])})), color=4) +
scale_x_continuous(breaks = seq(0, 19800, 200)) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

svg(paste0(OUTPATH, '/timeseries.svg'))
plot(plot_timeseries(df2))
dev.off()

best <- user_maxima(df2)
nulls$subject <- factor(nulls$subject, levels = levels(best$subject))

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
