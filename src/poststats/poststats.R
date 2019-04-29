## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

library(ggplot2)

PATH <- 'out/pymvpa/'

plot_timeseries <- function(df) {
    ggplot(df, aes(x = ms,
                       y = mean_accuracy,
                       group = subject,
                       color = subject)) +
    geom_line(aes(alpha=.01), show.legend = FALSE)
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

files <- list.files(path = PATH,
                    pattern = "result-time-series.txt",
                    full.names = TRUE,
                    recursive = TRUE)

df <- do.call(rbind, lapply(files, function(file) {
    df <- read.csv(file, header = FALSE, sep = " ")
    subject <- as.factor(regmatches(file, regexpr("\\d{3}", file)))
    cbind(df, rep(subject, nrow(df)), 1:nrow(df) * 100)
}))
names(df) <- c("sample_size", "mean_accuracy", "voxel_prop", "subject", "ms")

# subject ids <= 526 have missing events on eprime files. discard them
df2 <- df[as.numeric(as.character(df$subject)) > 526 &
          df$sample_size > 30, ]

svg('out/poststats/timeseries.svg')
plot(plot_timeseries(df2))
dev.off()

best <- user_maxima(df2)

svg('out/poststats/user_maxima.svg')
plot(plot_maxima_rank(best))
dev.off()

