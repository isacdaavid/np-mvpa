## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

library(ggplot2)

PATH <- 'out/feat/'
files <- list.files(path = PATH,
                    pattern = "result-time-series.txt",
                    full.names = TRUE,
                    recursive = TRUE)

df <- do.call(rbind, lapply(files, function(file) {
    df <- read.csv(file, header = FALSE, sep = " ")
    subject <- as.factor(regmatches(file, regexpr("\\d{3}", file)))
    run <- regmatches(file, regexpr("(1|2|3)(rep)?.feat", file))
    run <- as.factor(regmatches(run, regexpr("\\d", run)))
    cbind(df, rep(subject, nrow(df)), rep(run, nrow(df)), 1:nrow(df) * 100)
}))
names(df) <- c("sample_size", "mean_accuracy", "voxel_prop", "subject", "run", "ms")

svg("plt.svg")
ggplot(df, aes(x = ms,
               y = mean_accuracy,
               group = interaction(subject, run),
               color = subject)) + geom_line(aes(alpha=.01), show.legend = FALSE)
dev.off()

best_scores <- do.call(rbind, lapply(unique(df$subject), function(subject) {
    subdf <- df[df$subject == subject, ]
    subdf <- do.call(rbind, lapply(unique(subdf$run), function(run) {
        rundf <- subdf[subdf$run == run, ]
        best <- rundf[rundf$mean_accuracy == max(rundf$mean_accuracy), ]
        cbind(best[1, ], ocurrences = nrow(best))
    }))
    cbind(subdf, mean_mean_accuracy = rep(mean(subdf$mean_accuracy), nrow(subdf)))
}))
best_scores <-  best_scores[order(best_scores$mean_mean_accuracy,
                                  decreasing = TRUE), ]
best_scores$subject <- factor(best_scores$subject,
                              levels = unique(best_scores$subject))

colors = c()
col = TRUE
prev = best_scores[1, "subject"]
for (s in best_scores$subject) {
    if (s != prev) { col = !col }
    colors = c(colors, col)
    prev = s
}

svg("plt2.svg")
ggplot(best_scores, aes(x = subject, y = mean_accuracy)) +
    geom_point(aes(y = mean_mean_accuracy)) +
    geom_line(aes(color = subject), show.legend = FALSE) +
    geom_jitter(aes(color = subject,
                    size = ocurrences,
                    alpha = mean(ms) - abs(mean(ms) - ms)), show.legend = FALSE,
                shape = 21, stroke = 1.5, height = 0, width = 0) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    guides(colour = FALSE)
dev.off()
