## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

library(oro.nifti)
library(neurobase)

cosa <- readnii(paste0(PATH, "/mean-weights-T1.nii.gz"))
cosa2 <- cosa[cosa != 0]
threshold <- sort(cosa2)[length(cosa2) * .99]
cosa[cosa < threshold] <- 0
writenii(cosa, paste0(PATH, '/mean-weights-T1-99th.nii.gz'))
