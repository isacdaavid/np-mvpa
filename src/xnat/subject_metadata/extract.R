## author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
## license: GPLv3 or later

##  extract relevant XNAT subject metadata from all XML files. pour into single
##  CSV with the following format:
##
##     PATIENT_ID EXPERIMENT_ID NAME SATISFACTION #fMRI #T1 #T1_PU #DICOM #NIFTI
##
##  afterwards, extract PATIENT_IDs and EXPERIMENT_IDs with nonempty fMRI
##  sessions and eprime data for future use with ../images/xnat-download.sh

if (! "xml2" %in% rownames(installed.packages())) {
    install.packages("xml2")
}
library(xml2)

DATA_DIR <- './data/xnat/subject_metadata/experiments'
OUT <- './out/xnat/subject_metadata/subjects.csv'
OUT2 <- './out/xnat/subject_metadata/fmri_subject_ids.csv'

## per-file/subject extraction. subject('path-to-file')
subject <- function(xml) {
    file <- read_xml(xml)
    s <- list()
    ## s$pid <- xml_text(xml_find_all(file, '//xnat:dcmPatientId'))
    s$pid <- as.character(xml_find_all(file, '//xnat:experiment'))
    s$pid <- regmatches(s$pid, regexpr('label="[0-9]+"', s$pid))
    s$pid <- regmatches(s$pid, regexpr('[0-9]+', s$pid))
    s$expid <-  as.character(xml_find_all(file, '//xnat:experiment'))
    s$expid <- regmatches(s$expid,
                                  regexpr('ID="XNAT_E[0-9]+"', s$expid))
    s$expid <- regmatches(s$expid, regexpr('XNAT_E[0-9]+', s$expid))
    s$name <- xml_text(xml_find_all(file, '//xnat:dcmPatientName'))
    s$satisfaction <- xml_text(xml_find_all(file,
                                            '//xnat:field[@name="totalsatisfacción"]'))
    s$fmri <- length(xml_find_all(file,
                                  '//xnat:scan[@type="fMRI GazeCueing_1"]')) +
        length(xml_find_all(file,
                            '//xnat:scan[@type="fMRI GazeCueing_2"]')) +
        length(xml_find_all(file,
                            '//xnat:scan[@type="fMRI GazeCueing_3"]'))
    s$t1 <- length(xml_find_all(file,
                                '//xnat:scan[@type="Sag_T1_FSPGR_BRAVO"]'))
    s$t1_pu <- length(xml_find_all(file,
                                   '//xnat:scan[@type="PU:Sag_T1_FSPGR_BRAVO"]'))
    s$dicom <- length(xml_find_all(file,
                                   '//xnat:file[@label="DICOM"]'))
    s$nifti <- length(xml_find_all(file,
                                   '//xnat:file[@label="NIFTI"]'))
    return(s)
}

files <- list.files(DATA_DIR, full.names = TRUE)
subjects <- lapply(files, subject)

## output CSV
header <- c("pid", "expid", "name","satisfaction","fmri","t1","t1_PU","dicom","nifti")
write.table(as.matrix(t(header)),
            OUT,
            sep = ",",
            col.names = FALSE,
            row.names = FALSE,
            append = TRUE)
for (i in 1:length(subjects)) {
    ## remove NULL
    subjects[[i]][lapply(subjects[[i]], length) == 0] <- NA

    write.table(x = subjects[[i]],
                file = OUT,
                sep = ",",
                col.names = FALSE,
                row.names = FALSE,
                append = TRUE,
                fileEncoding = "UTF-8")
}

dataframe <- read.csv(OUT)

## subset of subject ids and experiment ids with full fMRI sessions
write.table(dataframe[dataframe[, "fmri"] >= 1 & dataframe[, "E.prime"] >= 1,
                      c("pid", "expid")],
            OUT2,
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
