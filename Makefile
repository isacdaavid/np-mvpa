# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

SHELL := /bin/bash
BUILD_DIR := out
SRC_DIR := src
DATA_DIR := data

IDS_FILE := $(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
# note the use of the lazy assignment operator (strict evaluation) to avoid
# memoization of IDS after $(IDS_FILE) is regenerated
IDS = $(shell cut -d ' ' -f 1 $(IDS_FILE) | sort)
DICOMS = $(shell find $(DATA_DIR)/xnat/images/ -type d -name DICOM -printf "%p\n" | \
                 grep -E '(fMRI_GazeCueing|T1|FSPGR|T2)' | sort)

.PHONY: build all
all: build
build: eprime nifti

.PHONY: nifti
nifti: $(DICOMS:DICOM=nifti.nii.gz)

%nifti.nii.gz: %DICOM
	@echo 'building $@'
	@dcm2niix -f 'nifti' -g y -i y -t y -z y "$</.." > /dev/null

.PHONY: eprime
eprime: $(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
	@rm -r "$(BUILD_DIR)/$@" # FIXME
# copy eprime event files into subject-specific directory structure
	@targets=($$(cut -d ' ' -f 1 "$<")) ; \
	for i in $${targets[@]}; do \
	    mkdir -p "$(BUILD_DIR)/$@/$${i}" ; \
	    find $(DATA_DIR)/$@/{victor,FEDERICA/TASK_Gaze\ Cueing} \
	         -type f -name "*$${i}*.txt" \
	         -exec cp {} "$(BUILD_DIR)/$@/$${i}" \; ; \
	done
# manually fix duplicates and repeats
	@rm $(BUILD_DIR)/$@/526/*RTs* \
	    $(BUILD_DIR)/$@/677/Gaze* \
	    "$(BUILD_DIR)/$@/682/Gaze Cueing_B Backup1 Backup1-682-1.txt" \
	    $(BUILD_DIR)/$@/678/{'Gaze Cueing_B Backup1-678-1.txt','Copia de Gaze Cueing_C Backup1-678-1.txt'} \
	    # $(BUILD_DIR)/$@/678/'Copia de Copia de Gaze Cueing_A Backup1 Backup1-678-1.txt'
	@mv $(BUILD_DIR)/$@/664/Copia\ de\ Gaze\ Cueing_{C,B}\ Backup1-664-1.txt
# homogenize disparate file names into {A,B,C}.txt
	@find $(BUILD_DIR)/$@ -type f -name '*.txt' -exec bash -c \
	    'mv "{}" $$(echo {} | sed -En "s/(.*)(\/)(.*)(_)(A|B|C)(.*)(.txt)/\1\2\5\7/p")' \;
# UTF-16 -> UTF-8, CRLF -> LF
	@find $(BUILD_DIR)/$@ -type f -name '*.txt' -exec bash -c \
	    'iconv -f UTF-16 -t UTF-8 "{}" | tr -d "\r" > "{}.new" && mv "{}.new" "{}"' \;
# eprime event list -> pyMVPA sample attribute matrix
	@find $(BUILD_DIR)/$@ -type f -name '*.txt' -exec bash -c \
	    'awk -f "$(SRC_DIR)/eprime/eprime-to-csv.awk" -- "{}" > "{}.csv"' \;

.PHONY: dicoms
dicoms: $(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
	@mkdir -p "$(DATA_DIR)/xnat/$@"
	@targets=($$(cut -d ' ' -f 1 "$<" | sort)) ; \
	for i in $${targets[@]}; do \
	    [[ -d "$(DATA_DIR)/xnat/$@/$$i" ]] && continue ; \
	    if [[ ! -f "$(DATA_DIR)/xnat/$@/$${i}.zip" ]]; then \
	        $(SRC_DIR)/xnat/$@/xnat-download.sh $$(grep "^$$i " $<) \
	                                            "$(DATA_DIR)/xnat/$@" ; \
	    fi ; \
	    unzip -q "$(DATA_DIR)/xnat/$@/$${i}.zip" \
	          -d "$(DATA_DIR)/xnat/$@/" && rm "$(DATA_DIR)/xnat/$@/$${i}.zip" ; \
	done ;

$(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv: $(SRC_DIR)/xnat/subject_metadata/extract.R
	@mkdir -p "$(BUILD_DIR)/xnat/subject_metadata"
	Rscript -e 'source("$(SRC_DIR)/xnat/subject_metadata/extract.R")'

.PHONY: clean
clean:
	@rm -rf "$(BUILD_DIR)"

