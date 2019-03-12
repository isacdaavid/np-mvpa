# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

SHELL := /bin/bash
BUILD_DIR := out
SRC_DIR := src
DATA_DIR := data

.PHONY: build all
all: build
build: images eprime

.PHONY: eprime
eprime: $(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
# copy eprime event files into subject-specific directory structure
	targets=($$(awk '{print $$1}' "$<")) ; \
	for i in $${targets[@]}; do \
	    mkdir -p "$(BUILD_DIR)/$@/$${i}" ; \
	    find $(DATA_DIR)/$@/{victor,FEDERICA/TASK_Gaze\ Cueing} \
	         -type f -name "*$${i}*.txt" \
	         -exec cp {} "$(BUILD_DIR)/$@/$${i}" \; ; \
	done
# manually fix duplicates and repeats
	rm $(BUILD_DIR)/$@/526/*RTs* \
	    $(BUILD_DIR)/$@/677/Gaze* \
	    $(BUILD_DIR)/$@/678/{'Copia de Copia de Gaze Cueing_A Backup1 Backup1-678-1.txt','Gaze Cueing_B Backup1-678-1.txt','Copia de Gaze Cueing_C Backup1-678-1.txt'} \
	    "$(BUILD_DIR)/$@/682/Gaze Cueing_B Backup1 Backup1-682-1.txt"
	mv $(BUILD_DIR)/$@/664/Copia\ de\ Gaze\ Cueing_{C,B}\ Backup1-664-1.txt
# homogenize disparate file names into {A,B,C}.txt
	find $(BUILD_DIR)/$@ -type f -exec bash -c \
	    'mv "{}" $$(echo {} | sed -En "s/(.*)(\/)(.*)(_)(A|B|C)(.*)(.txt)/\1\2\5\7/p")' \;

.PHONY: images
images: $(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
	@mkdir -p "$(DATA_DIR)/xnat/$@"
	@targets=($$(awk '{print $$1 ".zip"}' "$<" | sort)) ; \
	for i in $${targets[@]}; do \
	    if [[ ! -f "$(DATA_DIR)/xnat/$@/$$i" ]]; then \
	        $(SRC_DIR)/xnat/$@/xnat-download.sh $$(grep "^$${i%.zip} " $<) \
	                                            "$(DATA_DIR)/xnat/$@" ; \
	    fi ; \
	done ;

$(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv: $(SRC_DIR)/xnat/subject_metadata/extract.R
	@mkdir -p "$(BUILD_DIR)/xnat/subject_metadata"
	Rscript -e 'source("$(SRC_DIR)/xnat/subject_metadata/extract.R")'

.PHONY: clean
clean:
	@rm -rf "$(BUILD_DIR)"

