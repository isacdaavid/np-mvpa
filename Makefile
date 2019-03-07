# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

SHELL := /bin/bash
BUILD_DIR := out
SRC_DIR := src
DATA_DIR := data

.PHONY: build all
all: build
build: images

.PHONY: images
images: $(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
	@mkdir -p "$(DATA_DIR)/xnat/$@"
	targets=($$(awk '{print $$1 ".zip"}' "$<" | sort)) ; \
	for i in $${targets[@]}; do \
	    if [[ ! -f "$(DATA_DIR)/xnat/$@/$$i" ]]; then \
	        $(SRC_DIR)/xnat/$@/xnat-download.sh $$(grep "^$${i%.zip} " $<) \
	                                            "$(DATA_DIR)/xnat/$@" ; \
	    fi ; \
	done ;

$(BUILD_DIR)/xnat/subject_metadata/fmri_subject_ids.csv: $(BUILD_DIR)
	@mkdir -p "$(BUILD_DIR)/xnat/subject_metadata"
	Rscript -e 'source("$(SRC_DIR)/xnat/subject_metadata/extract.R")'

$(BUILD_DIR):
	@mkdir -p "$(BUILD_DIR)"

.PHONY: clean
clean:
	@rm -rf "$(BUILD_DIR)"

