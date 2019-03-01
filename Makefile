# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

SHELL := /bin/bash
BUILD_DIR := out
SRC_DIR := src

.PHONY: build all
all: build
build: subject_metadata

.PHONY: subject_metadata
subject_metadata: $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/xnat/$@
	Rscript -e 'source("$(SRC_DIR)/xnat/$@/extract.R")'

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

.PHONY: clean
clean:
	@\rm -rf $(BUILD_DIR)
