# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

SHELL := /bin/bash

BUILD_DIR := out
SRC_DIR := src
DATA_DIR := data

.PHONY : build all
all : build
build : concatenate_runs

IDS_FILE := $(DATA_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
# note the use of the lazy assignment operator (strict evaluation) to avoid
# memoization of IDS after $(IDS_FILE) is regenerated
IDS = $(shell cut -d ' ' -f 1 $(IDS_FILE) | sort)
DICOMS = $(shell find $(DATA_DIR)/xnat/images/ -type d -name DICOM | \
                 grep -E '(FMRI|RestState|T1|T2)' | sort)
VOLBRAIN_ZIPS = $(shell find $(DATA_DIR)/volbrain/ -type f -name '*.zip')
VOLBRAIN_IMAGES = $(shell find data/volbrain/ -type f -name 'n_mmni*')

################################################################################
# transform back to T1w space
################################################################################

.PHONY : register_results
register_results : $(addsuffix /all-weights-T1.nii.gz,  $(addprefix $(BUILD_DIR)/pymvpa/, $(IDS)))

%-weights-T1.nii.gz : %-weights.nii.gz
	@echo 'transforming $< to T1w space'
	@id=$(subst $(BUILD_DIR)/pymvpa/,,$<) ; \
	id=$${id%/*} ; \
	t1=$$(find "$(DATA_DIR)/volbrain/$$id" -name '*_brain.nii.gz') ; \
	mat=$$(find "$(DATA_DIR)/feat/$$id" -name '*func2highres.mat' | head -n1) ; \
	flirt -interp trilinear \
	      -in "$<" \
	      -ref "$$t1" \
	      -applyxfm \
	      -init "$$mat" \
	      -out "$@"

################################################################################
# pyMVPA rules
################################################################################

.PHONY : pymvpa
pymvpa : $(addprefix $(BUILD_DIR)/pymvpa/, $(IDS))

$(BUILD_DIR)/pymvpa/% : $(DATA_DIR)/pymvpa/%/concat-brain-norm.nii.gz
	@echo "running pyMVPA for $<"
	@mkdir -p "$@" ; \
	mask=$$(find "$(subst $(BUILD_DIR)/pymvpa,$(DATA_DIR)/feat,$@)" -name 'volbrain-mask.tmp.nii.gz') ; \
	id=$@ ; id=$${id##*/} ; \
	python2 "$(SRC_DIR)/pymvpa/pymvpa.py" "$(DATA_DIR)/psychopy/$${id}.csv" \
	        "$<" "$$mask" "$@" "$<" # > /dev/null 2>&1

.PHONY : detrend_normalize
detrend_normalize : $(addsuffix /concat-brain-norm.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

%/concat-brain-norm.nii.gz : %/concat-brain.nii.gz
	@echo 'detrending and normalizing into $@'
	@python2 "$(SRC_DIR)/pymvpa/detrend-normalize.py" "$@" "$<"  > /dev/null 2>&1

################################################################################
# post-FEAT gray matter extraction-related rules
################################################################################

.PHONY : feat_brains
feat_brains : $(addsuffix /concat-brain.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

# TODO: missing explicit dependency on mask
$(DATA_DIR)/pymvpa/%/concat-brain.nii.gz : $(DATA_DIR)/feat/%/feat.feat/filtered_func_data.nii.gz
	@echo 'extracting BOLD 4D brain to $@'
	@mask=$(subst concat-brain.nii.gz,volbrain-mask,$(subst pymvpa,feat,$@)) ; \
	fslmaths "$${mask}.nii.gz" -uthr 2 "$${mask}.tmp.nii.gz" ; \
	fslmaths "$${mask}.tmp.nii.gz" -thr 2 "$${mask}.tmp.nii.gz" ; \
	fslmaths "$${mask}.tmp.nii.gz" -div 2 "$${mask}.tmp.nii.gz" ; \
	fslmaths "$<" -mul "$${mask}.tmp.nii.gz" "$@"

%/filtered_func_data.nii.gz : % ;

.PHONY : feat_masks
feat_masks : $(addsuffix /volbrain-mask.nii.gz, $(addprefix $(DATA_DIR)/feat/, $(IDS)))
	@echo

#TODO: add prerequisite to _brain.nii.gz images, so as to complete dependency graph.
%/volbrain-mask.nii.gz : %/feat.feat/reg/highres2example_func.mat %/feat.feat/example_func.nii.gz
	@echo 'creating BOLD mask $@'
	@orig_mask_dir=$(subst volbrain-mask.nii.gz,,$(subst feat,volbrain,$@)) ; \
	orig_mask=$$(find "$$orig_mask_dir" -name '*crisp*') ; \
	flirt -interp nearestneighbour \
	      -in "$$orig_mask" \
	      -ref "$(subst volbrain-mask,feat.feat/example_func,$@)" \
	      -applyxfm \
	      -init "$<" \
	      -out "$@"

%/reg/highres2example_func.mat : % ;

%/example_func.nii.gz : % ;

################################################################################
# FSL FEAT preprocessing
################################################################################

.PHONY : feat_prepro
feat_prepro : $(addsuffix //feat.feat, $(addprefix $(DATA_DIR)/feat/, $(IDS)))
	@echo

# FIXME: fix dependency to be on pymvpa/concat.nii.gz, not just pymvpa/
$(DATA_DIR)/feat/%/feat.feat : $(DATA_DIR)/pymvpa/%
	@echo "FEAT preprocessing into $@"
	@featdir=$(subst pymvpa,feat,$<) ; mkdir -p "$$featdir" ; \
	t1dir=$(subst pymvpa,volbrain,$<) ; \
	t1=$$(find "$$t1dir" -name '*_brain.nii.gz') ; \
	nvols=$$(fslnvols "$</concat.nii.gz") ; \
	sed "s|MVPA_OUTPUTDIR|$$(pwd)/$${featdir}/feat| ; s|MVPA_FEAT_FILES|$$(pwd)/$</concat.nii.gz| ; s|MVPA_HIGHRES_FILES|$$(pwd)/$${t1}| ; s|MVPA_NVOLS|$$nvols|" \
	    "$(SRC_DIR)/feat/design.fsf" > "$${featdir}/design.fsf" ; \
	feat "$${featdir}/design.fsf"

.PHONY : concatenate
concatenate : $(addsuffix /concat.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

$(DATA_DIR)/pymvpa/%/concat.nii.gz : $(DATA_DIR)/xnat/images/%
	@echo 'concatenating runs into $@'
	@dir=$@ ; mkdir -p "$${dir%/*}" ; \
	fmris=($$(find "$<" -type f -name '*nifti.nii.gz' | grep 'tr_FMRI' | sort --version-sort)) ; \
	python2 "$(SRC_DIR)/pymvpa/concatenate.py" "$@" $${fmris[@]} > /dev/null 2>&1

################################################################################
# volbrain-related rules
################################################################################

.PHONY : t1w_brain_extraction
t1w_brain_extraction : $(VOLBRAIN_IMAGES:.nii=_brain.nii.gz)
	@echo

$(DATA_DIR)/volbrain/%_brain.nii.gz : $(DATA_DIR)/volbrain/%.nii
	@echo 'extracting brain to $@'
	@fslmaths $< -mul "$(subst n_,mask_n_,$<)" "$@"

.PHONY : volbrain_unzip
volbrain_unzip : $(VOLBRAIN_ZIPS:.zip=.volbrain)
	@echo

# FIXME: requires manual upload/download from https://volbrain.upv.es
%.volbrain : %.zip
	@echo 'unzip-ing $<'
	@mkdir "$@" && unzip -q "$<" -d "$@" && rm "$<"

# create directory tree where to put volbrain's results
.PHONY : volbrain_tree
volbrain_tree : nifti $(addsuffix /, $(addprefix $(DATA_DIR)/volbrain/, $(IDS)))
	@echo

$(DATA_DIR)/volbrain/%/ :
	@echo 'creating directory at $@'
	@mkdir -p "$@"

################################################################################
# convert xnat DICOMs to Nifti
################################################################################

.PHONY : nifti
nifti : images $(DICOMS:DICOM=nifti.nii.gz)
	@echo

%nifti.nii.gz : %DICOM
	@echo 'building $@'
	@dcm2niix -f 'nifti' -g y -i y -t y -z y "$</.." > /dev/null

################################################################################
# xnat DICOMs: download, unzip. DO NOT PARALLELIZE (don't run with -j )
################################################################################

.PHONY : images
images : $(IDS_FILE)
	@mkdir -p "$(DATA_DIR)/xnat/$@"
	@targets=($$(cut -d ' ' -f 1 "$<" | sort)) ; \
	for i in $${targets[@]}; do \
	    [[ -d "$(DATA_DIR)/xnat/$@/$$i" ]] && continue ; \
	    if [[ ! -f "$(DATA_DIR)/xnat/$@/$${i}.zip" ]]; then \
	        printf 'downloading DICOMs for subject %d\n\n' "$$i" ; \
	        $(SRC_DIR)/xnat/$@/xnat-download.sh $$(grep "^$$i " $<) \
	                                            "$(DATA_DIR)/xnat/$@" || \
	            rm "$(DATA_DIR)/xnat/$@/$${i}.zip" ; \
	    fi ; \
	    printf "unzip-ing DICOMs for subject %d\n\n" "$$i" ; \
	    unzip -q "$(DATA_DIR)/xnat/$@/$${i}.zip" \
	          -d "$(DATA_DIR)/xnat/$@/" && rm "$(DATA_DIR)/xnat/$@/$${i}.zip" ; \
	done ;

