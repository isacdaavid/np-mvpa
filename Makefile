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
FMRI_NIFTIS = $(shell find $(DATA_DIR)/xnat/images/ -type f \
                      -name '*nifti.nii.gz' | grep 'tr_FMRI' | sort)
FEAT_NIFTIS = $(subst /resources/nifti.nii.gz,.feat, \
                      $(subst xnat/images,feat,$(FMRI_NIFTIS)))
FEAT_CONCAT = $(shell find $(DATA_DIR)/pymvpa -type d -name 'mc.feat')

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
	mat=$$(find "$(DATA_DIR)/pymvpa/$$id" -name '*func2highres.mat' | head -n1) ; \
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

$(BUILD_DIR)/pymvpa/% : $(DATA_DIR)/pymvpa/%/mc.feat/filtered_func_data_norm_brain.nii.gz
	@echo "running pyMVPA for $<"
	@mkdir -p "$@" ; \
	mask=$$(find "$(subst $(BUILD_DIR),$(DATA_DIR),$@)" -name 'volbrain-mask.nii.gz' | head -n 1) ; \
	id=$@ ; id=$${id##*/} ; \
	python2 "$(SRC_DIR)/pymvpa/pymvpa.py" "$(DATA_DIR)/psychopy/$${id}.csv" \
	        "$<" "$$mask" "$@" "$<" # > /dev/null 2>&1

#.PHONY : detrend_normalize
#detrend_normalize : $(FEAT_CONCAT:%=%/filtered_func_data_brain_norm.nii.gz)

#$(DATA_DIR)/pymvpa/%/mc.feat/filtered_func_data_brain_norm.nii.gz : $(DATA_DIR)/pymvpa/%/mc.feat/filtered_func_data_brain.nii.gz
#	@echo 'Detrending and normalizing into $@'
#	@python2 "$(SRC_DIR)/pymvpa/detrend_normalize.py" "$@" "$<" > /dev/null 2>&1

################################################################################
# post-FEAT gray matter extraction-related rules
################################################################################

.PHONY : feat_brains
feat_brains : $(FEAT_CONCAT:%=%/filtered_func_data_norm_brain.nii.gz)
	@echo

%/filtered_func_data_norm_brain.nii.gz : %/volbrain-mask.nii.gz %/filtered_func_data_norm.nii.gz
	@echo 'extracting BOLD 4D brain to $@'
	@fslmaths "$<" -uthr 2 "$<"
	@fslmaths "$<" -thr 2 "$<"
	@fslmaths "$<" -div 2 "$<"
	@fslmaths "$(subst _brain,,$@)" -mul "$<" "$@"

%/filtered_func_data.nii.gz : % ;

.PHONY : feat_masks
feat_masks : $(FEAT_CONCAT:%=%/volbrain-mask.nii.gz)
	@echo

#TODO: add prerequisite to _brain.nii.gz images, so as to complete dependency graph.
#      deduplicate per-scan masks, since they are all the same for a given subj.
%/volbrain-mask.nii.gz : %/reg/highres2example_func.mat %/example_func.nii.gz
	@echo 'creating BOLD mask $@'
	@orig_mask_dir=$(subst mc.feat/volbrain-mask.nii.gz,,$(subst pymvpa,volbrain,$@)) ; \
	orig_mask=$$(find "$$orig_mask_dir" -name '*crisp*') ; \
	flirt -interp nearestneighbour \
	      -in "$$orig_mask" \
	      -ref "$(subst volbrain-mask,example_func,$@)" \
	      -applyxfm \
	      -init "$(subst volbrain-mask.nii.gz,reg/highres2example_func.mat,$@)" \
	      -out "$@"

%/reg/highres2example_func.mat : % ;

%/example_func.nii.gz : % ;

################################################################################
# FSL FEAT preprocessing
################################################################################

.PHONY : apply_motion_correction
apply_motion_correction : $(FEAT_CONCAT:%=%/filtered_func_data_norm.nii.gz)

$(DATA_DIR)/pymvpa/%/mc.feat/filtered_func_data_norm.nii.gz : $(DATA_DIR)/pymvpa/%/concat-norm.nii.gz
	@echo "applying motion correction into $@"
	@matdir=$@ ; matdir=$${matdir%/*}/mc/prefiltered_func_data_mcf.mat ; \
	applyxfm4D "$<" "$<" "$@" "$$matdir" -fourdigit

.PHONY : concatenate_detrend_norm
concatenate_detrend_norm : $(addsuffix /concat-norm.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

$(DATA_DIR)/pymvpa/%/concat-norm.nii.gz : $(DATA_DIR)/feat/%
	@echo 'detrending, normalizing and concatenating runs into $@'
	@dir=$@ ; mkdir -p "$${dir%/*}" ; \
	fmris=($$(find "$<" -type f -name 'filtered_func_data.nii.gz' | \
	          sort --version-sort)) ; \
	python2 "$(SRC_DIR)/pymvpa/concatenate-detrend-norm.py" "$@" \
	        $${fmris[@]} > /dev/null 2>&1

.PHONY : concat_motion_correction
concat_motion_correction : $(addsuffix //mc.feat/filtered_func_data.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))

$(DATA_DIR)/pymvpa/%/mc.feat/filtered_func_data.nii.gz : $(DATA_DIR)/pymvpa/%
	@echo "MCFlirt'ing into $@"
	@t1dir=$(subst pymvpa,volbrain,$<) ; \
	t1=$$(find "$$t1dir" -name '*_brain.nii.gz') ; \
	sed "s|MVPA_OUTPUTDIR|$$(pwd)/$</mc| ; s|MVPA_FEAT_FILES|$$(pwd)/$</concat.nii.gz| ; s|MVPA_HIGHRES_FILES|$$(pwd)/$$t1|" \
	    "$(SRC_DIR)/feat/concat-motion-correction.fsf" > "$</concat-motion-correction.fsf" ; \
	feat "$</concat-motion-correction.fsf" ; sleep 5m

.PHONY : concatenate
concatenate : $(addsuffix /concat.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

$(DATA_DIR)/pymvpa/%/concat.nii.gz : $(DATA_DIR)/feat/%
	@echo 'concatenating runs into $@'
	@dir=$@ ; mkdir -p "$${dir%/*}" ; \
	fmris=($$(find "$<" -type f -name 'filtered_func_data.nii.gz' | \
	          sort --version-sort)) ; \
	python2 "$(SRC_DIR)/pymvpa/concatenate.py" "$@" \
	        $${fmris[@]} > /dev/null 2>&1

.PHONY : feat_prepro
#TODO: add prerequisite to _brain.nii.gz images, so as to complete dependency graph
feat_prepro : $(addsuffix .feat,$(subst /resources/nifti.nii.gz,,$(subst xnat/images,feat,$(FMRI_NIFTIS))))
	@echo

$(DATA_DIR)/feat/%.feat : $(DATA_DIR)/xnat/images/%/resources/nifti.nii.gz $(SRC_DIR)/feat/design.fsf
	@echo 'FEAT preprocessing at $@'
	@t1dir=$(subst feat,volbrain,$@) ; \
	t1dir=$${t1dir%scans*} ; \
	t1=$$(find "$$t1dir" -name '*_brain.nii.gz') ; \
	featdir=$@ ; mkdir -p "$${featdir%/*}" ; \
	sed "s|MVPA_OUTPUTDIR|$$(pwd)/$@| ; s|MVPA_FEAT_FILES|$$(pwd)/$<| ; s|MVPA_HIGHRES_FILES|$$(pwd)/$$t1|" \
	    "$(SRC_DIR)/feat/design.fsf" > "$@.design.fsf" ; \
	feat "$@.design.fsf" ; sleep 5m

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

