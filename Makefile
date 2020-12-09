# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

SHELL := /bin/bash

BUILD_DIR := out
SRC_DIR := src
DATA_DIR := data

TASKNAME := emotionalfaces
RUNS := 5
PERMUTATIONS := 5000

IDS_FILE := $(DATA_DIR)/xnat/subject_metadata/fmri_subject_ids.csv
# note the use of the lazy assignment operator (strict evaluation) to avoid
# memoization of IDS after $(IDS_FILE) is regenerated
IDS = $(shell find $(DATA_DIR)/bids/ -maxdepth 1 -type d -name 'sub-*' | sort | cut -d '-' -f 2)
DICOMS = $(shell find $(DATA_DIR)/xnat/images/ -type d -name DICOM | \
                  grep -E '(FMRI|RestState|T1|T2)' | sort)
VOLBRAIN_ZIPS = $(shell find $(DATA_DIR)/volbrain/ -type f -name '*.zip')
VOLBRAIN_IMAGES = $(shell find data/volbrain/ -type f -name 'n_mmni*')
SENSITIVITY_MAPS = $(shell find $(BUILD_DIR)/pymvpa -type f -name '*.nii.gz' | grep -v T1)

################################################################################
# postpoststats
################################################################################

# summarize contrast results in a per-category CSV,
# compute mean p-val map per category
.PHONY : postpoststats
postpoststats :
	@for reduced in $$(ls $(BUILD_DIR)/pymvpa) ; do \
	    for category in $(SRC_DIR)/pymvpa/contrasts/* ; do \
	        cat=$$(basename "$$category") ; \
	        out="$(BUILD_DIR)/poststats/$$reduced/$$cat/stats.csv" ; \
	        printf "mean_p-val\tCohens_D\tcontrast\n" > "$$out" ; \
	        for f in $$(find $(BUILD_DIR)/poststats/$$reduced/$$cat/ -name 'stats.csv') ; do \
	            printf "%s\t%s\n" "$$(tail -n1 "$$f")" \
	                   "$$(echo $${f#*/*/*/*/} | sed 's|/stats.csv||')" ; \
	        done | sort >> "$$out" ; \
	        pvals_greater=($$(find "$(BUILD_DIR)/poststats/$$reduced/$$cat" -type f -name 'greater_zero_tfce_corrp_tstat1.nii.gz')) ; \
	        pvals_less=($$(find "$(BUILD_DIR)/poststats/$$reduced/$$cat" -type f -name 'less_zero_tfce_corrp_tstat1.nii.gz')) ; \
	        fslmerge -t "$(BUILD_DIR)/poststats/$$reduced/$$cat/greater_zero_tfce_corrp_tstat_MERGED.nii.gz" "$${pvals_greater[@]}" ; \
	        fslmaths "$(BUILD_DIR)/poststats/$$reduced/$$cat/greater_zero_tfce_corrp_tstat_MERGED.nii.gz" -Tmean "$(BUILD_DIR)/poststats/$$reduced/$$cat/greater_zero_tfce_corrp_tstat_MEAN.nii.gz" ; \
	        fslmerge -t "$(BUILD_DIR)/poststats/$$reduced/$$cat/less_zero_tfce_corrp_tstat_MERGED.nii.gz" "$${pvals_less[@]}" ; \
	        fslmaths "$(BUILD_DIR)/poststats/$$reduced/$$cat/less_zero_tfce_corrp_tstat_MERGED.nii.gz" -Tmean "$(BUILD_DIR)/poststats/$$reduced/$$cat/less_zero_tfce_corrp_tstat_MEAN.nii.gz" ; \
	    done ; \
	done ;

################################################################################
# poststats and group-level cluster inference
################################################################################

# compute FWE-corrected p-val map using cluster-informed randomise
.PHONY : group_level_glm
group_level_glm :
	@for contrast in $(SRC_DIR)/feat/3.level-3-lme-designs/* ; do \
	    outdir=$(BUILD_DIR)/feat/TFCE-randomise/$$(basename $${contrast%.fsf}) ; \
	    echo "running randomise for $$outdir" ; \
	    mkdir -p "$$outdir" ; \
	    paths=($$(grep -E cope.*.nii.gz "$$contrast" | cut -f 3 -d ' ' | tr -d '"' | sed 's|stats|reg_standard/stats|')) ; \
	    fslmerge -t "$${outdir}/merged_weights.nii.gz" "$${paths[@]}" ; \
	    fsl_sub /home/inb/lconcha/fmrilab_software/fsl_4.1.9/bin/randomise -i "$${outdir}/merged_weights.nii.gz" -o "$${outdir}/greater_zero" -1 -v 5 -T -n $(PERMUTATIONS) ; \
	    fslmaths "$${outdir}/merged_weights.nii.gz" -mul -1 "$${outdir}/merged_weights_neg.nii.gz" ; \
	    fsl_sub /home/inb/lconcha/fmrilab_software/fsl_4.1.9/bin/randomise -i "$${outdir}/merged_weights_neg.nii.gz" -o "$${outdir}/less_zero" -1 -v 5 -T -n $(PERMUTATIONS) ; \
	done

# compute FWE-corrected p-val map using cluster-informed randomise
.PHONY : group_level_mvpa
group_level_mvpa :
	@for reduced in "whole" ; do \
	    for category in $(SRC_DIR)/pymvpa/contrasts/* ; do \
	        while read contrast; do \
	            outdir=$(BUILD_DIR)/poststats/$${reduced}/$$(basename $$category)/$${contrast} ; \
	            echo "running randomise for $$outdir" ; \
	            mkdir -p "$$outdir" ; \
	            paths=($$(find "$(BUILD_DIR)/pymvpa/$$reduced" -type f -name "$${contrast}-weights-T1.nii.gz" | grep "/$${contrast}/" | sort)) ; \
	            fslmerge -t "$${outdir}/merged_weights.nii.gz" "$${paths[@]}" ; \
#	            fslmaths "$${outdir}/merged_weights.nii.gz" -abs "$${outdir}/merged_weights_abs.nii.gz" ; \
	            fslmaths "$${outdir}/merged_weights.nii.gz" -s 2.123 "$${outdir}/merged_weights_smooth.nii.gz" ; \
	            fsl_sub /home/inb/lconcha/fmrilab_software/fsl_4.1.9/bin/randomise -i "$${outdir}/merged_weights_smooth.nii.gz" -o "$${outdir}/greater_zero" -1 -v 5 -T -n $(PERMUTATIONS) ; \
	            fslmaths "$${outdir}/merged_weights_smooth.nii.gz" -mul -1 "$${outdir}/merged_weights_smooth_neg.nii.gz" ; \
	            fsl_sub /home/inb/lconcha/fmrilab_software/fsl_4.1.9/bin/randomise -i "$${outdir}/merged_weights_smooth_neg.nii.gz" -o "$${outdir}/less_zero" -1 -v 5 -T -n $(PERMUTATIONS) ; \
	        done < "$$category" ; \
	    done ; \
	done ;

# group-level hypothesis tests on classification accuracy, effect size, plots
.PHONY : poststats
poststats :
	@for reduced in $$(ls $(BUILD_DIR)/pymvpa) ; do \
	    for category in $(SRC_DIR)/pymvpa/contrasts/* ; do \
	        while read contrast; do \
	            outdir=$(BUILD_DIR)/poststats/$${reduced}/$$(basename $$category)/$${contrast} ; \
	            echo "running result statistics for $$outdir" ; \
	            mkdir -p "$$outdir" ; \
	            paths=$$(find "$(BUILD_DIR)/pymvpa/$$reduced" -type d -name "$${contrast}" -printf "'%p'\n" | sort | tr '\n' , ) ; \
	            nclasses=$$(awk -F , '{print NF}' <<< "$$contrast") ; \
	            Rscript -e "INPATH <- c($${paths::-1}) ; OUTPATH <- \"$$outdir\" ; NCLASSES <- $$nclasses ; source('$(SRC_DIR)/poststats/poststats.R')" ; \
	        done < "$$category" ; \
	    done ; \
	done ;

################################################################################
# univariate/GLM first-level analysis
################################################################################

.PHONY : feat_level1
feat_level1 :
	feat $(SRC_DIR)/feat/2.level-1-glm-design.fsf

################################################################################
# pyMVPA rules
################################################################################

# transform resulting sensitivity maps back to T1w space, then
# denoise to improve spatial detection at group analysis

.PHONY : register_results
register_results : $(SENSITIVITY_MAPS:.nii.gz=-T1.nii.gz)
	@echo

%-T1.nii.gz : %.nii.gz
	@echo 'transforming to T1w space $<'
	@id=$< ; id=$${id%/*/*/*} ; id=$${id##*/} ; \
	t1=$$(find "$(DATA_DIR)/volbrain/$$id" -name '*_brain.nii.gz') ; \
	mat=$$(find "$(DATA_DIR)/feat/$$id" -name '*func2standard.mat' | head -n1) ; \
	flirt -interp trilinear \
	      -in "$<" \
	      -ref "$$t1" \
	      -applyxfm \
	      -init "$$mat" \
	      -out "$(subst .nii.gz,-T1,$<)"

# train classifiers based on whole NIFTIs
.PHONY : pymvpa_whole
pymvpa_whole : $(addprefix $(BUILD_DIR)/pymvpa/whole/, $(IDS))
	@echo $<

$(BUILD_DIR)/pymvpa/whole/% : $(DATA_DIR)/pymvpa/%/concat-brain-norm.nii.gz
	@echo "running pyMVPA for $<"
	@mask=$$(find "$(subst $(BUILD_DIR)/pymvpa/whole,$(DATA_DIR)/feat,$@)" -name 'volbrain-mask.tmp.nii.gz') ; \
	id=$@ ; id=$${id##*/} ; \
	for category in $(SRC_DIR)/pymvpa/contrasts/* ; do \
	    while read contrast; do \
	        outdir=$@/$$(basename $${category})/$${contrast} ; \
	        mkdir -p "$$outdir" ; \
		tmpfile=s$${id}-$${contrast}.sh ; \
		echo python2 "$(SRC_DIR)/pymvpa/main.py" "$(DATA_DIR)/psychopy/$${id}.csv" "$<" "$$mask" "$$outdir" "$<" "$$contrast" $(PERMUTATIONS) >> $$tmpfile ; \
		chmod ugo+x $$tmpfile ; \
		qsub -l h_vmem=20G -l h='!(arwen.inb.unam.mx|giora.inb.unam.mx)' -V -cwd $$tmpfile ; \
	    done < "$$category" ; \
	done

# linearly detrend and z-normalize time series
.PHONY : detrend_normalize
detrend_normalize : $(addsuffix /concat-brain-norm.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

%/concat-brain-norm.nii.gz : %/concat-brain.nii.gz
	@echo 'detrending and normalizing into $@'
	@nvols=$$(fslnvols "$<") ; \
	python2 "$(SRC_DIR)/pymvpa/detrend-normalize.py" "$@" "$<" $$(($$nvols / $(RUNS))) > /dev/null 2>&1

################################################################################
# post-FEAT brain masking
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

# %/filtered_func_data.nii.gz : % ;

.PHONY : feat_masks
feat_masks : $(addsuffix /volbrain-mask.nii.gz, $(addprefix $(DATA_DIR)/feat/, $(IDS))) $(addsuffix /atlas.nii.gz, $(addprefix $(DATA_DIR)/feat/, $(IDS)))
	@echo

%/volbrain-mask.nii.gz : %/feat.feat/reg/highres2example_func.mat %/feat.feat/example_func.nii.gz
	@echo 'creating brain mask $@'
	@orig_mask_dir=$(subst volbrain-mask.nii.gz,,$(subst feat,volbrain,$@)) ; \
	orig_mask=$$(find "$$orig_mask_dir" -name '*crisp*') ; \
	flirt -interp nearestneighbour \
	      -in "$$orig_mask" \
	      -ref "$(subst volbrain-mask,feat.feat/example_func,$@)" \
	      -applyxfm \
	      -init "$<" \
	      -out "$@"

# TODO: feat.feat/ dirs are created as soon as feat starts. this confuses
#        dependents, which actually need it to finish creating all files

# %/reg/highres2example_func.mat : % ;

# %/example_func.nii.gz : % ;

################################################################################
# FSL FEAT preprocessing
################################################################################

.PHONY : feat_prepro
feat_prepro : $(addsuffix //feat.feat, $(addprefix $(DATA_DIR)/feat/, $(IDS)))
	@echo

$(DATA_DIR)/feat/%/feat.feat : $(DATA_DIR)/pymvpa/%/concat.nii.gz
	@echo "preprocessing with FEAT (high-pass, slice timing, motion correction and T1 registration) and saving into $@"
	@featdir=$(subst concat.nii.gz,,$(subst pymvpa,feat,$<)) ; mkdir -p "$$featdir" ; \
	t1dir=$(subst concat.nii.gz,,$(subst pymvpa,volbrain,$<)) ; \
	t1=$$(find "$$t1dir" -name '*_brain.nii.gz') ; \
	nvols=$$(fslnvols "$<") ; \
	sed "s|MVPA_OUTPUTDIR|$$(pwd)/$${featdir}/feat| ; s|MVPA_FEAT_FILES|$$(pwd)/$<| ; s|MVPA_HIGHRES_FILES|$$(pwd)/$${t1}| ; s|MVPA_NVOLS|$$nvols|" \
	    "$(SRC_DIR)/feat/1.prepro-design.fsf" > "$${featdir}/design.fsf" ; \
	feat "$${featdir}/design.fsf"

.PHONY : concatenate
concatenate : $(addsuffix /concat.nii.gz, $(addprefix $(DATA_DIR)/pymvpa/, $(IDS)))
	@echo

$(DATA_DIR)/pymvpa/%/concat.nii.gz : $(DATA_DIR)/bids/sub-%/func
	@echo 'concatenating fMRI series into $@'
	@dir=$@ ; mkdir -p "$${dir%/*}" ; \
	fmris=($$(find "$<" -type f -name '*bold.nii.gz' | grep 'task-$(TASKNAME)' | sort --version-sort)) ; \
	python2 "$(SRC_DIR)/pymvpa/concatenate.py" "$@" $${fmris[@]} > /dev/null 2>&1
	@fsledithd "$@" "$(SRC_DIR)/fix-nifti-units.sh"

################################################################################
# volbrain-related rules
################################################################################

# FIXME: requires manual upload/download from https://volbrain.upv.es
#        Reimplement with freesurfer?

.PHONY : t1w_brain_extraction
t1w_brain_extraction : $(VOLBRAIN_IMAGES:.nii=_brain.nii.gz)
	@echo

$(DATA_DIR)/volbrain/%_brain.nii.gz : $(DATA_DIR)/volbrain/%.nii
	@echo 'extracting brain to $@'
	@fslmaths $< -mul "$(subst n_,mask_n_,$<)" "$@"

.PHONY : volbrain_unzip
volbrain_unzip : $(VOLBRAIN_ZIPS:.zip=.volbrain)
	@echo

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
# convert DICOMs to Nifti
################################################################################

.PHONY : nifti
nifti : images # $(DICOMS:DICOM=nifti.nii.gz)
	@echo

%nifti.nii.gz : %DICOM
	@echo 'building $@'
	@dcm2niix -f 'nifti' -g y -i y -t n -z y "$</.." > /dev/null

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
