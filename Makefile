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
DICOMS = $(shell find $(DATA_DIR)/xnat/images/ -type d -name DICOM | \
                 grep -E '(fMRI_GazeCueing|FSPGR|T2)' | grep -v '00-PU' | sort)
VOLBRAIN_ZIPS = $(shell find $(DATA_DIR)/volbrain/ -type f -name '*.zip')
VOLBRAIN_IMAGES = $(shell find data/volbrain/ -type f -name 'native_n_*')
FMRI_NIFTIS = $(shell find $(DATA_DIR)/xnat/images/ -type f \
                      -name '*nifti.nii.gz' | grep GazeCueing | sort)
FEAT_NIFTIS = $(subst /resources/nifti.nii.gz,.feat, \
                      $(subst xnat/images,feat,$(FMRI_NIFTIS)))

.PHONY : build all
all : build
build : volbrain_tree volbrain_unzip concatenate_runs

################################################################################
# poststats
################################################################################

.PHONY : poststats
poststats : pymvpa
	@echo 'running result statistics'
	@Rscript -e 'source("$(SRC_DIR)/poststats/poststats.R")'

################################################################################
# pyMVPA rules
################################################################################

#.PHONY : pymvpa
#pymvpa : $(subst data,out,$(FEAT_NIFTIS:%=%/hap_vs_sad-weights-nn.nii.gz)) ;

# FIXME: missing RHS
%/hap_vs_sad-weights-nn.nii.gz :
	@outdir=$(subst hap_vs_sad-weights-nn.nii.gz,,$@) ; \
	nifti=$@ ; \
	nifti=$${nifti/out/data} ; \
	nifti=$${nifti/hap_vs_sad-weights-nn.nii.gz/filtered_func_data_brain.nii.gz} ; \
	eprime=$$outdir ; eprime=$${eprime/feat/eprime} ; eprime=$${eprime/scans*_/} ; \
	eprime=$$(echo "$$eprime" | sed -E 's!(1|2|3)(rep)?.feat/!\1.txt.csv!') ; \
	if [[ -e "$$eprime" ]]; then \
	    echo "running pyMVPA for $$outdir" ; \
	    mkdir -p "$$outdir" ; \
	    if [ $$((RANDOM % 2)) -eq 0 ]; then \
	        fsl_sub python2 "$(SRC_DIR)/pymvpa/pymvpa.py" "$$eprime" "$$nifti" "$$outdir" > /dev/null 2>&1 ; \
	    else \
	        python2 "$(SRC_DIR)/pymvpa/pymvpa.py" "$$eprime" "$$nifti" "$$outdir" > /dev/null 2>&1 ; \
	    fi ; \
	fi

.PHONY : pymvpa
pymvpa : $(addprefix $(BUILD_DIR)/pymvpa/, $(IDS))

$(BUILD_DIR)/pymvpa/% : $(DATA_DIR)/pymvpa/%
	@echo "running pyMVPA for $<" ; \
	mkdir -p "$@" ; \
	fsl_sub python2 "$(SRC_DIR)/pymvpa/pymvpa.py" "$</events.csv" \
	        "$</concat.nii.gz" "$@" > /dev/null 2>&1

.PHONY : concatenate_runs
concatenate_runs : $(addprefix $(DATA_DIR)/pymvpa/, $(IDS))
	@echo

$(DATA_DIR)/pymvpa/% : $(DATA_DIR)/feat/%
	@mkdir -p "$@"
	@echo 'concatenating runs to $@'
	@events=($$(find "$(subst $(DATA_DIR)/feat,$(BUILD_DIR)/eprime,$<)" \
	                 -type f -name '*.csv' | sort)) ; \
	fmris=($$(find "$<" -type f -name 'filtered_func_data_brain.nii.gz' | \
	          sed 's/GazeCueing_/ /' | sort -k 2 | sed 's/ /GazeCueing_/')) ; \
	python2 "$(SRC_DIR)/pymvpa/concatenate-events.py" "$@/events.csv" \
	        $${events[@]} > /dev/null 2>&1 & \
	python2 "$(SRC_DIR)/pymvpa/prepro.py" "$@/concat.nii.gz" \
	        $${fmris[@]} > /dev/null 2>&1

################################################################################
# post-FEAT brain extraction-related rules
################################################################################

.PHONY : feat-brains
feat-brains : $(FEAT_NIFTIS:%=%/filtered_func_data_brain.nii.gz)
	@echo

%/filtered_func_data_brain.nii.gz : %/volbrain-mask.nii.gz %/filtered_func_data.nii.gz
	@echo 'extracting EPI brain to $@'
	@fslmaths "$(subst _brain,,$@)" -mul "$<" "$@"

%/filtered_func_data.nii.gz : % ;

.PHONY : feat-masks
feat-masks : $(FEAT_NIFTIS:%=%/volbrain-mask.nii.gz)
	@echo

#TODO: add prerequisite to _brain.nii.gz images, so as to complete dependency graph.
#      deduplicate per-scan masks, since they are all the same for a given subj.
%/volbrain-mask.nii.gz : %/reg/highres2example_func.mat %/example_func.nii.gz
	@echo 'creating EPI mask $@'
	@orig_mask_dir=$(subst feat,volbrain,$@) ; \
	orig_mask=$$(find "$${orig_mask_dir/scans*/}" -name '*mask*') ; \
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

.PHONY : feat-prepro
#TODO: add prerequisite to _brain.nii.gz images, so as to complete dependency graph
feat-prepro : $(addsuffix .feat,$(subst /resources/nifti.nii.gz,,$(subst xnat/images,feat,$(FMRI_NIFTIS))))
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
	@fslmaths $< -mul "$(subst _n_,_mask_n_,$<)" "$@"

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
	        printf 'downloading DICOMS for subject %d\n\n' "$$i" ; \
	        $(SRC_DIR)/xnat/$@/xnat-download.sh $$(grep "^$$i " $<) \
	                                            "$(DATA_DIR)/xnat/$@" || \
	            rm "$(DATA_DIR)/xnat/$@/$${i}.zip" ; \
	    fi ; \
	    printf "unzip-ing DICOMs for subject %d\n\n" "$$i" ; \
	    unzip -q "$(DATA_DIR)/xnat/$@/$${i}.zip" \
	          -d "$(DATA_DIR)/xnat/$@/" && rm "$(DATA_DIR)/xnat/$@/$${i}.zip" ; \
	done ;
# delete duplicate T1 directories. we'll do smoothing and normalisation manually
	@find "$(DATA_DIR)/xnat/$@/" -type d -name '*00-PU*' -prune \
	         -exec bash -c 'echo "deleting derived images" {} ; rm -r {}' \;
# delete fMRI sequences with missing volumes (<260)
	@rm -rf "$(DATA_DIR)/xnat/$@/517/scans/5-fMRI_GazeCueing_1"
# TODO: rescue vols, unique sequence
	@rm -rf "$(DATA_DIR)/xnat/$@/812/scans/8-fMRI_GazeCueing_2"
	@echo
# delete fMRI sequences without corresponding emprime events file
	@rm -rf "$(DATA_DIR)/xnat/$@/535/scans/6-fMRI_GazeCueing_3"

################################################################################
# eprime events: find files, clean and convert into design matrices
################################################################################

.PHONY : eprime
eprime : $(IDS_FILE) $(DATA_DIR)/$@/
	@printf 'building design matrices from eprime event lists\n\n'
# FIXME: is it feasible to write a generic per-file or per-subject rule?
	@rm -rf "$(BUILD_DIR)/$@"
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
# homogenize disparate file names into {1,2,3}.txt
	@find $(BUILD_DIR)/$@ -type f -name '*.txt' -exec bash -c \
	    'mv "{}" $$(echo {} | sed -En "s/(.*)(\/)(.*)(_)(A|B|C)(.*)(.txt)/\1\2\5\7/ ; s/A/1/p ; s/B/2/p ; s/C/3/p")' \;
# UTF-16 -> UTF-8, CRLF -> LF
	@find $(BUILD_DIR)/$@ -type f -name '*.txt' -exec bash -c \
	    'iconv -f UTF-16 -t UTF-8 "{}" | tr -d "\r" > "{}.new" && mv "{}.new" "{}"' \;
# more manual fixes (file contents)
	@sed -i 's/Sex: male/Sex: female/' $(BUILD_DIR)/$@/559/3.txt
	@sed -i 's/Age: 0/Age: 35/' $(BUILD_DIR)/$@/575/3.txt
	@sed -i 's/Age: 22/Age: 23/' $(BUILD_DIR)/$@/590/{2,3}.txt
	@sed -i 's/Age: 0/Age: 31/' $(BUILD_DIR)/$@/672/{1,2}.txt
	@sed -i 's/Sex: male/Sex: female/' $(BUILD_DIR)/$@/678/2.txt
	@sed -i 's/Age: 23/Age: 22/' $(BUILD_DIR)/$@/678/3.txt
	@sed -i 's/Age: 28/Age: 27/' $(BUILD_DIR)/$@/696/3.txt
# delete files without corresponding fMRI sequence
	@rm $(BUILD_DIR)/$@/518/{1,2}.txt
	@rm $(BUILD_DIR)/$@/812/2.txt
# eprime event list -> pyMVPA sample attribute matrix
	@find $(BUILD_DIR)/$@ -type f -name '*.txt' -exec bash -c \
	    'awk -f "$(SRC_DIR)/eprime/eprime-to-csv.awk" -- "{}" > "{}.csv"' \;

################################################################################
# xnat's subject metadata DB: valid participants will be selected from there
################################################################################

$(IDS_FILE) : $(SRC_DIR)/xnat/subject_metadata/extract.R
	@printf '\building subject IDs list\n\n'
	@mkdir -p "$(BUILD_DIR)/xnat/subject_metadata"
	Rscript -e 'source("$(SRC_DIR)/xnat/subject_metadata/extract.R")'
