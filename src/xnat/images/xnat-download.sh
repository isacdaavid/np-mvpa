#!/usr/bin/env bash

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

# download all scan types (DICOM format only) from XNAT for a given
# PATIENT_ID EXPERIMENT_ID pair, using the REST API. one ${patient_id}.zip
# per patient, save at third argument.

readonly PID=$1
readonly EXPID=$2
readonly OUTDIR=$3
readonly XNAT_HOST='http://172.24.80.68:8080'

printf 'xnat username: '
read -r username
printf 'xnat password: '
read -rs password


readonly URL="${XNAT_HOST}/REST/archive/experiments/${EXPID}/scans/ALL/resources/DICOM/files?format=zip"

wget -O "${OUTDIR}/${PID}".zip \
         "--http-user=${username}" \
         "--http-password=${password}" \
         --auth-no-challenge \
         "$URL"
