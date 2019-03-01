#!/usr/bin/env bash

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

# download all scan types (DICOM format only) from XNAT for a
# given list $FILE with "PATIENT_ID EXPERIMENT_ID" pairs, using the REST API.
# one ${patient_id}.zip per patient.

readonly FILE='xnat-ids.csv'
readonly DELIMITER=' '
readonly XNAT_HOST='http://172.24.80.68:8080'

printf 'xnat username: '
read -r username
printf 'xnat password: '
read -rs password

while read -pr line; do
    pid=$(cut -d "$DELIMITER" -f 1 <<<"$line")
    expid=$(cut -d "$DELIMITER" -f 2 <<<"$line")
    url="${XNAT_HOST}/REST/archive/experiments/${expid}/scans/ALL/resources/DICOM/files?format=zip"

    wget -O "$pid".zip \
         "--http-user=${username}" \
         "--http-password=${password}" \
         --auth-no-challenge \
         "$url" || exit 0 # abort on error
done < "$FILE"
