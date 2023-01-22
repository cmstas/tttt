#!/bin/bash

declare -r chkdir="$1"

echo "Checking the integrity of ROOT files in ${chkdir}:"

for f in $(find ${chkdir} -name '*.root'); do
  CheckFileIntegrity ${f}
  INTEGRITY_STATUS=$?
  if [[ ${INTEGRITY_STATUS} -ne 0 ]]; then
    echo "${f} failed."
    rm ${f}
  else
    echo "${f} succeeded."
  fi
done

