#!/bin/bash
#
# Parse the command line inputs
if [ $# -ne 1 ]; then
    echo "$0: search directory is required"
    echo "     example: $0 /data/raw/dirname"
    exit 1
fi
SEARCHDIR=$1
echo "searching: $SEARCHDIR"
for file in ${SEARCHDIR}/*; do
  if [ -f "$file" ]; then
    echo "Found: ${file}"
    [ -n "$(tail -c1 "$file")" ] && printf '\r\n' >>$file
  fi
  done
