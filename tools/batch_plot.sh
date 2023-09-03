#!/bin/bash
# Read the cast numbers from a CSV file and feed them to the indicated plotting
# script.

# Parse the command line inputs
if [ $# -ne 3 ]; then
    echo "$0: required inputs are the parameter to plot, the full pathname of the"
    echo "calibration file (or \"none\") and the full pathname of the CSV file of"
    echo "cast numbers."
    echo "     example: $0 oxygen /path/to/calfile.csv /path/to/ssscc.csv"
    exit 1
fi

case ${1} in
  oxygen)
    SCRIPT="plot_rinko.py"
    ;;
  reft)
    SCRIPT="plot_reft.py"
    ;;
  *)
    echo "Unknown parameter. Exiting"
    exit 1
    ;;
esac

CALFILE=${2}
CSVFILE=${3}

# Set the tools directory path...
TOOLS="/Users/als026/code/ctdcal/tools/"

# Call the plotting script for each line in the CSV file...
while read ssscc; do
  python "${TOOLS}${SCRIPT}" "${CALFILE}" "${ssscc}"
done <"${CSVFILE}"
