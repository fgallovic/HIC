# Print shell input lines to standard error
set -v
#
# Synthesize final results calculated in directories 'dw' and 'dw2'
./KK
#
# Produce final set of output files
./analyze
#
echo "calculate3.sh done!"
