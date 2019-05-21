#!/bin/bash
# Submit all jobs based on list of input parameter files

FILES=arg_files/*.txt

for f in $FILES
do
    qsub run_general_bias_sim.sh -F $f
done
