#!/bin/bash -l

echo
echo "Preparing to submit dependent jobs..."
echo

# We will use this variable to keep track of the last job we submitted
latest_id=0

command="sbatch run.sh" # the first job, #0
latest_id=$($command | awk ' { print $4 }')

# For each step
for i in `seq 16`; do # bc already ran it by executing it (ln 10)
    echo "Submitting job $i depending from job $latest_id"
    command="sbatch --dependency=afterany:$latest_id run.sh"
    latest_id=$($command | awk ' { print $4 }')
done

echo
echo "Done submitting dependent jobs!"

# chmod +x submit_dependent_jobs.sh (for the first time only to make it executable)
#  ./submit_dependent_jobs.sh (run only once for the following sims)
