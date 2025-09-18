#!/bin/bash

. ~/.bashrc
source $CFDEM_PROJECT_DIR/etc/functions.sh
#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath="$casePath"
headerText="run_liggghts_init_DEM"
logfileName="log_$headerText"
solverName="init.in"
machineFileName="none"   # yourMachinefileName | none
#--------------------------------------------------------------------------------#
echo "checking liggghts restart file"
cd $casePath/DEM/restart
if find . -maxdepth 1 -type f -regex '.*/.*[0-9].*' | grep -q .; then
    echo "skip dem init"
else
    #- call function to run DEM case
    cd $casePath
    parDEMrun $logpath $logfileName $casePath $headerText $solverName 1 $machineFileName
fi
#--------------------------------------------------------------------------------#
#- define variables
headerText="run_parallel_cfdemSolverPiso_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs={numOfProcessor}
#--------------------------------------------------------------------------------#
#- liggghts, delete the last timestep since it might be corrupted, and there is no good ways to handle corruption.
echo "delete last time step"
cd $casePath/DEM/restart
if [ $? -eq 0 ]; then
	max=$(ls 2>/dev/null | grep -Eo '[0-9]+(\.[0-9]+)?' | sort -n | tail -n 1)
	echo "The last time step is:" $max
	if [[ -z "$max" ]]; then
	    echo "No files with numbers found."
	    exit 0
	fi
fi
# Find the file(s) that contain this max number
file_to_delete=$(ls | grep -F "$max" | head -n 1)

# Compare number > 0
if (( $(echo "$max > 0" | bc -l) )); then
    echo "Deleting file: $file_to_delete (because $max > 0)"
    rm -f "$file_to_delete"
else
    echo "Max number $max is not greater than 0. No file deleted."
fi

#- openfoam, delete the last timestep too..
echo "checking openfoam restart file"
target_base="$casePath/CFD/"
for proc_dir in "$target_base"/processor*/; do
    echo "Checking $proc_dir"
    cd "$proc_dir" || continue

    for folder in */; do
        num=$(echo "$folder" | grep -Eo '[0-9]+(\.[0-9]+)?')
        if [[ -n "$num" ]]; then
            if (( $(echo "$max > 0" | bc -l) )) && (( $(echo "$num >= $max" | bc -l) )); then
                echo "Deleting $proc_dir$folder (because $num >= $max)"
                rm -rf "$folder"
            fi
        fi
    done
done

#- liggghts restart file, not using liggghts built in directly since that only supports integer
echo "copying liggghts restart file"
cd $casePath/DEM/restart
max=$(ls 2>/dev/null | grep -Eo '[0-9]+(\.[0-9]+)?' | sort -n | tail -n 1)
if [[ -z "$max" ]]; then
    echo "No file names with numbers found."
    exit 1
fi
file_to_copy=$(ls | grep -F "$max" | head -n 1)
echo "file to copy: $file_to_copy (number: $max)"
cp "$file_to_copy" liggghts.restart
echo "copy : $file_to_copy to liggghts.restart"


# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

# check if decomposePar already done
if [ -d "$casePath/CFD/processor0" ]; then
    echo "decomposePar was done before"
else
    echo "decomposePar needed now"
    cd $casePath/CFD
    decomposePar
fi

# logging
rm $logpath/$logfileName
echo 2>&1 | tee -a /$logpath/$logfileName
echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
echo 2>&1 | tee -a $logpath/$logfileName
pwd 2>&1 | tee -a $logpath/$logfileName
echo 2>&1 | tee -a $logpath/$logfileName

# run
cd $casePath/CFD
mpirun -np $nrProcs $debugMode $solverName -parallel 2>&1 | tee -a $logpath/$logfileName

