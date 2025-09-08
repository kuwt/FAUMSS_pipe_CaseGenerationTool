#!/bin/bash

. ~/.bashrc
source $CFDEM_PROJECT_DIR/src/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath="$casePath"
headerText="run_liggghts_init_DEM"
logfileName="log_$headerText"
solverName="init.in"
machineFileName="none"   # yourMachinefileName | none
#--------------------------------------------------------------------------------#

#- call function to run DEM case
parDEMrun $logpath $logfileName $casePath $headerText $solverName 1 $machineFileName

#--------------------------------------------------------------------------------#
#- define variables
headerText="run_parallel_cfdemSolverPiso_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs={numOfProcessor}
#--------------------------------------------------------------------------------#

#- liggghts, delete the last timestep since it might be corrupted, and there is no good ways to handle corruption.
echo "checking liggghts restart file"
cd $casePath/DEM/restart
# Get the file with the highest number
file=$(ls | grep -E '[0-9]+' | while read f; do
    num=$(echo "$f" | grep -oE '[0-9]+' | tail -n1)
    echo "$num $f"
done | sort -k1,1n | tail -n1 | awk '{print $2}')

# Extract the number
num=$(echo "$file" | grep -oE '[0-9]+' | tail -n1)

# If the number > 0, delete the file
if [ -n "$num" ] && [ "$num" -gt 0 ]; then
    echo "Deleting: $file"
    rm -- "$file"
else
    echo "No file deleted."
fi

#- openfoam, delete the last timestep too..
echo "checking openfoam restart file"
cd $casePath
base="CFD/processor0"
# Find folder with greatest number
folder=$(ls -d "$base"/*/ 2>/dev/null | grep -Eo '[0-9]+/?$' | tr -d / | sort -n | tail -n1)

# Check if we got a number
if [ -n "$folder" ] && [ "$folder" -gt 0 ]; then
    echo "Deleting folder number: $folder"
    
    # Loop through all processor directories and delete matching folder
    for d in CFD/processor*/; do
        if [ -d "$d/$folder" ]; then
            echo "Deleting: $d/$folder"
            rm -rf -- "$d/$folder"
        fi
    done
else
    echo "No folder deleted."
fi


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

