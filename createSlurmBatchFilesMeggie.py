    
import sys
import os
import json
import pathlib
import shutil
import subprocess
import glob        

#####################################################
# Command example:
#
# (For LIGGGHTS) "mpirun -np $SLURM_NTASKS ~/code/liggghts/src/lmp_auto -in input.in"
# 
#
######################################################
meggie_node_limit = 20
if __name__ == "__main__":       
    try:
        targetDirs = sys.argv[1]
        numOfProcessors = sys.argv[2]
        command = sys.argv[3]
        print("targetDirs = {} ".format(targetDirs))
        print("numOfProcessors = {} ".format(numOfProcessors))
        print("command = {} ".format(command))
    except:
        raise Exception("no targetDirs and numOfProcessors specify as the 1st and 2nd argument")

    for targetDir in glob.glob(targetDirs): 
        print(targetDir)
        folder_path = pathlib.Path(targetDir)
        folder_name = folder_path.name
        sbatchName = folder_name
      
        target_path = targetDir + "/" + sbatchName + ".sbatch"
        with open(target_path, "w") as target_file:
            target_content = []
            target_content.append("#! /bin/bash\n")
            target_content.append("#SBATCH --nodes={}\n".format(min(1,numOfProcessors/meggie_node_limit)))
            target_content.append("#SBATCH --ntasks-per-node={}\n".format(max(meggie_node_limit,numOfProcessors)))
            target_content.append("#SBATCH --cpus-per-task=1\n")
            target_content.append("#SBATCH --time=23:00:00\n")
            target_content.append("#SBATCH --mail-user=wing.ku@fau.edu\n")
            target_content.append("#SBATCH --output=test.out\n")
            target_content.append("#SBATCH --dependency=singleton\n")
            target_content.append(command)
            target_file.writelines(target_content)
            print("write file: {}".format(target_path))