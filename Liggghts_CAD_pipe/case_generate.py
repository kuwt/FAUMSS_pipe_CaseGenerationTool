import sys
import os
import json
import pathlib
import shutil
import subprocess
import glob
import random
import math
import numpy as np
import subprocess

sys.path.append(os.path.abspath("../"))
import loadJSONPara
import helixWallTextureGenerator
import defaultWallTextureGenerator
import getA
import mathUtility
import liggghtsPipeUtility as lpu

def AddLiggghtsVariablePipeGeometry(variableContent):
    simulationBoxExtendFactor = 3.0
    simulationBoxEpsilon = 2e-15
    variableContent.append( "##### Geometry definition #######\n")
    variableContent.append("variable xmin equal {}\n".format(-pipeRadius * simulationBoxExtendFactor))
    variableContent.append("variable xmax equal {}\n".format(pipeRadius * simulationBoxExtendFactor))
    variableContent.append("variable ymin equal {}\n".format(-pipeRadius * simulationBoxExtendFactor))
    variableContent.append("variable ymax equal {}\n".format(pipeRadius * simulationBoxExtendFactor))
    variableContent.append("variable zmin equal {}\n".format(-simulationBoxEpsilon))
    variableContent.append("variable zmax equal {}\n".format(pipeLength+simulationBoxEpsilon))
    variableContent.append("variable pipelength equal {}\n".format(pipeLength))
    variableContent.append("variable piperadius equal {}\n".format(pipeRadius))
    variableContent.append( "\n")
    return variableContent
        
def AddLiggghtsVariableParticleDef(variableContent):
    variableContent.append( "##### particle definition #######\n")
    variableContent.append( "variable nparticletype equal {} #{} -> Pipe, #{} -> Helix, #{}..#n ->Particles\n".format(nparticletype,pipeSpecies,textureSpecies,particleSpecies))

    DensityArray = particleDensity
    RadiusArray = particleRadius
    kn_specified_Matrix,kt_specified_Matrix,gamman_specified_Matrix,gammat_specified_Matrix = \
        lpu.computeForceCoefficient(nparticletype, pipeSpecies, textureSpecies, youngModulus,poissionRatio, coefficientOfRestitution, DensityArray, RadiusArray)
    variableContent.append( "##### particle properties kn_specified #######\n")
    for i in range(nparticletype):
        for j in range(nparticletype): 
            variableContent.append("variable kn{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),kn_specified_Matrix[i][j]))
    variableContent.append( "##### particle properties kt_specified #######\n")
    for i in range(nparticletype):
        for j in range(nparticletype): 
            variableContent.append("variable kt{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),kt_specified_Matrix[i][j]))
    variableContent.append( "##### particle properties gamman_specified_Matrix #######\n")
    for i in range(nparticletype):
        for j in range(nparticletype): 
            variableContent.append("variable gamman{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),gamman_specified_Matrix[i][j]))
    variableContent.append( "##### particle properties gammat_specified_Matrix #######\n")
    for i in range(nparticletype):
        for j in range(nparticletype): 
            variableContent.append("variable gammat{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),gammat_specified_Matrix[i][j]))
    variableContent.append( "##### particle properties sliding friction #######\n")
    for i in range(nparticletype):
        for j in range(nparticletype): 
            if i == pipeSpecies or j == pipeSpecies or i == textureSpecies or j == textureSpecies:
                variableContent.append("variable sfc{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),particleWallfrictionCoeffient))
            else:
                variableContent.append("variable sfc{}{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),lpu.cvrtSpeciesID_PythonToLiggghts(j),particleFrictionCoeffient))

    for i in range(nparticletype):
        if i >= particleSpecies:
            variableContent.append("variable particleradii{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),particleRadius[i]))
            variableContent.append("variable particlefraction{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),particleSolidFrac[i]/sum(particleSolidFrac)))
            variableContent.append("variable particledensity{} equal {}\n".format(lpu.cvrtSpeciesID_PythonToLiggghts(i),particleDensity[i]))
    return variableContent

def AddLiggghtsVariableParticleInsertion(variableContent):
    particleNum = 0
    pipeVolume = np.pi * (pipeRadius **2 ) * pipeLength
    for radius,solidFrac in zip(particleRadius,particleSolidFrac):
        particleVolume = mathUtility.computeSphereVolume(radius)
        particleNum += solidFrac * pipeVolume/particleVolume
    particleNum = math.ceil(particleNum)
    print("Total particleNum : {}".format(particleNum))
    variableContent.append("variable realparticleinsert equal {}\n".format(particleNum))
    variableContent.append("variable maxparticleinsert equal {}\n".format(particleNum*2))
    variableContent.append("variable insertionparticlerate equal {}\n".format(math.ceil(particleNum/(simulationTime * 0.05)))) #using only 1/20 simulation time to insert particles
    variableContent.append( "\n")
    return variableContent

def roundToTen(value):
    return int(round(value,-1))

def AddLiggghtsTimeStepping(variableContent):
    variableContent.append( "##### time stepping #######\n")
    variableContent.append("variable dt equal {}\n".format(timeStepSize))
    variableContent.append("variable dumpstep equal {}\n".format(roundToTen(dumpTimeInteval/timeStepSize)))
    variableContent.append("variable Simulationstep equal {}\n".format(roundToTen(simulationTime/timeStepSize)))
    variableContent.append("variable detailDumpstep equal {}\n".format(roundToTen(detailDumpTimeInteval/timeStepSize)))
    variableContent.append("variable detailSimulationstep equal {}\n".format(roundToTen(detailSimulationTime/timeStepSize)))
    variableContent.append( "\n")
    return variableContent

def AddLiggghtsSeed(variableContent):
    variableContent.append( "##### seed #######\n")
    for i in range(numOfSeedsRequiredPerCase):
        variableContent.append("variable theseed{} equal {}\n".format(i,primeNumberDataBase[trialid * numOfSeedsRequiredPerCase + i]))   
    variableContent.append( "\n")
    return variableContent

if __name__ == "__main__":
    ###########################################################################
    # Read parameters                             
    ######################################################################################
    # read from json for input parameters
    try:
        json_input_file_path = sys.argv[1]
        print("load json path = {} ".format(json_input_file_path))
    except:
        raise Exception("no json path specify as the first argument")

    for json_file_path in glob.glob(json_input_file_path ):
        print("json path = {} ".format(json_file_path))
        numOfTrial = loadJSONPara.readwithdefault(json_file_path,"numOfTrial",1)
        timeStepSize = loadJSONPara.read(json_file_path,"timeStepSize")
        dumpTimeInteval = loadJSONPara.read(json_file_path,"dumpTimeInteval") 
        simulationTime = loadJSONPara.read(json_file_path,"simulationTime")
        detailSimulationTime = loadJSONPara.read(json_file_path,"detailSimulationTime")
        detailDumpTimeInteval = loadJSONPara.read(json_file_path,"detailDumpTimeInteval")
        pipeRadius = loadJSONPara.read(json_file_path,"pipeRadius") 
        pipeLength = loadJSONPara.read(json_file_path,"pipeLength") 
        particleFrictionCoeffient = loadJSONPara.read(json_file_path,"particleFrictionCoeffient") 
        particleWallfrictionCoeffient = loadJSONPara.read(json_file_path,"particleWallfrictionCoeffient") 
        youngModulus = loadJSONPara.read(json_file_path,"youngModulus") 
        poissionRatio = loadJSONPara.read(json_file_path,"poissionRatio") 
        coefficientOfRestitution = loadJSONPara.read(json_file_path,"coefficientOfRestitution") 
        nparticletype = loadJSONPara.read(json_file_path,"nparticletype") 
        pipeSpecies = loadJSONPara.read(json_file_path,"pipeSpecies") 
        pipeSpecies = cvrtSpeciesID_LiggghtsToPython(pipeSpecies)
        textureSpecies = loadJSONPara.readwithdefault(json_file_path,"textureSpecies",2) 
        textureSpecies = cvrtSpeciesID_LiggghtsToPython(textureSpecies)
        particleSpecies = loadJSONPara.read(json_file_path,"particleSpecies") 
        particleSpecies = cvrtSpeciesID_LiggghtsToPython(particleSpecies)

        particleDensity = loadJSONPara.read(json_file_path,"particleDensity") 
        particleRadius = loadJSONPara.read(json_file_path,"particleRadius") 
        particleSolidFrac = loadJSONPara.read(json_file_path,"particleSolidFrac") 

        salomePath = loadJSONPara.readwithdefault(json_file_path,"salomePath","/local_disk/tools/SALOME-9.13.0-native-UB20.04-SRC/salome")
        CADSurfaceMeshSize = loadJSONPara.read(json_file_path,"CADSurfaceMeshSize")
        helixSpineAmplitude = loadJSONPara.read(json_file_path,"helixSpineAmplitude")
        helixSpinePeriod = loadJSONPara.read(json_file_path,"helixSpinePeriod")
        numOfProcessor = loadJSONPara.readwithdefault(json_file_path,"numOfProcessor",4)

        ###########################################################################
        # Create random prime number base
        ######################################################################################
        numOfSeedsRequiredPerCase = 4
        primeNumberDataBase = generate_prime_list(numOfTrial * numOfSeedsRequiredPerCase, 31, 500000)
        print("primeNumberDataBase = {}".format(primeNumberDataBase))
        
        ###########################################################################
        # Create folder
        #####################################################################################
        for trialid in range(numOfTrial):
            output_directory = pathlib.Path(json_file_path).stem + "_{}".format(trialid)
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)

            ###########################################################################
            # Create init.in
            #####################################################################################
            target_content = []
            target_content = AddLiggghtsVariablePipeGeometry(target_content)
            target_content = AddLiggghtsVariableParticleDef(target_content)
            target_path = output_directory + "/init.in"
            with open(target_path, "w") as target_file:
                with open("init_template.in") as templ_file:
                    templ_content = templ_file.readlines()
                    target_content = target_content + templ_content
                    target_file.writelines(target_content)
                    print("write file: {}".format(target_path))
            ###########################################################################
            # Create run.in
            #####################################################################################
            target_content = []
            target_content = AddLiggghtsVariablePipeGeometry(target_content)
            target_content = AddLiggghtsVariableParticleDef(target_content)
            target_content = AddLiggghtsTimeStepping(target_content)
            target_content = AddLiggghtsSeed(target_content)
            target_path = output_directory + "/run.in"
            
            with open(target_path, "w") as target_file:
                with open("run_template.in") as templ_file:
                    templ_content = templ_file.readlines()
                    target_content = target_content + templ_content
                    target_file.writelines(target_content)
                    print("write file: {}".format(target_path))

            ###########################################################################
            # Create detailrun.in
            #####################################################################################
            target_content = []
            target_content = AddLiggghtsVariablePipeGeometry(target_content)
            target_content = AddLiggghtsVariableParticleDef(target_content)
            target_content = AddLiggghtsTimeStepping(target_content)
            target_content = AddLiggghtsSeed(target_content)

            target_path = output_directory + "/detailrun.in" 
            with open(target_path, "w") as target_file:
                with open("detailrun_template.in") as templ_file:
                    templ_content = templ_file.readlines()
                    target_content = target_content + templ_content
                    target_file.writelines(target_content)
                    print("write file: {}".format(target_path))
                    
            ###########################################################################
            # Create CAD
            #####################################################################################
            json_template_dictionary = {
                "helixSpineAmplitude":helixSpineAmplitude,
                "pipeLength":pipeLength,
                "helixSpinePeriod":helixSpinePeriod,
                "pipeRadius":pipeRadius,
                "surfaceMeshSize":CADSurfaceMeshSize,
                "output_directory" : output_directory
            }
            json_object = json.dumps(json_template_dictionary, indent=4)
            tmpCADcfgpath = output_directory + "/tmpCADcfg.json"
            with open(tmpCADcfgpath, "w") as outfile:
                outfile.write(json_object)
            print("running Salome")
            
            result = subprocess.run([salomePath,"-t","../constructHelixPipeSurfaceMeshSalome.py","args:{}".format(tmpCADcfgpath)])
            if result.returncode != 0:
                raise("running Salome fail")
                
            import vtkconvert
            pipeVolumeBinary_path = output_directory + '/pipeVolumeBinary.vtk'
            pipeVolumeAscii_path = output_directory + '/pipeVolumeAscii.vtk'
            vtkconvert.convert_vtk_binary_to_ascii(pipeVolumeBinary_path, pipeVolumeAscii_path)

             ###########################################################################
            # Create initNrun.sh
            #####################################################################################
            variables = {
            "numOfProcessor": numOfProcessor
            }
            target_path = output_directory + "/" + "initNrun.sh"
            with open(target_path, "w") as target_file:
                with open("initNrun_template.sh") as templ_file:
                    target_content = templ_file.read()
                    for key, val in variables.items():
                        target_content = target_content.replace(f"{{{key}}}", str(val))
                    target_file.write(target_content)
                    print("write file: {}".format(target_path))
        
