import sys
import os
import json
import pathlib
import shutil
import subprocess
import glob
import math
import numpy as np

script_path = os.path.abspath(__file__)
print("Script path:", script_path)
script_dir = os.path.dirname(script_path)
print("Script directory:", script_dir)
sys.path.append(os.path.abspath(script_dir + "/../"))
import loadJSONPara
import helixWallTextureGenerator
import defaultWallTextureGenerator
import mathUtility
import fileUtility
import liggghtsPipeUtility as lpu

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
        pipeRadius = loadJSONPara.read(json_file_path,"pipeRadius") 
        pipeLength = loadJSONPara.read(json_file_path,"pipeLength") 
        particleFrictionCoeffient = loadJSONPara.read(json_file_path,"particleFrictionCoeffient") 
        particleWallfrictionCoeffient = loadJSONPara.read(json_file_path,"particleWallfrictionCoeffient") 

        youngModulus = loadJSONPara.read(json_file_path,"youngModulus") 
        poissionRatio = loadJSONPara.read(json_file_path,"poissionRatio") 
        coefficientOfRestitution = loadJSONPara.read(json_file_path,"coefficientOfRestitution") 
        characteristicVelocity= loadJSONPara.readwithdefault(json_file_path,"characteristicVelocity",1.0) 
        nparticletype = loadJSONPara.read(json_file_path,"nparticletype") 
        pipeSpecies = loadJSONPara.read(json_file_path,"pipeSpecies") 
        pipeSpecies = lpu.cvrtSpeciesID_LiggghtsToPython(pipeSpecies)
        textureSpecies = loadJSONPara.read(json_file_path,"textureSpecies") 
        textureSpecies = lpu.cvrtSpeciesID_LiggghtsToPython(textureSpecies)
        particleSpecies = loadJSONPara.read(json_file_path,"particleSpecies") 
        particleSpecies = lpu.cvrtSpeciesID_LiggghtsToPython(particleSpecies)

        particleDensity = loadJSONPara.read(json_file_path,"particleDensity") 
        particleRadius = loadJSONPara.read(json_file_path,"particleRadius") 
        particleSolidFrac = loadJSONPara.read(json_file_path,"particleSolidFrac") 

        particleInsertionDiskSize = loadJSONPara.read(json_file_path,"particleInsertionDiskSize") 
        particleInsertionDiskVolumeFraction = loadJSONPara.read(json_file_path,"particleInsertionDiskVolumeFraction") 
        particleInsertionVelocity = loadJSONPara.read(json_file_path,"particleInsertionVelocity") 

        salomePath = loadJSONPara.readwithdefault(json_file_path,"salomePath","/local_disk/tools/SALOME-9.13.0-native-UB20.04-SRC/salome")
        CADSurfaceMeshSize = loadJSONPara.read(json_file_path,"CADSurfaceMeshSize")
        CADSurfaceMeshSizeMax = loadJSONPara.readwithdefault(json_file_path,"CADSurfaceMeshSizeMax",CADSurfaceMeshSize*2)
        helixSpineAmplitude = loadJSONPara.read(json_file_path,"helixSpineAmplitude")
        helixSpinePeriod = loadJSONPara.read(json_file_path,"helixSpinePeriod")

        fluidVelocity = loadJSONPara.read(json_file_path,"fluidVelocity")
        kinematicViscosity = loadJSONPara.read(json_file_path,"kinematicViscosity")
        gravity = loadJSONPara.read(json_file_path,"gravity")
        fluidGridSizeXY = loadJSONPara.read(json_file_path,"fluidGridSizeXY")
        fluidGridSizeZ = loadJSONPara.read(json_file_path,"fluidGridSizeZ")

        numOfProcessor = loadJSONPara.readwithdefault(json_file_path,"numOfProcessor",4) 
        ###########################################################################
        # Create random prime number base
        ######################################################################################
        numOfSeedsRequiredPerCase = 4
        primeNumberDataBase = mathUtility.generate_prime_list(numOfTrial * numOfSeedsRequiredPerCase, 31, 500000)
        print("primeNumberDataBase = {}".format(primeNumberDataBase))
        
        ###########################################################################
        # tmp working directory
        #####################################################################################
        output_directory = pathlib.Path(json_file_path).stem
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        output_directory = os.path.abspath(output_directory)
        ###########################################################################
        # Create CFD
        #####################################################################################
        CFD_dir = output_directory + "/CFD/"
        if not os.path.exists(CFD_dir):
            os.mkdir(CFD_dir)
        ###########################################################################
        # Create CFD/0
        #####################################################################################
        CFD_0_dir = CFD_dir + "/0/"
        if not os.path.exists(CFD_0_dir):
            os.mkdir(CFD_0_dir)

        variables = {
        "fluidVelocity": fluidVelocity
        }
        
        fileUtility.copyfile("CFD_template/0_template/U", CFD_0_dir + "U", variables)
        fileUtility.copyfile("CFD_template/0_template/Us", CFD_0_dir + "Us")
        fileUtility.copyfile("CFD_template/0_template/p", CFD_0_dir + "p")
        fileUtility.copyfile("CFD_template/0_template/voidfraction", CFD_0_dir + "voidfraction")
        fileUtility.copyfile("CFD_template/0_template/Ksl", CFD_0_dir + "Ksl")
        fileUtility.copyfile("CFD_template/0_template/rho", CFD_0_dir + "rho")
        ###########################################################################
        # Create CFD/system
        #####################################################################################
        CFD_system_dir = CFD_dir + "/system/"
        if not os.path.exists(CFD_system_dir):
            os.mkdir(CFD_system_dir)

        extension = 0.1
        variables = {
        "XMIN":- (pipeRadius + helixSpineAmplitude) * (1.0+extension),     
        "XMAX":(pipeRadius + helixSpineAmplitude) * (1.0+extension),     
        "YMIN": -(pipeRadius + helixSpineAmplitude) * (1.0+extension),      
        "YMAX": (pipeRadius + helixSpineAmplitude) * (1.0+extension),           
        "ZMIN": - pipeLength * extension,  # the pipe extend from 0 to the +ve x direction
        "ZMAX": pipeLength * 1.1,
        "fluidGridSizeXY": fluidGridSizeXY,
        "fluidGridSizeZ": fluidGridSizeZ
        }
        fileUtility.copyfile("CFD_template/system_template/blockMeshDict", CFD_system_dir + "blockMeshDict", variables)
        fileUtility.copyfile("CFD_template/system_template/fvSchemes", CFD_system_dir + "fvSchemes")
        fileUtility.copyfile("CFD_template/system_template/fvSolution", CFD_system_dir + "fvSolution")
        variables = {
        "numOfProcessor":numOfProcessor
        }
        fileUtility.copyfile("CFD_template/system_template/decomposeParDict", CFD_system_dir + "decomposeParDict",variables)
        variables = {
        "simulationTime":simulationTime,
        "timeStepSize":timeStepSize,
        "dumpTimeInteval":dumpTimeInteval
        }
        fileUtility.copyfile("CFD_template/system_template/controlDict", CFD_system_dir + "controlDict",variables)
        fileUtility.copyfile("CFD_template/system_template/meshQualityDict", CFD_system_dir + "meshQualityDict")
        fileUtility.copyfile("CFD_template/system_template/surfaceFeaturesDict", CFD_system_dir + "surfaceFeaturesDict")
        fileUtility.copyfile("CFD_template/system_template/snappyHexMeshDict", CFD_system_dir + "snappyHexMeshDict")
        
        ###########################################################################
        # Create CFD/constant
        #####################################################################################
        CFD_constant_dir = CFD_dir + "/constant/"
        if not os.path.exists(CFD_constant_dir):
            os.mkdir(CFD_constant_dir)
        fileUtility.copyfile("CFD_template/constant_template/turbulenceProperties", CFD_constant_dir + "turbulenceProperties")
        variables = {
        "kinematicViscosity":kinematicViscosity
        }
        fileUtility.copyfile("CFD_template/constant_template/transportProperties", CFD_constant_dir + "transportProperties",variables)

        variables = {
        "gravity":gravity
        }
        fileUtility.copyfile("CFD_template/constant_template/g", CFD_constant_dir + "g",variables)
        fileUtility.copyfile("CFD_template/constant_template/couplingProperties", CFD_constant_dir + "couplingProperties")
        fileUtility.copyfile("CFD_template/constant_template/liggghtsCommands", CFD_constant_dir + "liggghtsCommands")

        ###########################################################################
        # Create CAD
        #####################################################################################
        CADDir = CFD_constant_dir + "/triSurface/" 
        if not os.path.exists(CADDir):
            os.mkdir(CADDir)
            json_template_dictionary = {
                "helixSpineAmplitude":helixSpineAmplitude,
                "pipeLength":pipeLength,
                "helixSpinePeriod":helixSpinePeriod,
                "pipeRadius":pipeRadius,
                "surfaceMeshSize":CADSurfaceMeshSize,
                "surfaceMeshSizeMax":CADSurfaceMeshSizeMax,
                "output_directory" : CADDir
            }
            json_object = json.dumps(json_template_dictionary, indent=4)
            tmpCADcfgpath = output_directory + "/tmpCADcfg.json"
            with open(tmpCADcfgpath, "w") as outfile:
                outfile.write(json_object)
            print("running Salome")
            result = subprocess.run([salomePath,"-t","../constructHelixPipeSurfaceMeshSalome.py","args:{}".format(tmpCADcfgpath)])
            if result.returncode != 0:
                raise("running Salome fail")

            fromSalomeMeshToOpenfoamMesh = """
            cd {tmpDir}
            cat inlet.stl outlet.stl pipewall.stl  > surfacemesh.stl
            cd ../..
            blockMesh | tee blockMesh.log
            surfaceFeatures | tee surfaceFeatures.log
            snappyHexMesh -overwrite| tee log.snappyHexMesh
            """.format(tmpDir=CADDir)

            result = subprocess.run(fromSalomeMeshToOpenfoamMesh,shell=True, executable="/bin/bash", check=True)
            if result.returncode != 0:
                raise("fromSalomeMeshToOpenfoamMesh fail")
        else:
            print("CADDir:{} already exist. Skip remaking CAD since this takes a long time!".format(CADDir))
        ###########################################################################
        # of trials
        #####################################################################################
        for trialid in range(numOfTrial):
            output_directory = pathlib.Path(json_file_path).stem + "_{}".format(trialid)
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)
            
            ###########################################################################
            # Create DEM 
            #####################################################################################
            DEM_dir = output_directory + "/DEM/"
            if not os.path.exists(DEM_dir):
                os.mkdir(DEM_dir)
            
            ###########################################################################
            # Create DEM/init.in
            #####################################################################################
            target_content = []

            def AddLiggghtsVariablePipeGeometry(variableContent):
                simulationBoxExtendFactor = 3.0
                simulationBoxEpsilon = pipeRadius * 0.00001
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
                    lpu.computeForceCoefficient(nparticletype, pipeSpecies, textureSpecies, youngModulus,poissionRatio, coefficientOfRestitution, characteristicVelocity, DensityArray, RadiusArray)
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
            target_content = AddLiggghtsVariablePipeGeometry(target_content)
            target_content = AddLiggghtsVariableParticleDef(target_content)
            fileUtility.copyfilePrepend("DEM_template/init_template.in" ,DEM_dir + "/init.in",target_content)

            ###########################################################################
            # Create DEM/run.in
            #####################################################################################
            target_content = []
            def AddLiggghtsSeed(variableContent):
                variableContent.append( "##### seed #######\n")
                for i in range(numOfSeedsRequiredPerCase):
                    variableContent.append("variable theseed{} equal {}\n".format(i,primeNumberDataBase[trialid * numOfSeedsRequiredPerCase + i]))   
                variableContent.append( "\n")
                return variableContent

            target_content = AddLiggghtsVariableParticleDef(target_content)
            target_content = AddLiggghtsSeed(target_content)

            target_content.append( "##### particle insertion #######\n")
            target_content.append("variable particleInsertionDiskSize equal {}\n".format(particleInsertionDiskSize))
            target_content.append("variable particleInsertionDiskVolumeFraction equal {}\n".format(particleInsertionDiskVolumeFraction))
            target_content.append("variable particleInsertionVelocity equal {}\n".format(particleInsertionVelocity))

            target_content.append( "##### time #######\n")
            target_content.append("variable dt equal {}\n".format(timeStepSize))
            target_content.append("variable dumpstep equal {}\n".format(int(dumpTimeInteval/timeStepSize)))

            target_content.append( "##### gravity #######\n")
            target_content.append("variable definedGravity equal {}\n".format(gravity))
            
            target_content.append( "variable numOfProcessor equal {}\n".format(numOfProcessor))

            target_content.append( "################### variable definition finished #####################\n\n\n")
            fileUtility.copyfilePrepend("DEM_template/run_template.in" ,DEM_dir + "/run.in",target_content)

            ###########################################################################
            # CFD
            #####################################################################################
            CFD_targetDir = output_directory + "/CFD/"
            fileUtility.my_copytree(CFD_dir,CFD_targetDir)

            ###########################################################################
            # Create initNrun.sh
            #####################################################################################
            variables = {
            "numOfProcessor": numOfProcessor
            }
            fileUtility.copyfile("initNrun_template.sh", output_directory + "/" + "initNrun.sh", variables)


            
