import salome
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import os
import sys
import numpy as np
from pathlib import Path

def circle_points(center, normal, radius, n_points=20, include_endpoint=False):
    C = np.asarray(center, dtype=float)
    n = np.asarray(normal, dtype=float)
    if radius <= 0:
        raise ValueError("radius must be positive")
    if np.linalg.norm(n) == 0:
        raise ValueError("normal must be non-zero")
    # unit normal
    n = n / np.linalg.norm(n)
    # pick a vector not parallel to n to build an orthonormal basis
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    # u and v span the circle plane
    u = np.cross(n, a)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    # angles
    thetas = np.linspace(0.0, 2*np.pi, num=n_points, endpoint=include_endpoint)
    # points: C + r*(cosθ u + sinθ v)
    pts = C + radius*np.cos(thetas)[:,None]*u + radius*np.sin(thetas)[:,None]*v
    # find index of maximum x coordinate
    idx = np.argmax(pts[:,0])
    # rotate list so that the max-x point is first
    pts = np.roll(pts, -idx, axis=0)
    return pts



# helixradius = 0.001
# pitch = 0.12
# numOfPointsPerPitch = 15 # longer pitch requires more points, shorter pitch requires fewer points
# height = 0.12
# pipe_radius = 0.0067204
# helixDirectionFactor = 1
# surfaceMeshSize = 0.002
# gradualTurnNum = 1
# output_directory = "./aaa/"
# method = "robust"


import loadJSONPara
############ read from json for input parameters #############
try:
    json_file_path = sys.argv[1]
    print("load json path = {} ".format(json_file_path))
except:
    raise("no json path")
# # parameters (m)
numOfPointsPerPitch = loadJSONPara.readwithdefault(json_file_path,"numOfPointsPerPitch",30) 
helixradius = loadJSONPara.readwithdefault(json_file_path,"helixSpineAmplitude",0.005)  
height = loadJSONPara.readwithdefault(json_file_path,"pipeLength",7) 
pitch =loadJSONPara.readwithdefault(json_file_path,"helixSpinePeriod",0.159) 

helixDirectionFactor = loadJSONPara.readwithdefault(json_file_path,"helixDirectionFactor",1) 
pipe_radius = loadJSONPara.readwithdefault(json_file_path,"pipeRadius",0.02815) 
surfaceMeshSize = loadJSONPara.readwithdefault(json_file_path,"surfaceMeshSize",0.006) 
gradualTurnNum = loadJSONPara.readwithdefault(json_file_path,"gradualTurnNum",1) 
output_directory = loadJSONPara.readwithdefault(json_file_path,"output_directory",".") 
method = loadJSONPara.readwithdefault(json_file_path,"method","robust") 

############ Salome starts #############
geompy = geomBuilder.New()
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

t_step = np.linspace(0, height/pitch * 2*np.pi, int(numOfPointsPerPitch * height/pitch))
helixlocation = []
for t in t_step:
    if gradualTurnNum >= 1:
        r =  helixradius * min(t/(gradualTurnNum * 2*np.pi),(height/pitch * 2*np.pi - t)/(gradualTurnNum * 2*np.pi), 1)  
    else:
        r = helixradius
    x = r*math.cos(t * helixDirectionFactor)
    y = r*math.sin(t * helixDirectionFactor)
    z = t * pitch /(2*np.pi)
    helixlocation.append((x,y,z))

firstDerivativeList = []
for t in t_step:
    if gradualTurnNum >= 1:
        r =  helixradius * min(t/(gradualTurnNum * 2*np.pi),(height/pitch * 2*np.pi - t)/(gradualTurnNum * 2*np.pi), 1)  
    else:
        r = helixradius
    x = -r * np.sin(t)
    y = r * np.cos(t)
    z = pitch/(2*np.pi)
    firstDerivativeList.append(np.array([x,y,z]))

vertexList = []
for i,helixloc in enumerate(helixlocation):
    Vertex_tmp = geompy.MakeVertex(helixloc[0], helixloc[1], helixloc[2])
    vertexList.append(Vertex_tmp)

HelixCurve = geompy.MakePolyline(vertexList, False)
geompy.addToStudy( HelixCurve, 'HelixCurve' )


if method == "smooth": # maintain first order continuity  only for small theta(require the normal trick). For large theta many problems
    thetaConstrainstDegree = 15
    theta = math.atan(helixradius / (pitch/4)) * 180/np.pi
    if theta < thetaConstrainstDegree:
        print("using modified normal")
    else:
        print("using real normal")

    circleList =  []
    faceList =  []
    for i,vertex in enumerate(vertexList):
        if theta < thetaConstrainstDegree:
            O_axis_tmp = geompy.MakeVectorDXDYDZ(0, 0, 1)
        else:
            O_axis_tmp = geompy.MakeVectorDXDYDZ(firstDerivativeList[i][0], firstDerivativeList[i][1], firstDerivativeList[i][2])
        tmp_circle = geompy.MakeCircle(vertex, O_axis_tmp, pipe_radius)
        circleList.append(tmp_circle)
        geompy.addToStudy(tmp_circle, "circle_{}".format(i))
        tmp_face = geompy.MakeFaceWires([tmp_circle], 1)
        faceList.append(tmp_face)
    HelixPipeWall = geompy.MakePipeWithDifferentSections(circleList, vertexList, HelixCurve,False,False)
    geompy.addToStudy( HelixPipeWall, 'HelixPipeWall' )
    
    # make the first circle face
    first_circle_Face = faceList[0]
    geompy.addToStudy(first_circle_Face, 'first_circle_Face' )

    O_axis_lastCircle = geompy.MakeVectorDXDYDZ(firstDerivativeList[-1][0], firstDerivativeList[-1][1], firstDerivativeList[-1][2])
    last_circle = geompy.MakeCircle(vertexList[-1], O_axis_lastCircle, pipe_radius)
    last_circle_Face = geompy.MakeFaceWires([last_circle], 1)
    geompy.addToStudy( last_circle_Face, 'last_circle_Face' )

else:   #robust method but does not maintain first order continuity across different segment of the pipe
    circleList =  []
    faceList =  []
    for i,vertex in enumerate(helixlocation):
        normal = firstDerivativeList[i]
        #using manual points instead of build in circle since creating filling requires carefully designed order of helix 
        points = circle_points(vertex, normal, pipe_radius, n_points=20) 
        salome_points = [ geompy.MakeVertex(x,y,z) for (x,y,z) in points ]
        bspline = geompy.MakeInterpol(salome_points,True)    
        #geompy.addToStudy(bspline, "CircleBSpline_{}".format(i))
        circleList.append(bspline)
        tmp_face = geompy.MakeFaceWires([bspline], 1)
        faceList.append(tmp_face)

    fillingList = []
    for i,(circle,circleNext) in enumerate(zip(circleList,circleList[1:])):
        Filling_1 = geompy.MakeFilling([circle, circleNext])
        #geompy.addToStudy(Filling_1, "Filling_{}".format(i))
        fillingList.append(Filling_1)

    HelixPipeWall = geompy.MakeFuseList(fillingList, True, True)
    geompy.addToStudy( HelixPipeWall, 'HelixPipeWall' )

    # make the first circle face
    first_circle_Face = faceList[0]
    geompy.addToStudy(first_circle_Face, 'first_circle_Face' )

    O_axis_lastCircle = geompy.MakeVectorDXDYDZ(firstDerivativeList[-1][0], firstDerivativeList[-1][1], firstDerivativeList[-1][2])
    last_circle = geompy.MakeCircle(vertexList[-1], O_axis_lastCircle, pipe_radius)
    last_circle_Face = geompy.MakeFaceWires([last_circle], 1)
    geompy.addToStudy( last_circle_Face, 'last_circle_Face' )

HelixPipeSurface = geompy.MakeFuseList([HelixPipeWall, first_circle_Face, last_circle_Face], True, True)
geompy.addToStudy( HelixPipeSurface, 'HelixPipeSurface' )
HelixPipeVolume = geompy.MakeSolid([HelixPipeSurface])
geompy.addToStudy( HelixPipeVolume, 'HelixPipeVolume' )

# save the pipe
path = output_directory + '/rawHelixPipeWall.vtk'
try:
  geompy.ExportVTK(HelixPipeWall, path, 0.001)
  pass
except:
  print('ExportPartToSTL() failed. Invalid file name?')


###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

aFilterManager = smesh.CreateFilterManager()
Mesh_1 = smesh.Mesh(HelixPipeSurface,'Mesh_1')
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( surfaceMeshSize * 2 )
NETGEN_2D_Parameters_1.SetMinSize( surfaceMeshSize )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_1.SetUseDelauney( 160 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 3 )
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,first_circle_Face)
aCriteria.append(aCriterion)
aFilter_1 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_1.SetMesh(Mesh_1.GetMesh())
inlet = Mesh_1.GroupOnFilter( SMESH.FACE, 'inlet', aFilter_1 )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,last_circle_Face)
aCriteria.append(aCriterion)
aFilter_2 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_2.SetMesh(Mesh_1.GetMesh())
outlet = Mesh_1.GroupOnFilter( SMESH.FACE, 'outlet', aFilter_2 )

aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,first_circle_Face,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
aCriterion = smesh.GetCriterion(SMESH.FACE,SMESH.FT_BelongToGenSurface,SMESH.FT_Undefined,last_circle_Face,SMESH.FT_LogicalNOT)
aCriteria.append(aCriterion)
aFilter_3 = smesh.GetFilterFromCriteria(aCriteria)
aFilter_3.SetMesh(Mesh_1.GetMesh())
pipewall = Mesh_1.GroupOnFilter( SMESH.FACE, 'pipewall', aFilter_3 )

#### output ######
path = output_directory + '/pipewall.stl'
try:
  Mesh_1.ExportSTL(path, 1, pipewall)
  pass
except:
  print('ExportPartToSTL() failed. Invalid file name?')

path =  output_directory + '/inlet.stl'
try:
  Mesh_1.ExportSTL(path, 1, inlet)
  pass
except:
  print('ExportPartToSTL() failed. Invalid file name?')

path =  output_directory + '/outlet.stl'
try:
  Mesh_1.ExportSTL(path, 1, outlet)
  pass
except:
  print('ExportPartToSTL() failed. Invalid file name?')

## Set names of Mesh objects
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(inlet, 'inlet')
smesh.SetName(outlet, 'outlet')
smesh.SetName(pipewall, 'pipewall')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')

Mesh_2 = smesh.Mesh(HelixPipeVolume,'Mesh_2')
NETGEN_1D_2D_3D = Mesh_2.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
isDone = Mesh_2.Compute()
Mesh_2.CheckCompute()
smesh.SetName(Mesh_2, 'Mesh_2')

#### output ######
pipeVolumeBinary_path = output_directory + '/pipeVolumeBinary.vtk'
try:
  Mesh_2.ExportMESHIO(pipeVolumeBinary_path, 'VTK', Mesh_2)
  pass
except:
  print('ExportPartToSTL() failed. Invalid file name?')


  # old format is required for liggghts
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()



