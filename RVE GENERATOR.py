# Author: Mojtaba Abdolkhani
# Please contact me before using the codes: mojtababdolkhani@gmail.com
#
# This script was designed to model cortical bone tissues which contain matrix, fibers, and porosity.
# we modeled the crack initiation & propagation using the Phase-field model, However, for the debonding 
# between fibers and matrix we used the Cohesive zone model
############################################################################
##             Simple Tension Test: RVE GENERATOR                         ##
############################################################################
from abaqus import *
from abaqusConstants import *
from caeModules import *
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
from random import *
from math import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import numpy as np
import sys
import inspect
import os        
import shutil    
#--------------------------------------------Modifying-------------------------------------------------#
#!!! Change your working directory based on where your Umat subroutine is located.
#
#!!! IMPORTANT NOTE: BEFORE START MODIFYING THIS CODE IN ABAQUS CAE, copy paste this line below first in cae command: 
#session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
#-----------------------------------------Input parameters---------------------------------------------#
dis=numpy.zeros(1000)
rad1=0.1                # Radius of osteons 
rad0=0.03               # Radius of haversian canals
Max_iterations=1        # Set number of iterations 
max_incl = 8            # Set number of osteons required 
h=1.0                   # Height of the RVE 
w=1.0                   # Width of the RVE 
load1=0.03              # Load-Displacement on the Top edge of RVE
Hmg_th=0.05             # Homogenized area
preC=0.25               # Pre-crack  
#-----------------------------------------------------------------------------------------------------
for q in range (1,(Max_iterations+1)):
    # LET'S CREATE MODEL
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d' %(q))
    
    # LET'S CREATE PART
    mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', sheetSize=(h+w))
    mdb.models['Model-%d' %(q)].sketches['__profile__'].rectangle(point1=(-w/2, h/2), 
        point2=(w/2, -h/2))
    mdb.models['Model-%d' %(q)].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-%d' %(q)].parts['Part-1'].BaseShell(sketch=
        mdb.models['Model-%d' %(q)].sketches['__profile__'])

    #   Creating pre-crack
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=1.8, name='__profile__', 
        sheetSize=100*(h+w), transform=
        mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(-w/2, -h/2, 0.0)))
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
    
    s1 = mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', 
        sheetSize=(h+w), gridSpacing=0.07)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.Line(point1=(-w/2, 0.0005), point2=(((-w/2.0)+preC), 0.0))
    s1.Line(point1=(((-w/2.0)+preC), 0.0), point2=(-w/2.0, -0.0005))
    s1.Line(point1=(-w/2, -0.0005), point2=(-w/2, 0.0005))
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    p.Cut(sketch=s1)  

    num_incl = 0
    x_coordinate = []
    y_coordinate = []

    while (num_incl < max_incl):
        random_x=random.uniform(((-w/2.0)+rad1+(w/999999.0)), ((w/2.0)-rad1-(w/999999.0)))  #generate random x_coordinate within RVE
        random_y=random.uniform(((-h/2.0)+rad1+(w/999999.0)), ((h/2.0)-rad1-(w/999999.0)))  #generate random y_coordinate within RVE

        isPointIntersecting = False
        
        is_between = False
        is_between = ((-w/2.0) <= random_x <= ((-w/2.0)+preC+rad1+(w/999999.0))) and ((-rad1-0.0005-(h/999999.0))<= random_y <= (+rad1+0.0005+(h/999999.0)))
        # To check if new osteon intersects with any existing osteons
        for j in range (0,len(x_coordinate)):
    
    
            dis[j]=sqrt((random_x-x_coordinate[j])**2+(random_y-y_coordinate[j])**2)
            

                
            if dis[j] < (2.2*rad1):  

                isPointIntersecting = True
                break

        if (isPointIntersecting == False) and (is_between == False):
            x_coordinate.append(random_x)
            y_coordinate.append(random_y)
            num_incl = num_incl + 1  # count no of osteon       

    #   Creating Osteons
    for i in range(max_incl):    

        mdb.models['Model-%d' %(q)].sketches['__profile__'].CircleByCenterPerimeter(center=(
            x_coordinate[i], y_coordinate[i]), point1=((x_coordinate[i]-rad1), y_coordinate[i]))

        mdb.models['Model-%d' %(q)].parts['Part-1'].PartitionFaceBySketch(faces=
            mdb.models['Model-%d' %(q)].parts['Part-1'].faces.findAt(((((w/2)-(w/200.0)), 
            ((h/2)-(h/200.0)), 0.0), (0.0, 0.0, 1.0)), ), sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
    #    del mdb.models['Model-%d' %(q)].sketches['__profile__']
    #   Creating haversian canals
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=1.8, name='__profile__', 
        sheetSize=(h+w), transform=
        mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
    for i in range(max_incl):    

        mdb.models['Model-%d' %(q)].sketches['__profile__'].CircleByCenterPerimeter(center=(
            x_coordinate[i], y_coordinate[i]), point1=((x_coordinate[i]-rad0), y_coordinate[i]))

        mdb.models['Model-%d' %(q)].parts['Part-1'].Cut(sketch=
			mdb.models['Model-%d' %(q)].sketches['__profile__'])
                    
    del mdb.models['Model-%d' %(q)].sketches['__profile__']
    # 
    # LET'S CREATE MATERIAL-1 (MATRIX)
    mdb.models['Model-%d' %(q)].Material(name='Matrix')
    mdb.models['Model-%d' %(q)].materials['Matrix'].UserMaterial(mechanicalConstants=(19700.0, 0.25, 0.008, 0.52, 30.0, 41.0))
    mdb.models['Model-%d' %(q)].materials['Matrix'].Depvar(n=1)
    mdb.models['Model-%d' %(q)].materials['Matrix'].Conductivity(table=((1.0, ), ))    
    # LET'S CREATE MATERIAL-2 (ELASTIC Osteon)
    mdb.models['Model-%d' %(q)].Material(name='Elastic')
    mdb.models['Model-%d' %(q)].materials['Elastic'].UserMaterial(mechanicalConstants=(16600.0, 0.25, 0.008, 0.62, 30.0, 41.0))
    mdb.models['Model-%d' %(q)].materials['Elastic'].Depvar(n=1)
    mdb.models['Model-%d' %(q)].materials['Elastic'].Conductivity(table=((1.0, ), ))

    # LET'S CREATE SECTIONS    
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='Matrix', name='Matrix', 
        thickness=None)
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='Elastic', name='osteon', 
        thickness=None)
        
    # LET'S ASSIGN SECTIONS
        
    mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-%d' %(q)].parts['Part-1'].faces.findAt(((((w/2.0)-((w/9999999.0)/2.0)), 
        0.0, 0.0), (0.0, 0.0, 1.0)), )), sectionName='Matrix', 
        thicknessAssignment=FROM_SECTION)

    for i in range (num_incl):    
        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.2, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            faces=mdb.models['Model-%d' %(q)].parts['Part-1'].faces.findAt((((x_coordinate[i]-((rad0+rad1)/2.0)), 
            (y_coordinate[i]), 0.0), (0.0, 0.0, 1.0)), )), sectionName='osteon', 
            thicknessAssignment=FROM_SECTION)    
    # Creating two parts from one part
    p1 = mdb.models['Model-%d' %(q)].parts['Part-1']
    p = mdb.models['Model-%d' %(q)].Part(name='Part-1-Copy', 
        objectToCopy=mdb.models['Model-%d' %(q)].parts['Part-1'])
    p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
    f = p.faces
    p.RemoveFaces(faceList=(f.findAt(coordinates=(((w/2.0)), 0.0, 0.0)), 
        f.findAt(coordinates=(((w/2.0)-((w/9999999.0)/2.0)), 0.0, 0.0))), deleteCells=False)

    for i in range (num_incl):
     p = mdb.models['Model-%d' %(q)].parts['Part-1']
     f1 = p.faces
     p.RemoveFaces(faceList=(f1.findAt(coordinates=((x_coordinate[i]-((rad0+rad1)/2.0)), y_coordinate[i], 0.0)), ), 
         deleteCells=False)
    #Creating surfaces for the interface
    f45 = []
    f46 = []
    for i in range (num_incl):
        
        p = mdb.models['Model-%d' %(q)].parts['Part-1']
        s = p.edges
        s55 = s.findAt((((x_coordinate[i]+rad1), y_coordinate[i], 0.0), ))
        f45.append(s55)
        p.Surface(side1Edges=f45, name='Surf-Matrix')
    
    for i in range (num_incl):    
        p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
        s = p.edges
        s56 = s.findAt((((x_coordinate[i]+rad1), y_coordinate[i], 0.0), ))
        f46.append(s56)
    p.Surface(side1Edges=f46, name='Surf-Osteons')
    # Creating Assembly
    a = mdb.models['Model-%d' %(q)].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    a.Instance(name='Part-1-1', part=p, dependent=ON)
    p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
    a.Instance(name='Part-1-Copy-1', part=p, dependent=ON)

    # Creating All set 
    f11 = []    
    for i in range (num_incl):
        a = mdb.models['Model-%d' %(q)].rootAssembly
        f1 = a.instances['Part-1-Copy-1'].faces.findAt((((x_coordinate[i]-((rad0+rad1)/2.0)), 
            (y_coordinate[i]), 0.0), (0.0, 0.0, 1.0)), )
        f11.append(f1)
        
    a = mdb.models['Model-%d' %(q)].rootAssembly
    f12 = a.instances['Part-1-1'].faces
    f2 = f12.findAt(((((w/2.0)), 0.0, 0.0), ))
    f11.append(f2)
    
    a = mdb.models['Model-%d' %(q)].rootAssembly
    f13 = a.instances['Part-1-1'].faces
    f3 = f13.findAt(((((w/2.0)-((w/9999999.0)/2.0)), 0.0, 0.0), ))
    f11.append(f3)
    
    a.Set(faces=f11, name='All')
    
    # LET'S CREATE STEP
    mdb.models['Model-%d' %(q)].CoupledTempDisplacementStep(name='Step-1', 
        previous='Initial', response=STEADY_STATE, maxNumInc=1000000, 
        initialInc=0.01, minInc=0.01, maxInc=0.01, deltmx=None, cetol=None, 
        creepIntegration=None, amplitude=RAMP, extrapolation=NONE, matrixStorage=SYMMETRIC, 
        solutionTechnique=SEPARATED)    
    #Controls options
    mdb.models['Model-%d' %(q)].steps['Step-1'].control.setValues(allowPropagation=OFF, 
        resetDefaultValues=OFF, timeIncrementation=(5000.0, 5000.0, 5000.0, 
        5000.0, 5000.0, 5000.0, 12.0, 5.0, 6.0, 3.0, 50.0))
               
    # LET'S CREATE BOUNDARY CONDITIONS
    a = mdb.models['Model-%d' %(q)].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.findAt(((-w/4, -h/2, 0.0), ))
    region = a.Set(edges=edges1, name='Bot')
    mdb.models['Model-%d' %(q)].DisplacementBC(name='BC-bot_fixY', 
        createStepName='Initial', region=region, u1=SET, u2=SET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-%d' %(q)].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.findAt(((w/4, h/2, 0.0), ))
    region = a.Set(edges=edges1, name='Top1')
    mdb.models['Model-%d' %(q)].DisplacementBC(name='BC-top_fixX', 
        createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-%d' %(q)].rootAssembly
    e1 = a.instances['Part-1-1'].edges
    edges1 = e1.findAt(((w/4, h/2, 0.0), ))
    region = a.Set(edges=edges1, name='Top')
    mdb.models['Model-%d' %(q)].DisplacementBC(name='BC-load', createStepName='Step-1', 
        region=region, u1=UNSET, u2=load1, ur3=UNSET, amplitude=UNSET, 
        fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
        
    # PREDEFINED FIELDS
    a = mdb.models['Model-%d' %(q)].rootAssembly
    region = a.sets['All']
    mdb.models['Model-%d' %(q)].Temperature(name='Predefined Field-1', 
        createStepName='Initial', region=region, distributionType=UNIFORM, 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, 
        ))
    # CZM
    mdb.models['Model-%d' %(q)].ContactProperty('IntProp-1')
    mdb.models['Model-%d' %(q)].interactionProperties['IntProp-1'].TangentialBehavior(
     formulation=FRICTIONLESS)
    mdb.models['Model-%d' %(q)].interactionProperties['IntProp-1'].NormalBehavior(
     pressureOverclosure=HARD, allowSeparation=ON, 
     constraintEnforcementMethod=DEFAULT)
    mdb.models['Model-%d' %(q)].interactionProperties['IntProp-1'].CohesiveBehavior(
        defaultPenalties=OFF, table=((1000000.0, 1000000.0, 1000000.0), ))
    mdb.models['Model-%d' %(q)].interactionProperties['IntProp-1'].Damage(
        criterion=QUAD_TRACTION, initTable=((7.31, 7.31, 7.31), ), 
        useEvolution=ON, evolutionType=ENERGY, evolTable=((0.1629, ), ), 
        useStabilization=ON, viscosityCoef=0.005)
    a = mdb.models['Model-%d' %(q)].rootAssembly
    region1=a.instances['Part-1-Copy-1'].surfaces['Surf-Osteons']
    a = mdb.models['Model-%d' %(q)].rootAssembly
    region2=a.instances['Part-1-1'].surfaces['Surf-Matrix']
    mdb.models['Model-%d' %(q)].SurfaceToSurfaceContactStd(name='Int-1', 
        createStepName='Initial', master=region1, slave=region2, sliding=SMALL, 
        thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE, 
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
    #------------------mesh madules-----------    
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    p.seedPart(size=(w+h)/210.0, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
    p.seedPart(size=(w+h)/210.0, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    elemType1 = mesh.ElemType(elemCode=CPE4T, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=CPE3T, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    f = p.faces
    pickedRegions = f.findAt(((((w/2.0)-((w/9999999.0)/2.0)), 0.0, 0.0), ),)
    pickedRegions =(pickedRegions, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    elemType1 = mesh.ElemType(elemCode=CPE4T, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=CPE3T, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    f89 = []
    for i in range (num_incl):    
        p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
        f = p.faces
        faces = f.findAt((((x_coordinate[i]-((rad0+rad1)/2.0)), 
            (y_coordinate[i]), 0.0), (0.0, 0.0, 1.0)), )
        f89.append(faces)        
    pickedRegions =(f89, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p = mdb.models['Model-%d' %(q)].parts['Part-1-Copy']
    p.generateMesh()
    p = mdb.models['Model-%d' %(q)].parts['Part-1']
    p.generateMesh()     
    
    # LET'S CREATE HISTORY OUTPUT REQUESTS
    mdb.models['Model-%d' %(q)].FieldOutputRequest(name='F-Output-1', 
        createStepName='Step-1', variables=('U', 'NT', 'SDV'))
    regionDef=mdb.models['Model-%d' %(q)].rootAssembly.sets['Top']
    mdb.models['Model-%d' %(q)].HistoryOutputRequest(name='H-Output-1', 
        createStepName='Step-1', variables=('U2', 'RF2'), region=regionDef, 
        sectionPoints=DEFAULT, rebar=EXCLUDE)
     
    #LET'S CREATE JOBS 
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-%d' %(q), modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='Job-%d' %(q) , nodalOutputPrecision=SINGLE, 
        numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='C:\\temp\\test\\2D_Brittle.for', 
        waitHours=0, waitMinutes=0)
    mdb.jobs['Job-%d' %(q)].writeInput()
    mdb.jobs['Job-%d' %(q) ].submit(consistencyChecking=OFF)    
    mdb.jobs['Job-%d' %(q) ].waitForCompletion()
    ###################################################################################
    # session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=100., 
    # height=100.)
    # session.viewports['Viewport: 1'].makeCurrent()
    # session.viewports['Viewport: 1'].maximize()
    # my_inp = open('Job-1.inp', 'r')
    # string_list_inp = my_inp.readlines()
    # my_inp.close()
    # count = 0
    # NsetList_old=[]
    # for line in string_list_inp:
        # count += 1
        # tempII = '{}'.format(line.strip())
        # if (tempII == '*Nset, nset=Top, instance=Part-1-1'):
         # j=count
         # NsetList_old =string_list_inp[j].strip()
         # while (string_list_inp[j+1] != '*Elset, elset=Top, instance=Part-1-1\n'):
             # j += 1
             # temp2=string_list_inp[j].strip()
             # NsetList_old = NsetList_old+","+temp2
         # NsetList_old=NsetList_old.split(",")   
         # for ii in range(0,len(NsetList_old)):
          # NsetList_old[ii] = int(NsetList_old[ii])
         # break  
    # NsetList = NsetList_old
    # ########################################################################################
    # TopNode_iterations = len(NsetList)
    # TopLeftNode = NsetList[1]
    # ######################################################################
    # def OpenODBandPLOT(odbPath, curveName, sessionPlot):
        # o1 = session.openOdb(name = odbPath)
        # session.viewports['Viewport: 1'].setValues(displayedObject=o1)
        # odb = session.odbs[odbPath]

        # for qqqq in range(0,TopNode_iterations):
            # session.XYDataFromHistory(name='XY-%d' %(qqqq+1), odb=odb, 
                # outputVariableName='Reaction force: RF2 PI: PART-1-1 Node %d in NSET TOP' %(NsetList[qqqq]), )

        # session.XYDataFromHistory(name='XY-%d' %(TopNode_iterations+1), odb=odb, 
            # outputVariableName='Spatial displacement: U2 PI: PART-1-1 Node %d in NSET TOP' %(TopLeftNode), )
            
        # this = sys.modules[__name__]
        # vList = []
        # for qq in range(0,TopNode_iterations+1):
            # setattr(this, 'xy%d' %(qq+1), session.xyDataObjects['XY-%d' %(qq+1)])
        # for qq in range(0,TopNode_iterations):
            # vList.append(session.xyDataObjects['XY-%d' %(qq+1)])
    # ###########  
        # vList1=session.xyDataObjects['XY-%d' %(TopNode_iterations+1)]
        # l5=combine(vList1, sum(vList))
        # xl5=[]
        # yl5=[]
        # for s11 in l5:
         # xl5.append(s11[0])
        # for s12 in l5:
         # yl5.append(s12[1])     
        # output = np.trapz(yl5,xl5)
        # f3 = open("Job-1.dat", "a")
        # f3.write('\noutput:\n')
        # f3.write(str(output))
        # f3.close()
        # setattr(this, 'xy%d' %(TopNode_iterations+2), l5)
        # l5.setValues(sourceDescription=curveName)
        # tmpName = l5.name
        # session.xyDataObjects.changeKey(tmpName, curveName)
    # ###########
        # for qqq in range(0,TopNode_iterations+1):
         # del session.xyDataObjects['XY-%d' %(qqq+1)]
        

    # cwd = os.getcwd()

    # for file in os.listdir(cwd):
        # if file.endswith(".odb"):
            # SessionPlotName = 'Session-'+str(1)
            # ODBname = os.path.join(cwd,file)
            # OpenODBandPLOT(ODBname, file, SessionPlotName)    
         
