#> \author Elias Ghadam Soltani
#> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. If you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. If you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>





#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,math,cmath,csv,time,sys,os,pdb
from opencmiss.iron import iron

# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
#iron.OutputSetOn("Testing")
numberOfRandomSeeds = iron.RandomSeedsSizeGet()
randomSeeds = [0]*numberOfRandomSeeds
randomSeeds[0] = 100
iron.RandomSeedsSet(randomSeeds)
# Get the computational nodes info
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber    = iron.ComputationalNodeNumberGet()

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================
numberOfDimensions     = 1  #(One-dimensional)

derivIdx = 1

ProgressDiagnostics = False   # Set to diagnostics


#================================================================================================================================
#  Start Program
#================================================================================================================================

# Set program variables
CoordinateSystemUserNumber                = 1
RegionUserNumber                          = 2
BasisUserNumberSpace                      = 3
MeshUserNumber                            = 4
DecompositionUserNumber                   = 5
GeometricFieldUserNumber                  = 6
EquationsSetUserNumberNavierStokes        = 7
EquationsSetFieldUserNumberNavierStokes   = 8
DependentFieldUserNumber                  = 9
materialsFieldUserNumber                  = 10
sourceFieldUserNumber                     = 11
IndependentFieldUserNumber                = 12
ProblemUserNumber                         = 13

#================================================================================================================================
#  Mesh Reading
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> Reading geometry from files... << == ")

# Read the node file
with open('input/nodes.csv','r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    rownum=0
    for row in reader:
        if (rownum==0):
            header = row
        else:
            if (rownum==1):
                totalNumberOfNodes = int(row[4])
                xValues = numpy.zeros((totalNumberOfNodes+1,1),dtype = numpy.float)
                yValues = numpy.zeros((totalNumberOfNodes+1,1),dtype = numpy.float)
                zValues = numpy.zeros((totalNumberOfNodes+1,1),dtype = numpy.float)
                A0  = numpy.zeros((totalNumberOfNodes+1,1),dtype = numpy.float)
            xValues[rownum] = float(row[0])*1000
            yValues[rownum] = float(row[1])*1000
            zValues[rownum] = float(row[2])*1000
            A0[rownum]  = float(row[3])*1000000
        rownum+=1

 # Read the element file
with open('input/elements.csv','r') as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
    rownum=0
    for row in reader:
        if (rownum==0):
            header = row
        else:
            if (rownum==1):
                totalNumberOfElements=int(row[3])
                elementNodes = (totalNumberOfElements+1)*[3*[0]]
            elementNodes[rownum]=[int(row[0]),int(row[1]),int(row[2])]
        rownum+=1

#================================================================================================================================
#  Initial Data & Default Values
#================================================================================================================================

#-------------------=========
Alpha = 1.57e-7*1000000        # mm2/s Diffusivity
U    = 0.07                                # m/s flow velocity

# Set the time parameters
timeIncrement   = 0.1
startTime       = 0.0
stopTime  = 100.1

# Set the output parameters
DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY = 100

# Set the solver parameters
#relativeTolerance   = 1.0E-05  # default: 1.0E-05
#absoluteTolerance   = 1.0E-08  # default: 1.0E-10
#DIVERGENCE_TOLERANCE = 1.0e+10  # default: 1.0e+05
MAXIMUM_ITERATIONS   = 1000   # default: 100000
#RESTART_VALUE        = 3000     # default: 30


# Navier-Stokes solver
EquationsSetSubtype = iron.EquationsSetSubtypes.ADVECTION_DIFFUSION
ProblemSubtype      = iron.ProblemSubtypes.LINEAR_SOURCE_ADVECTION_DIFFUSION

#================================================================================================================================
#  Coordinate System
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> COORDINATE SYSTEM << == ")

# Start the creation of RC coordinate system
CoordinateSystem = iron.CoordinateSystem()
CoordinateSystem.CreateStart(CoordinateSystemUserNumber)
CoordinateSystem.DimensionSet(3)
CoordinateSystem.CreateFinish()

#================================================================================================================================
#  Region
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> REGION << == ")

# Start the creation of SPACE region
region = iron.Region()
region.CreateStart(RegionUserNumber,iron.WorldRegion)
region.label = "ArterialSystem"
region.coordinateSystem = CoordinateSystem
region.CreateFinish()

#================================================================================================================================
#  Bases
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> BASIS << == ")

# Start the creation of SPACE bases
basisXiGaussSpace = 3
BasisSpace = iron.Basis()
BasisSpace.CreateStart(BasisUserNumberSpace)
BasisSpace.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
BasisSpace.numberOfXi = numberOfDimensions
BasisSpace.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
BasisSpace.quadratureNumberOfGaussXi = [basisXiGaussSpace]
BasisSpace.CreateFinish()

#================================================================================================================================
#  Nodes
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> NODES << == ")

# Start the creation of mesh nodes
Nodes = iron.Nodes()
Nodes.CreateStart(region,totalNumberOfNodes)
Nodes.CreateFinish()

#================================================================================================================================
#  Mesh
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MESH << == ")

# Start the creation of SPACE mesh
Mesh = iron.Mesh()
Mesh.CreateStart(MeshUserNumber,region,numberOfDimensions)
Mesh.NumberOfElementsSet(totalNumberOfElements)
Mesh.NumberOfComponentsSet(1)
MeshElementsSpace = iron.MeshElements()

# Specify the SPACE mesh component
meshComponentNumberSpace = 1
MeshElementsSpace.CreateStart(Mesh,meshComponentNumberSpace,BasisSpace)

for elemIdx in range(1,totalNumberOfElements+1):
  MeshElementsSpace.NodesSet(elemIdx,elementNodes[elemIdx])
MeshElementsSpace.CreateFinish()

# Finish the creation of the meshGeometricFieldUserNumber
Mesh.CreateFinish()

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MESH DECOMPOSITION << == ")

# Start the creation of SPACE mesh decomposition
Decomposition = iron.Decomposition()
Decomposition.CreateStart(DecompositionUserNumber,Mesh)
Decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
Decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
Decomposition.CreateFinish()

#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> GEOMETRIC FIELD << == ")

# Start the creation of SPACE geometric field
GeometricField = iron.Field()
GeometricField.CreateStart(GeometricFieldUserNumber,region)
GeometricField.meshDecomposition = Decomposition
GeometricField.NumberOfVariablesSet(1)
GeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'Coordinates')
GeometricField.TypeSet = iron.FieldTypes.GEOMETRIC
GeometricField.ScalingTypeSet = iron.FieldScalingTypes.NONE

for componentNumber in range(1,CoordinateSystem.dimension+1):
    GeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentNumber,1)

GeometricField.CreateFinish()

# Set the geometric field values for version 1
versionIdx = 1
for nodeIdx in range(1,totalNumberOfNodes+1):
    nodeDomain = Decomposition.NodeDomainGet(nodeIdx,1)
    if (nodeDomain == computationalNodeNumber):
        GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx,1,xValues[nodeIdx][0])
        GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx,2,yValues[nodeIdx][0])
        GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx,3,zValues[nodeIdx][0])         

# Finish the parameter update
GeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
GeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#================================================================================================================================
#  Equations Sets
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> EQUATIONS SET << == ")

# Create the equations set for NAVIER-STOKES
EquationsSetNavierStokes = iron.EquationsSet()
EquationsSetFieldNavierStokes = iron.Field()
# Set the equations set to be a dynamic nonlinear problem
NavierStokesEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
				      iron.EquationsSetTypes.ADVECTION_EQUATION,
				      EquationsSetSubtype]
EquationsSetNavierStokes.CreateStart(EquationsSetUserNumberNavierStokes,region,GeometricField,
    NavierStokesEquationsSetSpecification,EquationsSetFieldUserNumberNavierStokes,EquationsSetFieldNavierStokes)
EquationsSetNavierStokes.CreateFinish()

#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> DEPENDENT FIELD << == ")

# Create the equations set dependent field variables
DependentFieldNavierStokes = iron.Field()
EquationsSetNavierStokes.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)  
DependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U,'Temperature')
DependentFieldNavierStokes.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
DependentFieldNavierStokes.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
EquationsSetNavierStokes.DependentCreateFinish()

# Initialise dependent field
DependentFieldNavierStokes.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,37.0)

DependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
DependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#================================================================================================================================
#  Materials Field
#================================================================================================================================
 
if (ProgressDiagnostics):
    print( " == >> MATERIALS FIELD << == ")

#!  !Create the equations set material field variables
materialsField = iron.Field()
EquationsSetNavierStokes.MaterialsCreateStart(materialsFieldUserNumber,materialsField)
materialsField.VariableLabelSet(iron.FieldVariableTypes.U,'Materials')
materialsField.ComponentLabelSet(iron.FieldVariableTypes.U,1,'Diffusivity')
materialsField.ComponentLabelSet(iron.FieldVariableTypes.U,2,'Source T coeff.')
EquationsSetNavierStokes.MaterialsCreateFinish()


# diffusivity=1.57e-7+U*beta*le/2 #U*beta*le/2=0.000416667 almost 3000 times of the real diffusivity Pe=Ule/2a=0.2*0.05/12/2/0.0004=1
diffusivity=Alpha
materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,diffusivity)
materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    math.pi*4*diffusivity) # mm2/s. b-cT. b=pi*Nu*alpha/A * Tw and c = pi*Nu*alpha/A. We still need to divide by cross-section area.

materialsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
materialsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#================================================================================================================================
#  Source Field
#================================================================================================================================
  

if (ProgressDiagnostics):
    print( " == >> SOURCE FIELD << == ")

sourceField = iron.Field()
EquationsSetNavierStokes.SourceCreateStart(sourceFieldUserNumber,sourceField)
EquationsSetNavierStokes.SourceCreateFinish()  

sourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
sourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


#================================================================================================================================
# Independent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print (" == >> INDEPENDENT FIELD << == ")

# Create the equations set independent field variables
IndependentFieldNavierStokes = iron.Field()
#IndependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U,'flow velocity')
# Set the mesh component to be used by the field components.
#IndependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)

# NAVIER-STOKES
EquationsSetNavierStokes.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
EquationsSetNavierStokes.IndependentCreateFinish()

# Set the velocity
for nodeIdx in range(1,totalNumberOfNodes+1):
  nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
  if (nodeDomain == computationalNodeNumber):
    for timeStep in range(1,12):
        comp2=timeStep*2
        comp1=comp2-1
        IndependentFieldNavierStokes.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
        versionIdx,derivIdx,nodeIdx,comp1,U*A0[nodeIdx][0])
        IndependentFieldNavierStokes.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
        versionIdx,derivIdx,nodeIdx,comp2,A0[nodeIdx][0])     

# Finish the parameter update
IndependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
IndependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#================================================================================================================================
#  Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> EQUATIONS << == ")

# Create equations
equations = iron.Equations()
EquationsSetNavierStokes.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
EquationsSetNavierStokes.EquationsCreateFinish()

# I want to solve this type of equation, dT/dt+udT/dx-alpha d2T/dx2-(b-cT)=0. 

#================================================================================================================================
#  Problems
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> PROBLEM << == ")

# Start the creation of a problem.
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                        iron.ProblemTypes.ADVECTION_DIFFUSION_EQUATION,
                        ProblemSubtype]
problem.CreateStart(ProblemUserNumber,problemSpecification)
problem.CreateFinish()

#================================================================================================================================
#  Control Loops
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> PROBLEM CONTROL LOOP << == ")

# Create control loops
problem.ControlLoopCreateStart()
TimeLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],TimeLoop)
TimeLoop.LabelSet('Time Loop')
TimeLoop.TimesSet(startTime,stopTime,timeIncrement)
TimeLoop.TimeOutputSet(DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY)
problem.ControlLoopCreateFinish()

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> SOLVERS << == ")

# Create problem solver
solver = iron.Solver()
LinearSolver = iron.Solver()

problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
#solver.outputType = iron.SolverOutputTypes.SOLVER
solver.DynamicLinearSolverGet(LinearSolver)
#solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.linearIterativeAbsoluteTolerance = 1.0E-12
#solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> SOLVER EQUATIONS << == ")

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(EquationsSetNavierStokes)
problem.SolverEquationsCreateFinish()

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> BOUNDARY CONDITIONS << == ")

boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

nodes = iron.Nodes()
region.NodesGet(nodes)

# for nodeNumber in boundary:
nodeDomain = Decomposition.NodeDomainGet(1,1)
if nodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(DependentFieldNavierStokes,iron.FieldVariableTypes.U,1,1,1,1,
    iron.BoundaryConditionsTypes.FIXED,[37.0])


DependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
DependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
#!  !Finish the creation of the equations set boundary conditions
#!  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)	

solverEquations.BoundaryConditionsCreateFinish()

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

# Solve the problem
print( "Solving problem...")

# Solve the problem
start = time.time()
problem.Solve()
end = time.time()
elapsed = end - start
print( "Total Number of Elements = %d " %totalNumberOfElements)
print( "Calculation Time = %3.4f" %elapsed)
print( "Problem solved!")
print( "#")

# Export results
#baseName = "laplace"
#dataFormat = "PLAIN_TEXT"
#fml = iron.FieldMLIO()
#fml.OutputCreate(Mesh, "", baseName, dataFormat)
#fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, GeometricField,
#    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#fml.OutputAddFieldNoType(baseName+".phi", dataFormat, DependentFieldNavierStokes,
#    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#fml.OutputWrite("LaplaceExample.xml")
#fml.Finalise()

iron.Finalise()

#================================================================================================================================
#  Finish Program
#================================================================================================================================













