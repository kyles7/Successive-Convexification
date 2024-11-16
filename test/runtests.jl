# MainModule MUST BE LOADED FIRST, THEN THE TEST MODULES
include("../src/MainModule.jl")
include("unit/test_6dofdynamics.jl")
include("unit/test_dynamics3dof.jl")
include("unit/test_convex_subproblem.jl")
# include("test_scvx.jl")

using .MainModule
using .Test6DOFDynamics
using .Test3DOFDynamics
using .TestConvexSubproblem
using LinearAlgebra
using Test
using YAML

# run the tests 
# runTest6DOFDynamics()
# runTest3DOFDynamics()
runTestConvexSubproblem()
