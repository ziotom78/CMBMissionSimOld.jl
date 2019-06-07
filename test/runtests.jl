import CMBMissionSim
const CMBSim = CMBMissionSim
using Test

@testset "Quaternions" begin
    include("test_quaternions.jl")
end

@testset "Pointing generation" begin
    include("test_genpointings.jl")
end

@testset "Dipole analysis" begin
    include("test_dipole.jl")
end

@testset "Beam functions" begin
    include("test_beams.jl")
end
