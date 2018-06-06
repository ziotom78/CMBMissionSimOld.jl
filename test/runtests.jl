using CMBMissionSim
using Base.Test

#######################################################################
# Test case

dirs, ψ = genpointings(0:3600:86400, [0, 0, 1], 0.0,
                       spinsunang=0,
                       borespinang=0,
                       spinrpm=0,
                       precrpm=0,
                       hwprpm=0)

# Check the colatitudes
for idx in 1:size(dirs)[1]
    # Check that we're on the Ecliptic plane (colatitude = 90°)
    @test dirs[idx, 1] ≈ π / 2
end

# Check that the longitude increases as expected
@test dirs[end, 2] - dirs[1, 2] ≈ 2π / DAYS_PER_YEAR

for idx in eachindex(ψ)
    @test ψ[idx] ≈ -π / 2
end

#######################################################################
# Test case

dirs, ψ = genpointings(0:0.1:120, [0, 0, 1], 0.0,
                       spinsunang=0,
                       borespinang=0,
                       spinrpm=1,
                       precrpm=0,
                       hwprpm=1)

@test maximum(ψ) ≈ π / 2
@test minimum(ψ) ≈ -π / 2

#######################################################################
# Test case

dirs = rad2deg.(genpointings(0:1:60, [0, 0, 1], 0.0,
                             borespinang=deg2rad(15),
                             spinsunang=deg2rad(0),
                             spinrpm=1,
                             precrpm=0,
                             yearlyrpm=0.1)[1])

# Colatitudes should depart no more than ±15° from the Ecliptic 
@test minimum(dirs[:, 1]) ≈ 90 - 15
@test maximum(dirs[:, 1]) ≈ 90 + 15

# Note that 
@test dirs[1, 2] ≈ 0.0
@test dirs[end, 2] ≈ 36

#######################################################################
# Test case

dirs = rad2deg.(genpointings(0:1:60, [0, 0, 1], 0.0,
                             borespinang=0,
                             spinsunang=deg2rad(15),
                             spinrpm=0,
                             precrpm=0,
                             yearlyrpm=0.1)[1])
for idx in 1:size(dirs, 1)
    @test dirs[idx, 1] ≈ 90 - 15
end

#######################################################################
# Test case

dirs = rad2deg.(genpointings(0:1:60, [0, 0, 1], 0.0,
                             borespinang=0,
                             spinsunang=deg2rad(15),
                             spinrpm=0,
                             precrpm=1,
                             yearlyrpm=0.1)[1])
@test minimum(dirs[:, 1]) ≈ 90 - 15
@test maximum(dirs[:, 1]) ≈ 90 + 15

#######################################################################
# Test case

#######################################################################
# Test case

#######################################################################
# Test case
