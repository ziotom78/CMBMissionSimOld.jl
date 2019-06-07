# Check of handedness
@test CMBSim.rotationmatrix_normalized(CMBSim.qrotation_x(π/2)) * [0, 1, 0] ≈ [0, 0, 1]
@test CMBSim.rotationmatrix_normalized(CMBSim.qrotation_y(π/2)) * [1, 0, 0] ≈ [0, 0, -1]
@test CMBSim.rotationmatrix_normalized(CMBSim.qrotation_z(π/2)) * [1, 0, 0] ≈ [0, 1, 0]

# Simple composition of rotations
@test CMBSim.qrotation_x(π/3) * CMBSim.qrotation_x(π/3) ≈ CMBSim.qrotation_x(2π/3)
@test CMBSim.qrotation_y(π/3) * CMBSim.qrotation_y(π/3) ≈ CMBSim.qrotation_y(2π/3)
@test CMBSim.qrotation_z(π/3) * CMBSim.qrotation_z(π/3) ≈ CMBSim.qrotation_z(2π/3)

# More complex compositions
comp = CMBSim.qrotation_y(π/2) * CMBSim.qrotation_x(π/2)
@test CMBSim.rotationmatrix_normalized(comp) * [0, 1, 0] ≈ [1, 0, 0]

comp = CMBSim.qrotation_z(π/2) * CMBSim.qrotation_y(π/2)
@test CMBSim.rotationmatrix_normalized(comp) * [0, 0, 1] ≈ [0, 1, 0]
