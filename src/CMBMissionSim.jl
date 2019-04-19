module CMBMissionSim

using LinearAlgebra

import Healpix
import Quaternions
using StaticArrays
import Printf
using LinearAlgebra: dot, ×

export MINUTES_PER_DAY,
       DAYS_PER_YEAR,
       SOLSYSDIR_ECL_THETA,
       SOLSYSDIR_ECL_PHI,
       SOLSYSSPEED_M_S,
       SPEED_OF_LIGHT_M_S,
       T_CMB,
       SPEED_OF_LIGHT_M_S,
       PLANCK_H_MKS,
       BOLTZMANN_K_MKS,
       rpm2angfreq,
       period2rpm,
       genpointings,
       DipoleParameters,
       doppler_temperature,
       dipole_temperature,
       PointingInfo,
       time2pointing

# These functions are faster than Quaternions.qrotation2

function qrotation_x(theta)
    Quaternions.Quaternion(cos(theta / 2), sin(theta / 2), 0.0, 0.0, true)
end

function qrotation_y(theta)
    Quaternions.Quaternion(cos(theta / 2), 0.0, sin(theta / 2), 0.0, true)
end

function qrotation_z(theta)
    Quaternions.Quaternion(cos(theta / 2), 0.0, 0.0, sin(theta / 2), true)
end

"""
    qrotation_x(theta)
    qrotation_y(theta)
    qrotation_z(theta)

Return a `Quaternions.Quaternion` object representing a rotation
around the ``e_x``, ``e_y``, or ``e_z`` axis by an angle `theta` (in
radians).
"""
qrotation_x, qrotation_y, qrotation_z


"""
    rpm2angfreq(v)

Convert rotations per minute into angular frequency 2πν (in Hertz).
"""
rpm2angfreq(v) = 2π * v / 60

"""
    period2rpm(p)

Convert a period (time span) in seconds into a number of rotations per minute
"""
period2rpm(p) = 60 / p

const MINUTES_PER_DAY = 60 * 24
const DAYS_PER_YEAR = 365.25

struct PointingInfo
    ω_spin::Float64
    ω_prec::Float64
    ω_year::Float64
    ω_hwp::Float64
    spinsunang::Float64
    borespinang::Float64
    # First quaternion used in the rotation
    q1::Quaternions.Quaternion
    # Third quaternion used in the rotation
    q3::Quaternions.Quaternion
    
    PointingInfo(spin_rpm = 0,
        prec_rpm = 0,
        yearly_rpm = 1  / (MINUTES_PER_DAY * DAYS_PER_YEAR),
        hwp_rpm = 0,
        spinsunang_rad = deg2rad(45.0),
        borespinang_rad = deg2rad(50.0),
    ) = new(rpm2angfreq(spin_rpm),
        rpm2angfreq(prec_rpm),
        rpm2angfreq(yearly_rpm),
        rpm2angfreq(hwp_rpm),
        spinsunang_rad,
        borespinang_rad,
        qrotation_y(borespinang_rad),
        qrotation_y(π / 2 - spinsunang_rad))
end

function time2pointing!(pinfo::PointingInfo, time_s, dir, polangle, resultvec)
    curpolang = mod2pi(polangle + 4 * pinfo.ω_hwp * time_s)
    # The polarization vector lies on the XY plane; if polangle=0 then
    # the vector points along the X direction at t=0.
    poldir = SVector(cos(curpolang), sin(curpolang), 0.0)
    
    q2 = qrotation_z(pinfo.ω_spin * time_s)
    q4 = qrotation_x(pinfo.ω_prec * time_s)
    q5 = qrotation_z(pinfo.ω_year * time_s)
    
    qtot = q5 * (q4 * (pinfo.q3 * (q2 * pinfo.q1)))
    rot = Quaternions.rotationmatrix(qtot)
    # Direction in the sky of the beam main axis
    resultvec[4:6] = rot * dir
    # Direction in the sky of the beam polarization axis
    poldir = rot * poldir
    
    # The North for a vector v is just -dv/dθ, as θ is the
    # colatitude and moves along the meridian
    (θ, ϕ) = Healpix.vec2ang(resultvec[4:6]...)
    northdir = SVector(-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ))
    
    cosψ = clamp(dot(northdir, poldir), -1, 1)
    crosspr = northdir × poldir
    sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
    resultvec[3] = atan(cosψ, sinψ)

    resultvec[1], resultvec[2] = θ, ϕ
end

function time2pointing(pinfo::PointingInfo, time_s, dir, polangle)
    resultvec = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    time2pointing!(pinfo, time_s, dir, polangle, resultvec)
    resultvec
end

function genpointings!(pinfo::PointingInfo, timerange_s, dir, polangle, dirs)
    @assert size(dirs)[1] == length(timerange_s)

    @inbounds for (idx, t) in enumerate(timerange_s)
        time2pointing!(pinfo, t, dir, polangle, view(dirs, idx, :))
    end
end

function genpointings(pinfo::PointingInfo, timerange_s, dir, polangle)
    dirs = Array{Float64}(undef, length(timerange_s), 6)
    genpointings!(pinfo, timerange_s, dir, polangle, dirs)

    dirs
end

function genpointings!(timerange_s, dir, polangle, dirs;
    spinsunang = deg2rad(45.0),
    borespinang = deg2rad(50.0),
    spinrpm = 0.0,
    precrpm = 0.0,
    yearlyrpm = 1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR),
    hwprpm = 0.0
)

    ω_spin = rpm2angfreq(spinrpm)
    ω_prec = rpm2angfreq(precrpm)
    ω_year = rpm2angfreq(yearlyrpm)
    ω_hwp = rpm2angfreq(hwprpm)
    
    # These quaternions never change
    q1 = qrotation_y(borespinang)
    q3 = qrotation_y(π / 2 - spinsunang)

    @inbounds for (idx, time_s) in enumerate(timerange_s)
        curpolang = mod2pi(polangle + 4 * ω_hwp * time_s)
        # The polarization vector lies on the XY plane; if polangle=0 then
        # the vector points along the X direction at t=0.
        poldir = SVector(cos(curpolang), sin(curpolang), 0.0)

        q2 = qrotation_z(ω_spin * time_s)
        q4 = qrotation_y(ω_prec * time_s)
        q5 = qrotation_z(ω_year * time_s)
        
        qtot = q5 * (q4 * (q3 * (q2 * q1)))
        rot = Quaternions.rotationmatrix(qtot)
        # Direction in the sky of the beam main axis
        dirs[idx, 4:6] = rot * dir
        # Direction in the sky of the beam polarization axis
        poldir = rot * poldir
        
        # The North for a vector v is just -dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, ϕ) = Healpix.vec2ang(dirs[idx, 4:6]...)
        dirs[idx, 1] = θ
        dirs[idx, 2] = ϕ

        northdir = SVector(-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ))
        
        cosψ = clamp(dot(northdir, poldir), -1, 1)
        crosspr = northdir × poldir
        sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
        dirs[idx, 3] = atan(cosψ, sinψ)
    end
end

"""
    genpointings(timerange_s, dir, polangle; spinsunang=deg2rad(45.0), borespinang=deg2rad(50.0), spinrpm=0, precrpm=0, yearlyrpm=1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR), hwprpm=0)

Generate a set of pointing directions and angles for a given orientation
of the boresight beam. The function returns a N×5 matrix containing the following fields:

1. The colatitude (in radians)
2. The longitude (in radians)
3. The X component of the one-length pointing vector
4. The Y component
5. The Z component
"""
function genpointings(timerange_s, dir, polangle;
    spinsunang = deg2rad(45.0),
    borespinang = deg2rad(50.0),
    spinrpm = 0.0,
    precrpm = 0.0,
    yearlyrpm = 1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR),
    hwprpm = 0.0)
    
    dirs = Array{Float64}(undef, length(timerange_s), 6)

    genpointings!(timerange_s, dir, polangle, dirs;
                  spinsunang = spinsunang,
                  borespinang = borespinang,
                  spinrpm = spinrpm,
                  precrpm = precrpm,
                  yearlyrpm = yearlyrpm,
                  hwprpm = hwprpm)

    dirs
end

# These data have been taken from Planck 2015 I, Table 1 (10.1051/0004-6361/201527101)
const SOLSYSDIR_ECL_THETA = 1.7656051330336222
const SOLSYSDIR_ECL_PHI = 2.9958842149922833
const SOLSYSSPEED_M_S = 370082.2332
const T_CMB = 2.72548

const SPEED_OF_LIGHT_M_S = 2.99792458e8
const PLANCK_H_MKS = 6.62606896e-34
const BOLTZMANN_K_MKS = 1.3806504e-23

@doc raw"""
This structure encodes the information needed to estimate the contribution of the
temperature of the CMB dipole caused by the motion of the Solar System with respect
to the CMB rest frame.

# Fields

- `solsysdir_theta_rad`: colatitude (in radians) of the dipole axis
- `solsysdir_phi_rad`: longitude (in radians) of the dipole axis
- `solsys_speed_m_s`: speed (in m/s) of the reference frame with respect to the
  CMB
- `solsys_velocity_m_s`: velocity 3-vector (in m/s) of the reference frame with
  respect to the CMB
- `tcmb_k`: thermodynamic temperature (in Kelvin) of the CMB

# Object creation

The simplest way to create a `DipoleParameters` object is to call the
constructor without arguments. In this case, the dipole axis will be provided in
Ecliptic coordinates, using the estimate by Planck 2015, and the best COBE
estimate for the CMB monopole temperature will be used. Otherwise, you can pass
any of the following keywords to set the parameters:

- `theta_rad` defaults to `SOLSYSDIR_ECL_THETA`
- `phi_rad` defaults to `SOLSYSDIR_ECL_PHI`
- `speed_m_s` defaults to `SOLSYSSPEED_M_S`
- `t_k` defaults to `T_CMB`

Here is an example where we produce a dipole that is 10% stronger than Planck's:

````julia
dip = DipoleParameters(speed_m_s = SOLSYS_SPEED_VEC_M_S * 1.10)
````
"""
struct DipoleParameters
    solsysdir_theta_rad::Float64
    solsysdir_phi_rad::Float64
    solsys_speed_m_s::Float64
    solsys_velocity_m_s::StaticArrays.StaticVector
    tcmb_k::Float64

    DipoleParameters(; 
        theta_rad = SOLSYSDIR_ECL_THETA, 
        phi_rad = SOLSYSDIR_ECL_PHI, 
        speed_m_s = SOLSYSSPEED_M_S, 
        t_k = T_CMB) = new(theta_rad, 
        phi_rad, 
        speed_m_s, 
        SVector{3}(speed_m_s * [
            sin(theta_rad) * cos(phi_rad),
            sin(theta_rad) * sin(phi_rad),
            cos(theta_rad),
        ]),
        t_k)
end

function Base.show(io::IO, params::DipoleParameters)
    if get(io, :compact, false)
        Printf.@printf(io, "DipoleParameters([%.4f, %.4f, %.4f], %.4f)",
            params.solsys_velocity_m_s[1],
            params.solsys_velocity_m_s[2],
            params.solsys_velocity_m_s[3],
            params.tcmb_k,
        )
    else
        Printf.@printf(io, """DipoleParameters(
    theta_rad=%.6e,
    phi_rad=%.6e,
    speed_m_s=%.6e,
    t_k=%.6e,
)""",
            params.solsysdir_theta_rad,
            params.solsysdir_phi_rad,
            params.solsys_speed_m_s,
            params.tcmb_k,
        )
    end
end

function doppler_temperature(velocity_m_s, dir, tcmb_k)
    betavec = velocity_m_s / SPEED_OF_LIGHT_M_S
    gamma = 1 / sqrt(1 - dot(betavec, betavec))
    
    tcmb_k * (1 / (gamma * (1 - dot(betavec, dir))) - 1)
end

function doppler_temperature(velocity_m_s,
    dir,
    tcmb_k, freq_hz,
)
    fact = PLANCK_H_MKS * freq_hz / (BOLTZMANN_K_MKS * tcmb_k)
    expfact = exp(fact)
    q = (fact / 2) * (expfact + 1) / (expfact - 1)

    betavec = velocity_m_s / SPEED_OF_LIGHT_M_S
    dotprod = dot(betavec, dir)
    tcmb_k * (dotprod + q * dotprod^2)
end

@doc raw"""
    doppler_temperature(velocity_m_s, dir, tcmb_k)
    doppler_temperature(velocity_m_s, dir, tcmb_k, freq_hz)

Compute the temperature caused by the motion by `velocity_m_s` (a 3-vector, in
m/s) through the Doppler effect. If `freq_hz` is specified, the computation
includes quadrupolar relativistic corrections.

# See also

If you need to compute the solar dipole caused by the motion of the Solar System
with respect to the CMB rest frame, it is easier to use
[`dipole_temperature`](@ref).

# Example

````jldoctest
julia> 

julia> 
````
"""
doppler_temperature

################################################################################

function dipole_temperature(dir;
    params::DipoleParameters = DipoleParameters(),)
    doppler_temperature(params.solsys_velocity_m_s, dir, params.tcmb_k)
end

function dipole_temperature(dir,
    freq_hz::Number; 
    params::DipoleParameters = DipoleParameters(),)
    doppler_temperature(params.solsys_velocity_m_s, dir, params.tcmb_k, freq_hz)
end

@doc raw"""
    dipole_temperature(dir; params::DipoleParameters = DipoleParameters())
    dipole_temperature(dir, freq_hz; params::DipoleParameters = DipoleParameters())
    
Compute the temperature of the dipole along the direction `dir`, which should be
a 3 - element array containing the XYZ components of the pointing vector.The
vector must be expressed in the same coordinates as the vector used to specify
the dipole in `params`, which is an object of type `DipoleParameters`.(Hint:if
you created this object calling the constructor without arguments, as in
`DipoleParameters()`, Ecliptic coordinates will be used.)

If `freq_hz` is specified, a relativistic frequency - dependent correction is
applied.

# See also

If you need to compute the temperature caused by any other kinetic component
(e.g., the orbital dipole caused by the motion of the spacecraft), use
[`doppler_temperature`](@ref).

# Example

````jldoctest
julia> dipole_temperature([0, 0, 1])
-0.0006532169921239991

julia> dipole_temperature([0, 0, 1], params=DipoleParameters(speed_m_s=371e3))
-0.0006548416782232279

julia> dipole_temperature([0, 0, 1], 30e9, params=DipoleParameters(speed_m_s=371e3))
-0.0006527515338213527
````
"""
dipole_temperature

end # module
