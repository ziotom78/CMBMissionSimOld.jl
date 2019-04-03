module CMBMissionSim

using LinearAlgebra

import Healpix
import Quaternions
using StaticArrays
using LinearAlgebra: dot, ×

export MINUTES_PER_DAY,
       DAYS_PER_YEAR,
       SOLSYSDIR_ECL_θ,
       SOLSYSDIR_ECL_φ,
       SOLSYSSPEED_M_S,
       SPEED_OF_LIGHT_M_S,
       T_CMB,
       SOLSYS_SPEED_VEC_M_S,
       rpm2angfreq,
       period2rpm,
       genpointings,
       dipoletemperature,
       PointingInfo,
       time2pointing

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
    
    PointingInfo(
        spin_rpm=0,
        prec_rpm=0,
        yearly_rpm = 1  / (MINUTES_PER_DAY * DAYS_PER_YEAR),
        hwp_rpm = 0,
        spinsunang_rad = deg2rad(45.0),
        borespinang_rad = deg2rad(50.0),
    ) = new(
        rpm2angfreq(spin_rpm),
        rpm2angfreq(prec_rpm),
        rpm2angfreq(yearly_rpm),
        rpm2angfreq(hwp_rpm),
        spinsunang_rad,
        borespinang_rad,
        Quaternions.qrotation(Float64[0, 1, 0], borespinang_rad),
        Quaternions.qrotation(Float64[0.0, 1.0, 0.0], π / 2 - spinsunang_rad)
    )
end

function time2pointing(pinfo::PointingInfo, time_s, dir, polangle)
    curpolang = mod2pi(polangle + 4 * pinfo.ω_hwp * time_s)
    # The polarization vector lies on the XY plane; if polangle=0 then
    # the vector points along the X direction at t=0.
    poldir = SVector(cos(curpolang), sin(curpolang), 0.0)
    
    q2 = Quaternions.qrotation(Float64[0.0, 0.0, 1.0], pinfo.ω_spin * time_s)
    q4 = Quaternions.qrotation(Float64[1.0, 0.0, 0.0], pinfo.ω_prec * time_s)
    q5 = Quaternions.qrotation(Float64[0.0, 0.0, 1.0], pinfo.ω_year * time_s)
    
    qtot = q5 * (q4 * (pinfo.q3 * (q2 * pinfo.q1)))
    rot = Quaternions.rotationmatrix(qtot)
    # Direction in the sky of the beam main axis
    curvec = rot * dir
    # Direction in the sky of the beam polarization axis
    poldir = rot * poldir
    
    # The North for a vector v is just -dv/dθ, as θ is the
    # colatitude and moves along the meridian
    (θ, ϕ) = Healpix.vec2ang(curvec...)
    northdir = SVector(-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ))
    
    cosψ = clamp(dot(northdir, poldir), -1, 1)
    crosspr = northdir × poldir
    sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
    ψ = atan(cosψ, sinψ)

    (θ, ϕ, ψ, curvec...)
end

function genpointings(pinfo::PointingInfo, timerange_s, dir, polangle; usedirs=true)
    dirs = Array{Float64}(undef, length(timerange_s), usedirs ? 2 : 3)
    ψ = Array{Float64}(undef, length(timerange_s))

    @inbounds for (idx, t) in enumerate(timerange_s)
        curθ, curϕ, curψ, curx, cury, curz = time2pointing(pinfo, t, dir, polangle)

        if usedirs
            dirs[idx, 1] = curθ
            dirs[idx, 2] = curϕ
        else
            dirs[idx, 1] = curx
            dirs[idx, 2] = cury
            dirs[idx, 3] = curz
        end
        ψ[idx] = curψ
    end

    dirs
end

function genpointings!(timerange_s, dir, polangle, dirs, ψ;
    spinsunang=deg2rad(45.0),
    borespinang=deg2rad(50.0),
    spinrpm=0.0,
    precrpm=0.0,
    yearlyrpm=1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR),
    hwprpm=0.0,
    usedirs=true)

    ω_spin = rpm2angfreq(spinrpm)
    ω_prec = rpm2angfreq(precrpm)
    ω_year = rpm2angfreq(yearlyrpm)
    ω_hwp = rpm2angfreq(hwprpm)
    
    # These quaternions never change
    q1 = Quaternions.qrotation(Float64[0, 1, 0], borespinang)
    q3 = Quaternions.qrotation(Float64[0.0, 1.0, 0.0], π / 2 - spinsunang)

    @inbounds for (idx, time_s) in enumerate(timerange_s)
        curpolang = mod2pi(polangle + 4 * ω_hwp * time_s)
        # The polarization vector lies on the XY plane; if polangle=0 then
        # the vector points along the X direction at t=0.
        poldir = SVector(cos(curpolang), sin(curpolang), 0.0)
        
        q2 = Quaternions.qrotation(Float64[0.0, 0.0, 1.0], ω_spin * time_s)
        q4 = Quaternions.qrotation(Float64[1.0, 0.0, 0.0], ω_prec * time_s)
        q5 = Quaternions.qrotation(Float64[0.0, 0.0, 1.0], ω_year * time_s)
        
        qtot = q5 * (q4 * (q3 * (q2 * q1)))
        rot = Quaternions.rotationmatrix(qtot)
        # Direction in the sky of the beam main axis
        curvec = rot * dir
        # Direction in the sky of the beam polarization axis
        poldir = rot * poldir
        
        # The North for a vector v is just -dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, ϕ) = Healpix.vec2ang(curvec...)
        northdir = SVector(-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ))
        
        cosψ = clamp(dot(northdir, poldir), -1, 1)
        crosspr = northdir × poldir
        sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
        ψ[idx] = atan(cosψ, sinψ)

        if usedirs
            dirs[idx, 1] = θ
            dirs[idx, 2] = ϕ
        else
            dirs[idx, :] = curvec
        end
    end
    
    (dirs, ψ)
end

"""
    genpointings(timerange_s, dir, polangle; spinsunang=deg2rad(45.0), borespinang=deg2rad(50.0), spinrpm=0, precrpm=0, yearlyrpm=1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR), hwprpm=0, usedirs=true)

Generate a set of pointing directions and angles for a given orientation
of the boresight beam. Depending on the value of "usedirs", the following
data will be returned:

- if `usedirs` is true (the default), return the tuple `(dirs, ψ)`, where
  `dirs` is a N×2 matrix containing the values of θ (first column) and ϕ
  (second column)
- if `usedirs` is false, the `dirs` element in the tuple `(dirs, ψ)` will
  be a N×3 matrix containing vectors of length one.
"""
function genpointings(timerange_s, dir, polangle;
    spinsunang=deg2rad(45.0),
    borespinang=deg2rad(50.0),
    spinrpm=0.0,
    precrpm=0.0,
    yearlyrpm=1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR),
    hwprpm=0.0,
    usedirs=true)
    
    dirs = Array{Float64}(undef, length(timerange_s), usedirs ? 2 : 3)
    ψ = Array{Float64}(undef, length(timerange_s))

    genpointings!(timerange_s, dir, polangle, dirs, ψ;
                  spinsunang=spinsunang,
                  borespinang=borespinang,
                  spinrpm=spinrpm,
                  precrpm=precrpm,
                  yearlyrpm=yearlyrpm,
                  hwprpm=hwprpm,
                  usedirs=usedirs)
end

# These data have been taken from Planck 2015 I, Table 1 (10.1051/0004-6361/201527101)
const SOLSYSDIR_ECL_θ = 1.7656051330336222
const SOLSYSDIR_ECL_φ = 2.9958842149922833
const SOLSYSSPEED_M_S = 370082.2332
const SPEED_OF_LIGHT_M_S = 2.99792458e8
const T_CMB = 2.72548

const SOLSYS_SPEED_VEC_M_S = SOLSYSSPEED_M_S * [sin(SOLSYSDIR_ECL_θ) * cos(SOLSYSDIR_ECL_φ),
                                                sin(SOLSYSDIR_ECL_θ) * sin(SOLSYSDIR_ECL_φ),
                                                cos(SOLSYSDIR_ECL_θ)]

"""
    dipoletemperature(ecldir; solsys_speed_vec_m_s=SOLSYS_SPEED_VEC_M_S)

Compute the temperature of the dipole along the direction "ecldir", which should
be a 3-element array containing the XYZ components of the pointing vector. The
vector must be expressed in Ecliptic coordinates.
"""
function dipoletemperature(ecldir; solsys_speed_vec_m_s=SOLSYS_SPEED_VEC_M_S)
    betavec = SOLSYS_SPEED_VEC_M_S / SPEED_OF_LIGHT_M_S
    gamma = 1 / sqrt(1 - dot(betavec, betavec))
    
    T_CMB * (1 / (gamma * (1 - dot(betavec, ecldir))) - 1)
end

end # module
