import Healpix
import Quaternions
import LinearAlgebra: dot, ×

export MINUTES_PER_DAY,
       DAYS_PER_YEAR,
       rpm2angfreq,
       period2rpm,
       genpointings,
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
    hwprpm = 0.0)

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
        poldir = StaticArrays.SVector(cos(curpolang), sin(curpolang), 0.0)

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

        northdir = StaticArrays.SVector(-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ))
        
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
