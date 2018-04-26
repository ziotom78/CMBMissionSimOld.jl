module CMBMissionSim

import Healpix
import Quaternions

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
       dipoletemperature

doc"""
    rpm2angfreq(v)

Convert rotations per minute into angular frequency 2πν (in Hertz).
"""
rpm2angfreq(v) = 2 * π * v / 60.0

doc"""
    period2rpm(p)

Convert a period (time span) in seconds into a number of rotations per minute
"""
period2rpm(p) = 60.0 / p

const MINUTES_PER_DAY = 60 * 24
const DAYS_PER_YEAR = 365.25

doc"""
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
                      spinrpm=0,
                      precrpm=0,
                      yearlyrpm=1.0 / (MINUTES_PER_DAY * DAYS_PER_YEAR),
                      hwprpm=0,
                      usedirs=true)
    
    ω_spin = rpm2angfreq(spinrpm)
    ω_prec = rpm2angfreq(precrpm)
    ω_year = rpm2angfreq(yearlyrpm)
    ω_hwp = rpm2angfreq(hwprpm)
    
    dirs = Array{Float64}(length(timerange_s), usedirs ? 2 : 3)
    ψ = Array{Float64}(length(timerange_s))

    for (idx, time_s) in enumerate(timerange_s)
        curpolang = mod2pi(polangle + 4 * ω_hwp * time_s)
        # The polarization vector lies on the XY plane; if polangle=0 then
        # the vector points along the X direction at t=0.
        poldir = [cos(curpolang), sin(curpolang), 0]
        
        q1 = Quaternions.qrotation([0, 1, 0], borespinang)
        q2 = Quaternions.qrotation([0, 0, 1], ω_spin * time_s)
        q3 = Quaternions.qrotation([0, 1, 0], π / 2 - spinsunang)
        q4 = Quaternions.qrotation([1, 0, 0], ω_prec * time_s)
        q5 = Quaternions.qrotation([0, 0, 1], ω_year * time_s)
        
        qtot = q5 * (q4 * (q3 * (q2 * q1)))
        rot = Quaternions.rotationmatrix(qtot)
        
        # Direction in the sky of the beam main axis
        curvec = rot * dir
        # Direction in the sky of the beam polarization axis
        curpol = rot * poldir
        
        # The North for a vector v is just -dv/dθ, as θ is the
        # colatitude and moves along the meridian
        (θ, φ) = Healpix.vec2ang(curvec...)
        northdir = [cos(θ) * cos(φ); cos(θ) * sin(φ); sin(θ)]
        
        cosψ = clamp(dot(northdir, curpol), -1, 1)
        crosspr = northdir × curpol
        sinψ = clamp(dot(crosspr, crosspr), -1, 1)
        ψ[idx] = atan2(cosψ, sinψ)

        if usedirs
            dirs[idx, :] = [θ, φ]
        else
            dirs[idx, :] = curvec
        end
    end
    
    (dirs, ψ)
end

const SOLSYSDIR_ECL_θ = 1.7656131194951572
const SOLSYSDIR_ECL_φ = 2.995889600573578
const SOLSYSSPEED_M_S = 370082.2332
const SPEED_OF_LIGHT_M_S = 2.99792458e8
const T_CMB = 2.72548

const SOLSYS_SPEED_VEC_M_S = SOLSYSSPEED_M_S * [sin(SOLSYSDIR_ECL_θ) * cos(SOLSYSDIR_ECL_φ),
                                                sin(SOLSYSDIR_ECL_θ) * sin(SOLSYSDIR_ECL_φ),
                                                cos(SOLSYSDIR_ECL_θ)]

doc"""
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
