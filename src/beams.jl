import Healpix
export std2fwhm, fwhm2std, gaussian_beam

@doc raw"""
    std2fwhm(stddev)

Convert a standard deviation into a Full Width Half Maximum (FWHM) value. The
measure unit for `stddev` can be arbitrary, as the transformation is linear.
See also `fwhm2std`.
"""
std2fwhm(stddev) = 2 * sqrt(2 * log(2)) * stddev

@doc raw"""
    std2fwhm(stddev)

Convert a standard deviation into a Full Width Half Maximum (FWHM) value. The
measure unit for `stddev` can be arbitrary, as the transformation is linear.
See also `std2fwhm`.
"""
fwhm2std(fwhm) = fwhm / (2 * sqrt(2 * log(2)))

@doc raw"""
    gaussian_beam!(beam_map::Healpix.Map{T,O}, angle_std_rad; normalization = 1.0) where {T,O <: Healpix.Order}

Save a map of a circual Gaussian beam in `beam_map`. The width of the beam is
provided through the parameter `angle_std_rad`, which is expressed in radians.
If you want to specify the FWHM instead, use `fwhm2std` like in the following
example:

````julia
import Healpix
nside = 1024
beam_map = Healpix.Map{Float64, RingOrder}(nside)

# Assume 2.5° of FWHM
gaussian_beam!(beam_map, 2.5 |> fwhm2std)
````

See also `gaussian_beam`.
"""
function gaussian_beam!(beam_map::Healpix.Map{T,O}, angle_std_rad; normalization = 1.0) where {T,O <: Healpix.Order}
    for pixidx in 1:beam_map.resolution.numOfPixels
        theta, _ = Healpix.pix2ang(beam_map, pixidx)
        beam_map[pixidx] = normalization * exp(-theta^2 / (2 * angle_std_rad))
    end
end

function gaussian_beam(nside::Integer, angle_std_rad; normalization = 1.0)
    beam_map = Healpix.Map{Float64,Healpix.RingOrder}(nside)
    gaussian_beam!(beam_map, angle_std_rad, normalization = normalization)

    beam_map
end

function gaussian_beam(res::Healpix.Resolution, angle_std_rad; normalization = 1.0)
    gaussian_beam(res.nside, angle_std_rad, normalization = normalization)
end

@doc raw"""
    gaussian_beam(res::Healpix.Resolution, angle_std_rad; normalization = 1.0) where {T,O <: Healpix.Order}
    gaussian_beam(nside::Integer, angle_std_rad; normalization = 1.0) where {T,O <: Healpix.Order}

Save a map of a circual Gaussian beam in `beam_map`. The width of the beam is
provided through the parameter `angle_std_rad`, which is expressed in radians.
If you want to specify the FWHM instead, use `fwhm2std` like in the following
example:

````julia
# Assume 2.5° of FWHM
beam_map = gaussian_beam(1024, 2.5 |> fwhm2std)
````

See also `gaussian_beam!`.
"""
gaussian_beam