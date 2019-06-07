# Pointing generation

One of the most basic facilities provided by CMBMissionSim is the generation of
pointing timelines. CMBMissionSim assumes that the spacecraft orbits around the
L2 Lagrangean point of the Sun-Earth system, and that it performs a monotonous
spin around some axis, optionally with some precession period. This is a quite
rigid scheme, but it is the one followed by Planck, CORE, and LiteBIRD.

The basic object used to create pointings is the [`ScanningStrategy`](@ref)
structure. You can build a `ScanningStrategy` object using one of its
constructors. The most basic one takes the following syntax:

```julia
ScanningStrategy(;
    spin_rpm = 0,
    prec_rpm = 0,
    yearly_rpm = 1  / (MINUTES_PER_DAY * DAYS_PER_YEAR),
    hwp_rpm = 0,
    spinsunang_rad = deg2rad(45.0),
    borespinang_rad = deg2rad(50.0),
)
```

Note that this constructor only takes keywords, and that the measurement units
used for them are different from the ones used in the fields in
`ScanningStrategy`. Units in the constructors have been chosen in order to be
easy to use, while units in `ScanningStrategy` allow to make computations more
efficient.

Another kind of constructor allows to load the scanning strategy definition
from a JSON file. It comes in two flavours:

```julia
# First version
open("my_scanning_strategy.json", "r") do inpf
    sstr = ScanningStrategy(inpf)
end

# Second version (easier)
sstr = ScanningStrategy("my_scanning_strategy.json")
```

These JSON files can be created using the function [`save`](@ref):

```julia
sstr = ScanningStrategy(spinsunang_rad = deg2rad(35.0))
save()
```

Once you have a [`ScanningStrategy`](@ref) object, you generate your pointings
through one of the following functions:

- `genpointings` and `genpointings!` generate a set of pointing directions that
  encompass a range of time; most of the time you will use this.
- `time2pointing` and `time2pointing!` generate one pointing direction at a
  time; they are useful if you are working on a system with limited memory
  resources.

The version with the `!` and the end save their result in a preallocated block,
while the other ones use the return value of the function. The former is useful
if you are using some strategy to pre-allocate memory in order to optimize
running times. For instance, the following code is not optimal, as
`genpointings` is re-creating the same result matrix over and over again:

```julia
const SECONDS_IN_A_HOUR = 3600.
sstr = ScanningStrategy()
start_time = 0.0
# Simulate an observation lasting 1000 hours
for hour_num in 1:1000
    pnt = genpointings(sstr, start_time:(start_time + SECONDS_IN_A_HOUR),
        Float64[0, 0, 1], 0.0)
    # Use "pnt" here
    # ...
end
```

The code below is more efficient, as the allocation is done only once:

```julia
const SECONDS_IN_A_HOUR = 3600.
sstr = ScanningStrategy()
start_time = 0.0

# Allocate this variable once and for all
pnt = Array{Float64}(undef, SECONDS_IN_A_HOUR, 6)

# Simulate an observation lasting 1000 hours
for hour_num in 1:1000
    genpointings!(sstr, start_time:(start_time + SECONDS_IN_A_HOUR),
        Float64[0, 0, 1], 0.0, pnt)
    # Use "pnt" here
    # ...
end
```


## `ScanningStrategy`

```@docs
ScanningStrategy
update_scanning_strategy
```

## Loading and saving scanning strategies

```@docs
load_scanning_strategy
to_dict
save
```

## Pointing generation

```@docs
time2pointing
time2pointing!
genpointings
genpointings!
```

## Utility functions

```@docs
rpm2angfreq
angfreq2rpm
period2rpm
rpm2period
```
