```@meta
DocTestSetup = quote
    using CMBMissionSim
    const CMBSim = CMBMissionSim
end
```

# CMBMissionSim User's Manual

A simulator for space CMB missions.

To install Stripeline, start Julia and type the following command:
```julia
using Pkg
Pkg.add("https://github.com/ziotom78/CMBMissionSim")
```

To run the test suite, type the following command:
```julia
using Pkg; Pkg.test("CMBMissionSim")
```

In this manual, we will often assume that `CMBMissionSim` has been imported
using the following commands:

```julia
import CMBMissionSim
const CMBSim = CMBMissionSim
```

## Documentation

The documentation was built using
[Documenter.jl](https://github.com/JuliaDocs).

```@example
import Dates: now #hide
println("Documentation built $(now()) with Julia $(VERSION).") # hide
```

## Index

```@index
```
