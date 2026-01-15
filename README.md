# Cartan

[![Build Status](https://github.com/ooalshei/Cartan.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ooalshei/Cartan.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://app.travis-ci.com/ooalshei/Cartan.jl.svg?branch=main)](https://app.travis-ci.com/ooalshei/Cartan.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ooalshei/Cartan.jl?svg=true)](https://ci.appveyor.com/project/ooalshei/Cartan-jl)
[![Coverage](https://codecov.io/gh/ooalshei/Cartan.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ooalshei/Cartan.jl)
[![Coverage](https://coveralls.io/repos/github/ooalshei/Cartan.jl/badge.svg?branch=main)](https://coveralls.io/github/ooalshei/Cartan.jl?branch=main)

## Documentation

Full API documentation and examples are available in the `docs/` folder and
can be built with Documenter.jl. The algorithms implemented follow:

- "Cartan decompositions for Pauli operator algebras", arXiv:2512.06070

To build the docs locally:

```julia
]
activate .
cd("docs")
using Pkg; Pkg.instantiate()
include("make.jl")
```
