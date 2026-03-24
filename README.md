[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/lennoxkeeble/NeutronStarOscillations/blob/main/LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/lennoxkeeble/NeutronStarOscillations)
[![DOI](https://zenodo.org/badge/1136267109.svg)](https://doi.org/10.5281/zenodo.19207181)
[![Build Status](https://github.com/lennoxkeeble/NeutronStarOscillations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lennoxkeeble/NeutronStarOscillations.jl/actions/workflows/CI.yml?query=branch%3Amain)

# NeutronStarOscillations
Julia code for solving the equations of motion governing linearized radial perturbations of spherically symmetric neutron stars in the frequency and time domains. Methods are provided for inviscid perfect fluid (PF) stars and viscous stars within the frameworks of Eckart and Bemfica-Disconzi-Noronha-Kovtun (BDNK) hydrodynamics. Frequency-domain methods are provided for PF and Eckart stars. Time-domain methods are provided for all three fluid models. Currently, the code is restricted to cold, polytropic neutron stars with an equation of state of the form $p=\kappa \epsilon^{1+1/n}$. Example code is provided in the /src/Examples/ directory. For details about the numerical methods employed, see [arxiv:](arxiv.org).

[Code in constant development.]

## Citation 

If you use this code in your work, please cite the following reference:

```
@article{Keeble:2026
    author = "Keeble, Lennox and Redondo-Yuste, Jaime",
    title = "{}",
    eprint = "26",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    year = "2026"
}
```

## Limitations and known performance issues ##

* This code has been developed and tested on Julia version 1.12.0 on macOS (x86_64-apple-darwin24.0.0).
  
## Authors ##

- [Lennox S. Keeble](https://lennoxkeeble.github.io)
- [Jaime Redondo-Yuste](https://jredondoyuste.github.io)

## MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.