# Dynamical Diffraction using Wie algorithm

Mostly written and maintained by [Eric Landahl, DePaul University Physics Dept.](https://sites.google.com/site/elandahl/)

## Run benchmark
Wie_Adapt_Test.m will compare a simple strained crystal to a GID result

## For help
In MatLab or Octave type `help filename`

## Update 12/7/2017
Fixed incompatibilities with standard MatLAB functions
Wrote README.md and added To Do List

## To Do List
1. Handle a few simple bulk crystals over a reasonable energy range (7 - 15 keV)
  * GaAs
  * Si
  * InSb
  * Ge
2. Develop an easy to use driver to TRXD calling various strain functions
3. Develop simple strain functions
  * Thermal diffusion
  * Simple mean free path modified thermal diffusion model
  * Electron diffusion with Auger and radiative decay
  * Thomsen model
  * Arbitrary strain propogation (e.g. Diffusion drives strain)
  * Noninear absorption (e.g. Saturable Absorption, Two Photon Absorption)
4. Transverse and shear strain, benchmarked to GID
5. More sophisticated models for nanoscale thermal and electrical transport


