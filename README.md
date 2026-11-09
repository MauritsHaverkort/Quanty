# Quanty

Quanty is a script language for quantum many-body problems written in second quantization. The user writes a Lua script in which operators, wave functions and spectra are ordinary objects, so that `H = a*OppA + b*OppB` or `psi * H * psi` does the numerical work; there is no input file, the script is the calculation. Quanty started as a code for core-level spectroscopy of correlated materials (XAS, XMCD, RIXS, NIXS and photoemission) and grew into a general toolkit for exact diagonalization, configuration interaction, ligand-field and Anderson impurity models, tight-binding and dynamical mean-field theory (DMFT). Fast routines in C (Lanczos, Green's functions, spectra) do the work, and Lua controls everything around them: loops, temperature averages, file output and plotting. User documentation, tutorials and a forum can be found at [www.quanty.org](https://www.quanty.org).

## Getting Quanty

- Ready-made binaries, documentation and tutorials are at [www.quanty.org](https://www.quanty.org).
- The stable releases (`master` branch) are at [git.quanty.org/Quanty/release](https://git.quanty.org/Quanty/release).
- The current development version (`develop` branch) is at [git.quanty.org/Quanty/develop](https://git.quanty.org/Quanty/develop).


## Questions and bug reports

All communication goes through the forum at [www.quanty.org](https://www.quanty.org): questions on how to set up a calculation, bug reports, suggestions for new features, and discussions about the code.

## Building

Quanty compiles as a single file and only needs BLAS and LAPACK. On Linux with GNU g++:

```bash
g++ LuaInterface/main.cpp -D LUA_USE_LINUX -DBLASunderscore -O3 -fopenmp -o Quanty -lblas -llapack
```

On macOS, `cmake -S . -B build && cmake --build build` does everything. See [BUILDING.md](BUILDING.md) for the Intel compiler, macOS details, universal binaries and the compiler flags to avoid.

## Contributing

If you intend to work on part of the code, please discuss this with Maurits on the forum first; he coordinates the different parts and tasks. [CONTRIBUTING.md](CONTRIBUTING.md) explains the branches, the coding style and how the code is organized. Please also read the [code of conduct](CODE_OF_CONDUCT.md). [SECURITY.md](SECURITY.md) explains why a Quanty script should be treated like any other program.

## License and citation

Quanty is dual-licensed under your choice of the Creative Commons Attribution 4.0 International License (CC BY 4.0) or the MIT License; see [LICENSE](LICENSE). If Quanty contributes to scientific work, please cite the appropriate Quanty publications. A current list can be found at [www.quanty.org](https://www.quanty.org) and can also be obtained by running Quanty.

---

All in all, I hope we made something useful. Have fun calculating and coding, be respectful, and we appreciate feedback and bug reports on the forum.

Best wishes,
Maurits
