This repository hosts the source code for the paper:

> *Designing Molecular RNA Switches with Restricted Boltzmann Machines*
> 
> by Jorge Fernandez-de-Cossio-Diaz, Pierre Hardouin, Francois-Xavier Lyonnet du Moutier, Andrea Di Gioacchino, Bertrand Marchand, Yann Ponty, Bruno Sargueil, Rémi Monasson, Simona Cocco
> 
> bioRxiv preprint: https://www.biorxiv.org/content/10.1101/2023.05.10.540155v2.abstract

If you use this code, please cite this paper (or you can use the included [CITATION.bib](https://github.com/cossio/SamApp2024.jl/blob/main/CITATION.bib)).

# Setup

The code is organized as a [Julia](https://julialang.org) package. If you haven't setup Julia on your computer, you can install it following the simple instructions from the [official website](https://julialang.org/downloads/). We have tested our package on the latest Julia version at the time of writing (`v1.10.4`). It should work on future non-breaking versions (`v1.x.y`) of Julia as well.

After having installed Julia, clone this repository locally using `git`, on a local directory, by running the following command from the terminal:

```bash
git clone https://github.com/cossio/SamApp2024.jl.git
```

Then, navigate to this directory in your terminal, start Julia (by running `julia` in the terminal), and then activate the included `Project.toml` by running the command:

```julia
julia> import Pkg
julia> Pkg.activate(pwd()) # activate the project in the current directory
julia> Pkg.instantiate() # install all the dependencies
julia> import SamApp2024 # import the package
```

The second to last command will install all the dependencies of this package. The last command, imports the package, so you can use it in the Julia REPL. The first time you import it, it might some time for precompilation (but future imports should be faster).

At this point you can call any of the functions in the package.

## Running the Pluto notebooks

The included [Pluto notebooks](https://github.com/cossio/SamApp2024.jl/tree/main/pluto) contain the code for the downstream analysis done in the paper (which depend on the `SamApp2024` package). For more information about Pluto notebooks, see: https://plutojl.org.

To setup Pluto, start a fresh Julia session, and then run the following commands:

```julia
julia> import Pkg
julia> Pkg.add("Pluto")
julia> import Pluto; Pluto.run()
```

This starts a Pluto server. Now open your browser, and navigate to `localhost:1234`, which is the default address where Pluto is listening. You can then open the notebooks by navigating to the `pluto` directory in the repository, and opening the `.jl` files.

# Related packages 

The code implementing Restricted Boltzmann machines (training, sampling, and other functions) is provided in a separate package: https://github.com/cossio/RestrictedBoltzmannMachines.jl. Note that this package will be installed automatically as a dependency of this repository by the Julia package manager.

The code for simulations of other riboswitch families (reported in Supplementary Materials of the paper cited above) can be found at the following repository: https://github.com/cossio/RiboswitchesSimulations20240620.jl.

# Issues

If you encounter any problems, please [open an issue](https://github.com/cossio/SamApp2024.jl/issues/new/choose) in this repository.