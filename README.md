This repository hosts the source code for the paper:

> *Designing Molecular RNA Switches with Restricted Boltzmann Machines*
> 
> by Jorge Fernandez-de-Cossio-Diaz, Pierre Hardouin, Francois-Xavier Lyonnet du Moutier, Andrea Di Gioacchino, Bertrand Marchand, Yann Ponty, Bruno Sargueil, Rémi Monasson, Simona Cocco
> 
> bioRxiv preprint: https://www.biorxiv.org/content/10.1101/2023.05.10.540155v2.abstract

If you use this code, please cite this paper (or you can use the included [CITATION.bib](https://github.com/cossio/SamApp2024.jl/blob/main/CITATION.bib)).

# Usage

The code is organized as a [Julia](https://julialang.org) package. After having installed Julia, clone this repository locally and activate the included `Project.toml`.

The included [Jupyter notebooks](https://github.com/cossio/SamApp2024.jl/tree/main/notebooks) and [Pluto notebooks](https://github.com/cossio/SamApp2024.jl/tree/main/pluto) contain the all the code for the analysis done in the paper (for more information about Pluto notebooks, see: https://plutojl.org).

The code implementing Restricted Boltzmann machines (training, sampling, and other functions) is provided in a separate package: https://github.com/cossio/RestrictedBoltzmannMachines.jl. Note that this package will be installed automatically as a dependency.