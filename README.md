![Amburana](recipe/amburana.png)

Amburana generates sketches of genomic sequences using weighted and unweighted minhash signatures, and calculates
distance matrices using several pre-determined sets of k-mer sizes. 
It also implements several clustering algorithms which are used in these pairwise distance matrices. 
The algorithms implemented include hierarchical clustering (UPGMA, WPGMA, median, etc.), OPTICS, and affinity
propagation.

*Amburana* is a genus of endangered South American trees which are used in the production of cachaça casks.
## License 
Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

Amburana is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

This software relies on the low-level library [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib), 
which is defined as a submodule &mdash; so don't forget to git recursively. 
