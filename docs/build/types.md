
<a id='Informations-about-types-1'></a>

# Informations about types


As there are many ways to represent ecological networks, and the correct way to measure a given property varies in function of the representation, `EcologicalNetwork` has a series of built-in types.


<a id='Type-hierarchy-1'></a>

## Type hierarchy


The `EcologicalNetwork` package has four main types that are meant to be used by the user: `BipartiteNetwork`, `BipartiteProbaNetwork`, `UnipartiteNetwork`, and `UnipartiteProbaNetwork`. All types with `Proba` in their names are meant to represent probabilistic networks, and are also part of the union type `ProbabilisticNetwork`. All types *without* `Proba` in their name are part of the union type `DeterministicNetwork`, and represent networks in which interactions are either present or absent. All types starting with `Bipartite` are also part of the abstract type `Bipartite`, and types with `Unipartite` in their names are part of the abstract type `Unipartite`. Finally, both `Unipartite` and `Bipartite` are part of the abstract type `EcoNetwork`.


Although this may seem convoluted, this is important to understand: when writing functions, you can restrict them to any type of network you want by using the right type in their declaration. You can also check properties of a network just by looking at its type. For example, one can check whether a network `N` is bipartite with:


```julia
typeof(N) <: Bipartite
```


Networks are represented as two-dimensional matrices. All types are simply wrappers around an adjacency matrix, stored as the `A` property of the object. To look at the adjacency matrix of a network `N`, one therefore uses `N.A`. These matrices must be read as: the existence/probability of an interaction *from* the species of the *i*-th row *to* the species in the *j*-th column. This implies that the networks are, by default, directed.


<a id='Data-types-1'></a>

## Data types


Interactions in deterministic networks are represented as boolean (`true`/`false`) values. This is memory efficient, so large networks can be represented (one interactions represented as a boolean uses 8 times fewer memory than the same interaction represented as an integer). This being said, all networks of the `DeterministicNetwork` type can be read from matrices of integers, *as long as these matrices only contain 0 and 1*.


In probabilistic networks, interactions are stored as floating point (`Float64`) numbers. These values have to be between 0.0 and 1.0, as they represent probabilities.


Networks of the `Unipartite` class must have the same number of rows and columns. The species in the rows and columns are the same. Networks of the `Bipartite` class are expected to have different numbers of rows and columns, as the species in rows and columns are different species. It is possible to convert a network from `Bipartite` to `Unipartite` using the `make_unipartite` function:


```julia
B = BipartiteNetwork(rand(Bool, (3, 5)))
richnes(B) # 8
U = make_unipartite(B)
richness(U) # 8
```


The documentation for `make_unipartite` gives additional explanations about how the conversion is done. In the overwhelming majority of cases, applying any measure to a bipartite network, and to the same network made unipartite, should give the same results (connectance is one notable example).


<a id='Type-reference-1'></a>

## Type reference

<a id='EcologicalNetwork.Bipartite' href='#EcologicalNetwork.Bipartite'>#</a>
**`EcologicalNetwork.Bipartite`** &mdash; *Type*.



Bipartite type

All bipartite networks

<a id='EcologicalNetwork.BipartiteNetwork' href='#EcologicalNetwork.BipartiteNetwork'>#</a>
**`EcologicalNetwork.BipartiteNetwork`** &mdash; *Type*.



BipartiteNetwork

A bipartite deterministic network.

<a id='EcologicalNetwork.BipartiteNetwork-Tuple{Array{Int64,2}}' href='#EcologicalNetwork.BipartiteNetwork-Tuple{Array{Int64,2}}'>#</a>
**`EcologicalNetwork.BipartiteNetwork`** &mdash; *Method*.



Construct a bipartite network from a matrix of integer

<a id='EcologicalNetwork.EcoNetwork' href='#EcologicalNetwork.EcoNetwork'>#</a>
**`EcologicalNetwork.EcoNetwork`** &mdash; *Type*.



EcoNetwork type

This is an abstract type that allows to generate functions for all sorts of networks. All other types are derived from this one.

<a id='EcologicalNetwork.Unipartite' href='#EcologicalNetwork.Unipartite'>#</a>
**`EcologicalNetwork.Unipartite`** &mdash; *Type*.



Unipartite type

All unipartite networks

<a id='EcologicalNetwork.UnipartiteNetwork' href='#EcologicalNetwork.UnipartiteNetwork'>#</a>
**`EcologicalNetwork.UnipartiteNetwork`** &mdash; *Type*.



UnipartiteNetwork

An unipartite deterministic network.

<a id='EcologicalNetwork.UnipartiteNetwork-Tuple{Array{Int64,2}}' href='#EcologicalNetwork.UnipartiteNetwork-Tuple{Array{Int64,2}}'>#</a>
**`EcologicalNetwork.UnipartiteNetwork`** &mdash; *Method*.



Construct an unipartite network from a matrix of integer

