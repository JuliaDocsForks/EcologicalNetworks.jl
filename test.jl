include("./src/EcologicalNetwork.jl")
using EcologicalNetwork
using Traceur
using StatsBase
using NamedTuples
using StatPlots
using Base.Test

N = web_of_life("A_HP_002")
Traceur.@trace connectance(N)
Traceur.@trace degree(N, 1)

Traceur.@trace each_species_its_module(N)

N = convert(BinaryNetwork, web_of_life("M_PL_010"))

partitions = [brim(lp(N)...) for i in 1:40]
Qs = partitions .|> x -> Q(x...)
bestQ = maximum(Qs)
bestPos = first(filter(i -> Qs[i] == bestQ, eachindex(Qs)))
m = partitions[bestPos][2]

m1 = [m[s] for s in species(N,1)]
m2 = [m[s] for s in species(N,2)]

N.A[sortperm(m1), sortperm(m2)] |> x -> heatmap(x, c=:Greys, leg=false, frame=:none)
