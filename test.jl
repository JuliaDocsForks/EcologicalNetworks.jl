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

function modular_roles(N::BipartiteNetwork, m::Dict{T,Int64}) where {T <: AllowedSpeciesTypes}

    z_scores = Dict([s => 0.0 for s in species(N)])
    c_scores = Dict([s => 0.0 for s in species(N)])

    modules = unique(collect(values(m)))

    for mod in modules
        modsp = collect(keys(filter((k,v) -> v == mod, m)))
        s1m = filter(x -> x ∈ species(N,1), modsp)
        s2m = filter(x -> x ∈ species(N,2), modsp)
        modN = N[s1m, s2m]
        degN = degree(modN)
        km, ksd = mean(values(degN)), std(values(degN))
        for s in species(modN)
            z_scores[s] = ksd == 0.0 ? 0.0 : (degN[s] - km)/ksd
        end
    end

    c_temp = Dict([s => zeros(Int64, length(modules)) for s in species(N)])
    for i in interactions(N)
        f = i.from
        t = i.to
        c_temp[f][m[t]] += 1
        c_temp[t][m[f]] += 1
    end

    for s in species(N)
        mod_prof = c_temp[s]
        c_scores[s] = 1.0-sum((mod_prof./sum(mod_prof)).^2)
    end

    roles = Dict([s => (z_scores[s], c_scores[s]) for s in species(N)])

    return roles
end
