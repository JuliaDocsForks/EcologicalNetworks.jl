include("./src/EcologicalNetwork.jl")
using EcologicalNetwork
using Traceur
using StatsBase
using NamedTuples
using StatPlots

N = nz_stream_foodweb()[5]

# Trophic species
@progress for i in 1:(richness(N)-1)
    s1 = species(N)[i]
    for j in (i+1):richness(N)
        s2 = species(N)[j]
        same_in = N[s1,:] == N[s2,:]
        same_out = N[:,s1] == N[:,s2]
        if same_in  && same_out
            println(s1, "\t", s2, "\t", same_in, "/", same_out, "\t", same_in&same_out)
        end
    end
end
