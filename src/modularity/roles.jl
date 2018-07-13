"""
    modular_roles(N::T, m::Dict{K,Int64}) where {T <: BinaryNetwork, K <: AllowedSpeciesTypes}

Returns the modular roles of species in a network. The returned object is a
dictionary, where values are a tuple with `(z,c)`, where `z` is the connectivity
*within* a species' module, and `c` is the connectivity *across* modules,
including a species' own.
"""
function modular_roles(N::T, m::Dict{K,Int64}) where {T <: BinaryNetwork, K <: AllowedSpeciesTypes}

    z_scores = Dict([s => 0.0 for s in species(N)])
    c_scores = Dict([s => 0.0 for s in species(N)])

    modules = unique(collect(values(m)))

    for mod in modules
        modN = extract_sub_module(N, m, mod)
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

function extract_sub_module(N::BipartiteNetwork, M::Dict{T,Int64}, m::Int64) where {T <: AllowedSpeciesTypes}
    modsp = collect(keys(filter((k,v) -> v == m, M)))
    s1m = filter(x -> x ∈ species(N,1), modsp)
    s2m = filter(x -> x ∈ species(N,2), modsp)
    modN = N[s1m, s2m]
    return modN
end

function extract_sub_module(N::UnipartiteNetwork, M::Dict{T,Int64}, m::Int64) where {T <: AllowedSpeciesTypes}
    modsp = collect(keys(filter((k,v) -> v == m, M)))
    sm = filter(x -> x ∈ species(N,1), modsp)
    modN = N[sm]
    return modN
end
