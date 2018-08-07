"""
    brim(N::NT, L::Dict{E,Int64}) where {NT<:AbstractEcologicalNetwork,E<:AllowedSpeciesTypes}

Uses BRIM to optimize the modularity of an ecological network. The `L` argument
is a dictionary mapping every species in the network to its module. This
function returns a tuple of the network and its module assignment.
"""
function brim(N::NT, L::Dict{E,Int64}) where {NT<:AbstractEcologicalNetwork,E<:AllowedSpeciesTypes}
  @assert all(species(N) .∈ keys(L))
  EcologicalNetworks.tidy_modules!(L)

  old_Q = Q(N, L)
  new_Q = old_Q+0.00001

  m = links(N)

  c = length(unique(collect(values(L))))

  R = zeros(Int64, (richness(N; dims=1),c))
  for si in eachindex(species(N; dims=1))
    R[si,L[species(N; dims=1)[si]]] = 1
  end

  T = zeros(Int64, (richness(N; dims=2),c))
  for sj in eachindex(species(N; dims=2))
    T[sj,L[species(N; dims=2)[sj]]] = 1
  end

  d_in, d_out = degree(N; dims=2), degree(N; dims=1)
  s_d_in = [d_in[s] for s in species(N; dims=2)]
  s_d_out = [d_out[s] for s in species(N; dims=1)]
  P = kron(s_d_out, s_d_in')./m
  B = N.A .- P

  while old_Q < new_Q

    t_tilde = B*T
    R = map(Int64, t_tilde .== maximum(t_tilde; dims=2))
    r_tilde = B'*R
    T = map(Int64, r_tilde .== maximum(r_tilde; dims=2))
    S = vcat(R, T)
    L = Dict(zip(species(N), vec(mapslices(r -> StatsBase.sample(findall(x->x!=0, r)), S; dims=2))))
    EcologicalNetworks.tidy_modules!(L)
    old_Q = new_Q
    new_Q = Q(N,L)

  end

  EcologicalNetworks.tidy_modules!(L)
  return (N, L)

end
