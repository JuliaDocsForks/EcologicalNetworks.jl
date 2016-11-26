module EcologicalNetwork

# Dependencies
using StatsBase
using Cairo
using Combinatorics

VERSION < v"0.4-dev" && import Lexicon

export EcoNetwork,

    # General types for all bipartite / unipartite
    Bipartite, Unipartite,

    # Types
    BipartiteNetwork, UnipartiteNetwork,
    BipartiteProbaNetwork, UnipartiteProbaNetwork,
    BipartiteQuantiNetwork, UnipartiteQuantiNetwork,
    
    # Union types for all proba or deterministic
    ProbabilisticNetwork, DeterministicNetwork, QuantitativeNetwork,
    NonProbabilisticNetwork,

    # Richness
    richness,

    # Nestedness
    η,

    # Measures of centrality
    centrality_katz,

    # Links and connectances
    links, links_var, connectance, connectance_var,
    linkage_density, link_number,

    # Measures of degree
    degree_out, degree_in, degree_out_var, degree_in_var, degree, degree_var,
    specificity,

    # Expected number of species with 0 interactions
    species_has_no_successors, species_has_no_predecessors,
    species_is_free, free_species,

    # Matrix manipulation utilities
    make_unipartite, make_threshold, make_binary, make_bernoulli, nodiag,
    adjacency,

    # Probability algebra utilities
    checkprob, i_esp, i_var, a_var, m_var,

    # Null models and testing
    null1, null2, null3out, null3in, nullmodel,
    test_network_property,

    # Draw
    plot_network,

    # Modularity
    Partition, Q, Qr, label_propagation, modularity,

    # Paths
    number_of_paths,

    # Beta-diversity
    betadiversity,
    whittaker, sorensen, jaccard, gaston,
    williams, lande, ruggiero, hartekinzig,
    harrison,

    # Motifs
    motif_p, motif_v, count_motifs, motif, motif_var,

    # Data
    stony, mcmullen, ollerton, bluthgen, robertson

include(joinpath(".", "types.jl"))
include(joinpath(".", "centrality.jl"))
include(joinpath(".", "connectance.jl"))
include(joinpath(".", "degree.jl"))
include(joinpath(".", "free_species.jl"))
include(joinpath(".", "matrix_utils.jl"))
include(joinpath(".", "nestedness.jl"))
include(joinpath(".", "proba_utils.jl"))
include(joinpath(".", "nullmodels.jl"))
include(joinpath(".", "draw.jl"))
include(joinpath(".", "modularity.jl"))
include(joinpath(".", "paths.jl"))
include(joinpath(".", "betadiversity.jl"))
include(joinpath(".", "motifs.jl"))
include(joinpath(".", "data.jl"))
include(joinpath(".", "test.jl"))

end