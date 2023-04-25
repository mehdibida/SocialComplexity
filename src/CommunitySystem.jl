using StructArrays


###Pair of maxiumm and current values (used for resources)
mutable struct BoundedVal
    bound::Float64 ##Bound of the value
    current::Float64 ##Current value
end

struct CommunitySystem
    ###Size of the spatial grid
    grid_size_a::Int64 #First coordinate
    grid_size_b::Int64 #Second coordinate

    span_of_control::Int64

    ###resources of each community
    com_resource::Matrix{BoundedVal}

    ##Gives the index of the comunity that is the Superior
    com_superior::Matrix{Tuple{Int64,Int64}}
    ##Gives the index of the community that is chief
    com_polity::Matrix{Tuple{Int64,Int64}}
    ##Set of communities having a subordinate community
    superior_set::Set{Tuple{Int64,Int64}}

    ##power of each polity
    pol_power::Dict{Tuple{Int64, Int64}, BoundedVal}
end

function recover_resource(cs::CommunitySystem, diff_decay)
    for b = 1:cs.grid_size_b
        for a = 1:cs.grid_size_a
            cs.com_resource[a, b].current = cs.com_resource[a, b].bound - 
                                            diff_decay * (cs.com_resource[a, b].bound - cs.com_resource[a, b].current)
        end
    end
end

###Payment of tribute
function pay_tribute(cs::CommunitySystem, prop_res)
    ###Communities at the bottom of the hierarchy are the next payers
    next_payers = Set((a, b) for a = 1:cs.grid_size_a for b = 1:cs.grid_size_b if (a, b) ∉ cs.superior_set)
    while (length(next_payers) > 1) | ((0, 0) ∉ next_payers)
        for p in next_payers
            ##Resource transfer
            cs.com_resource[cs.com_superior[p...]].currrent += prop_res * cs.com_resource[p...].current
            cs.com_resource[p...].current *= (1.0 - prop_res)
            ##Remove current payer from list and get next payer
            pop!(next_payers, p)
            push!(next_payers, cs.com_superior[p...])
        end
    end
    return cs
end


####Update the power of each polity
function update_power(cs::CommunitySystem)
    ##Reinitialize communities powers
    ###########################################Check if polities are up to date
    map!(x -> (x.curent = 0.0 ; x.bound = 0.0), values(cs.pol_power))
    ##Calculate power
    for b = 1:cs.grid_size_b
        for a = 1:cs.grid_size_a
            cs.pol_power[cs.com_polity[a, b]].current += cs.com_resource[a, b].current
            cs.pol_power[cs.com_polity[a, b]].bound += cs.com_resource[a, b].bound
        end
    end
    return cs
end


### Space is topologically a torus, so it has no borders
###Convert any pair of coordinates so that they are on grid
on_grid(cs::CommunitySystem, a::Int64, b::Int64) = (a % cs.grid_size_a) + 1, (b % cs.grid_size_b) + 1

### Neighbors are connected points on a hexagonal grid
### We use the following coordinate system designate points: a * (1, 0) + b * (cos(pi/3), sin(pi/3)) 
### where (1, 0) and (0, 1) are the canonical bese vectors of R^2
### Space is topologically a torus, so it has no borders
function get_comm_neighbors(cs::CommunitySystem, a::Int64, b::Int64)
    (1 <= a <= cs.grid_size_a) & (1 <= b <= cs.grid_size_b) && error("Incorrect community coordinates.")
    return on_grid(cs, a, b - 1), on_grid(cs, a + 1, b - 1),
    on_grid(cs, a - 1, b), on_grid(cs, a + 1, b),
    on_grid(cs, a - 1, b + 1), on_grid(cs, a, b + 1)
end

###Returns a dictionary with indices (2d indices) of polities as keys and a
###set of neighboring communities
function get_polities_neighbors(cs::CommunitySystem)
    ###Dictionary with polities as keys and set of communities as values
    pol_neigh_comm = Dict{Tuple{Int64,Int64},Set{Tuple{Int64,Int64}}}()
    ###Go through communities and get thei neighbors
    ###if their neighbors are from other polities, add those neighbors to the set
    ###of neighboring communities of polities
    for b = 1:cs.grid_size_b
        for a = 1:cs.grid_size_a
            n_ab = get_comm_neighbors(cs, a, b)
            for (n_a, n_b) in n_ab
                if cs.com_polity[a, b] != cs.com_polity[n_a, n_b]
                    if haskey(pol_neigh_comm, cs.com_polity[a, b])
                        push!(pol_neigh_comm[cs.com_polity[a, b]], (n_a, n_b))
                    else
                        pol_neigh_comm[cs.com_polity[a, b]] = Set{Tuple{Float64,Float64}}((n_a, n_b))
                    end
                end
            end
        end
    end
    return pol_neigh_comm
end

###Returns a dictionary with polities and their initial communities targets
###pol_niegh_comm is a dictionary with indices (2d indices) of polities as keys and a
###set of neighboring communities (i.e. belonging to other polities)
function polity_attack(cs::CommunitySystem, pol_neigh_comm)
    ###Dictionary with attacking polities as keys and initial target community
    pol_attck_target = Dict{Tuple{Int64,Int64},Set{Tuple{Int64,Int64}}}()
    for p in pol_neigh_comm
        ###Find weakest neighbor; initialize with 
        c_com = (0, 0) ##current target community
        t_pol = (0, 0) ##current target polity
        t_pol_pow = -1.0 ##current polity power
        for c_p_n in pol_neigh_comm[p]
            ###Case a weaker polity is found
            if (cs.communities.power[cs.com_polity[c_p_n...]] < t_pol_pow) || ##
               ((cs.communities.power[cs.com_polity[c_p_n...]] == t_pol_pow) & (rand() < 0.5)) ## case: powers are the same, choose randomly
                t_pol = cs.com_polity[c_p_n...]
                t_pol_pow = cs.communities.power[t_pol...]
                c_com = c_p_n
            ##Case the polity is already a target but subordinate community is wealthier
            elseif (t_pol == cs.com_polity[c_p_n...]) && 
                   ((cs.communities.resource[c_p_n...] > cs.communities.resource[c_com...]) | ##case considered communities has more resources
                    ((cs.communities.resource[c_p_n...] == cs.communities.resource[c_com...]) & (rand() < 0.5))) ##case: both communities resources are the same, choose randomly
                c_com = c_p_n
            end
        end
        ##polity p decides weather to attack or not

        pol_attck_target[p] = c_com
    end
end




