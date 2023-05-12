using DataStructures
using Random

###Type alias for CartesianIndex{2} to be used as coordinates for communities
const ComCoord = CartesianIndex{2}

####Big problem: the model is not flow consistent
####The resource accounting is flawed
#### It seems that there is no accumulation of resources

###Pair of maxiumm and current values (used for resources)
mutable struct BoundedVal
    max_lvl::Float64 ##Max resource
    cur_lvl::Float64 ##Current level
    ##Constructor
    BoundedVal(m, c) = c <= m ? new(m, c) : error("Current value cannot be greater than bound.")
end

function set_func_cur!(bv::BoundedVal, func::F, f_args...) where {F}
    cv = func(f_args...)
    if cv > bv.max_lvl
        error("Current value cannot be greater than bound.")
    else
        bv.cur_lvl = cv
    end
    return bv
end

function set_func_max!(bv::BoundedVal, func::F, f_args...) where {F}
    mv = func(f_args...)
    if mv < bv.cur_lvl
        error("Bound cannot be less than current value.")
    else
        bv.max_lvl = mv
    end
    return bv
end
##Getters
get_cur(bv::BoundedVal) = bv.cur_lvl
get_max(bv::BoundedVal) = bv.max_lvl
##Setters for maximum and current values
set_cur!(bv::BoundedVal, cv::Real) = set_func_cur!(bv, identity, cv)
set_max!(bv::BoundedVal, mv::Real) = set_func_max!(bv, identity, mv)
##Add a value to maximum and current levels
add_to_cur!(bv::BoundedVal, val::Real) = set_func_cur!(bv, +, val, bv.cur_lvl)
add_to_max!(bv::BoundedVal, val::Real) = set_func_max!(bv, +, val, bv.max_lvl)
##Multiply maximum and current levels by a value
mult_cur!(bv::BoundedVal, val::Real) = set_func_cur!(bv, *, val, bv.cur_lvl)
mult_max!(bv::BoundedVal, val::Real) = set_func_max!(bv, *, val, bv.max_lvl)


struct CommunitySystem ##A is grid_size_a, B is grid_size_b
    ###Size of the spatial grid
    grid_size_a::Int64 #First coordinate
    grid_size_b::Int64 #Second coordinate

    span_of_control::Int64

    power_importance::Float64
    will_to_attack::Float64
    cost_factor::Float64

    tax_rate::Float64

    ###resources of each community
    com_resource::Matrix{BoundedVal}
    ###amount of resources received by each community from tribute
    com_power::Matrix{BoundedVal}
    ###Cost of war (in proportion to the power)
    com_war_cost::Matrix{Float64}
    
    ##Gives the index of the comunity that is the Superior, superior of chiefs are themselves
    com_superior::Matrix{ComCoord}
    ##Gives the index of the community that is chief
    com_polity::Matrix{ComCoord}
    
    ##Gives the rank of each community
    ##Rank is defined as the max depth of the tree of vassals
    ##Communities with no vassals have rank 0
    com_rank::Matrix{Int64}

    com_list_cache::Vector{ComCoord}
    ##Power of subchiefs, i.e. first order vassals (only current power is taken into account)
    #subc_power::Dict{Tuple{Int64, Int64}, Float64}

    ##power of each polity
    #pol_power::Dict{Tuple{Int64, Int64}, Resource}
end


function CommunitySystem(size_A::Integer, size_B::Integer, span_control, pow_importance, will_attack)
    return CommunitySystem(size_A, size_B,
                           span_control, pow_importance, will_attack,
                           reshape([BoundedVal(1.0, 0.0) for i=1:(size_A * size_B)], size_A, size_B), ##com_resource
                           reshape([BoundedVal(1.0, 0.0) for i=1:(size_A * size_B)], size_A, size_B), #com_power
                           collect(CartesianIndices((size_A, size_B))), ##com_superior (no hierarchy)
                           collect(CartesianIndices((size_A, size_B))), ##com_polity
                           zeros(Int64, size_A, size_B))
end

###Returns an iterator over all communities 2d indices (rows vary faster than columns)
all_comm(cs::CommunitySystem) = CartesianIndices((cs.grid_size_a, cs.grid_size_b))

function recover_resource!(cs::CommunitySystem, diff_decay)
    for c in eachindex(cs.com_resource)
            cs.com_resource[c].cur_lvl = cs.com_resource[c].max_lvl - 
                                         diff_decay * (cs.com_resource[c].max_lvl - cs.com_resource[c].cur_lvl)
    end
    return cs
end

##Updates the rank of each community
function update_rank!(cs::CommunitySystem)
    cs.com_rank .= 0
    for c in all_comm(cs)
        this_c = c
        sup_c = cs.com_superior[this_c]
        while this_c != sup_c
            ###Break if this branch was already considered with a higher rank
            cs.com_rank[sup_c] > cs.com_rank[this_c] && break
            ##Increase rank of superior
            cs.com_rank[sup_c] = cs.com_rank[this_c] + 1
            ##Go to superior
            this_c = cs.com_superior[this_c]
            sup_c = cs.com_superior[this_c]
        end
    end
    return cs
end


function update_power!(cs::CommunitySystem)
    ##Initialize power to 0
    set_cur!.(cs.com_power, 0.0)
    set_max!.(cs.com_power, 0.0)
    #####/!\ To check if there is non double counting
    for c in all_comm(cs)
        last_sup = c
        curr_sup = c
        hier_dist = 0 ##Distance in the hierarchy
        while curr_sup != last_sup
            ##Update current and maximum theoretical powers
            add_to_max!(cs.com_power[curr_sup], cs.tax_rate^hier_dist * get_max(cs.com_resource[c]))
            add_to_cur!(cs.com_power[curr_sup], cs.tax_rate^hier_dist * get_cur(cs.com_resource[c]))
            ##Break if we reached top of hierarchy
            cs.com_superior[curr_sup] == curr_sup && break
            last_sup = curr_sup
            curr_sup = cs.com_superior[curr_sup]
            hier_dist += 1
        end
    end
    return cs.com_power
end


###Payment of tribute
#=
function pay_tribute!(cs::CommunitySystem, tax_rate::Float64)
    ##Amount of tribute received by each community
    cs.com_res_received .= 0.0
    ###Iterate over communities and pay tribute to superior
    for c in all_comm(cs)
        if cs.com_superior[c...] != c
            last_sup = (0, 0)
            curr_sup = cs.com_superior[c...]
            hier_dist = 1 ##Distance in the hierarchy
            while curr_sup != last_sup
                cs.com_res_received[curr_sup...] += tax_rate^hier_dist * cs.com_resource[c...].cur_lvl * #part received
                                            (1.0 - (curr_sup != cs.com_superior[curr_sup...]) * tax_rate) #part given to superior
                last_sup = curr_sup
                curr_sup = cs.com_superior[curr_sup...]
                hier_dist += 1
            end
            cs.com_resource[c...].cur_lvl *= (1.0 - tax_rate)
        end
    end
    for c in all_comm(cs)
        cs.com_resource[c...].cur_lvl += cs.com_res_received[c...]
    end
    return cs
end
=#


### Space is topologically a torus, so it has no borders
###Convert any pair of coordinates so that they are on grid
on_grid(cs::CommunitySystem, t::NTuple{2, Int64}) = (t[1] % cs.grid_size_a) + 1, (t[2] % cs.grid_size_b) + 1

### Neighbors are connected points on a hexagonal grid
### We use the following coordinate system designate points: a * (1, 0) + b * (cos(pi/3), sin(pi/3)) 
### where (1, 0) and (0, 1) are the canonical bese vectors of R^2
### Space is topologically a torus, so it has no borders
function get_comm_neighbors(cs::CommunitySystem, com_index::ComCoord)
    a, b = com_index[1], com_index[2]
    (1 <= a <= cs.grid_size_a) & (1 <= b <= cs.grid_size_b) || error("Incorrect community coordinates.")
    return CartesianIndex.(
                           on_grid.(Ref(cs), ((a, b - 1), (a + 1, b - 1),
                                              (a - 1, b), (a + 1, b),
                                              (a - 1, b + 1), (a, b + 1))
                                    )
                           )
end


###Returns a dictionary with indices (2d indices) of polities as keys and a
###set of neighboring communities
function get_pol_neighbors(cs::CommunitySystem)
    ###Dictionary with polities as keys and set of communities as values
    pol_neigh_comm = Dict{ComCoord, Set{ComCoord}}(p => Set{ComCoord}() for p in all_comm(cs) if p == cs.com_polity[p])
    ###Go through communities and get thei neighbors
    ###if their neighbors are from other polities, add those neighbors to the set
    ###of neighboring communities of polities
    for c_i in all_comm(cs) ##First argument changes the fastest
        n_c = get_comm_neighbors(cs, c_i)
        for nc_i in n_c
            if (cs.com_polity[c_i] != cs.com_polity[nc_i]) & (nc_i ∉ pol_neigh_comm[cs.com_polity[c_i]])
                push!(pol_neigh_comm[cs.com_polity[c_i]], nc_i)
            end
        end
    end
    return pol_neigh_comm
end


function get_pol_neighbors_array(cs::CommunitySystem)
    ###Dictionary with polities as keys and set of communities as values
    pol_neigh_comm = Dict{ComCoord, Vector{ComCoord}}(p => ComCoord[] for p in all_comm(cs) if p == cs.com_polity[p])
    ###Go through communities and get their neighbors
    ###if their neighbors are from other polities, add those neighbors to the set
    ###of neighboring communities of polities
    for c_i in all_comm(cs) ##First argument changes the fastest
        n_c = get_comm_neighbors(cs, c_i)
        for nc_i in n_c
            (cs.com_polity[c_i] != cs.com_polity[nc_i]) && push!(pol_neigh_comm[cs.com_polity[c_i]], nc_i)
        end
    end
    ###Make communities in vector unique
    map!(sort!, values(pol_neigh_comm))
    map!(unique!, values(pol_neigh_comm))
    return pol_neigh_comm
end

###probability for an attacker or a rebel to win
prob_att_reb_win(cs::CommunitySystem, pow_att, pow_def) = pow_att^cs.power_importance / 
    (pow_att^cs.power_importance + pow_def^cs.power_importance)

###Probability of a rebel or an attacker to attack
function prob_attack_rebel(cs::CommunitySystem, pot_att_reb, target)
    ##potential attacker or rebel evaluates chances of winning
    p_pot_att_reb_win = prob_att_reb_win(cs, cs.com_power[pot_att_reb], 
                                         cs.com_power[target] - 
                                         ifelse(cs.com_polity[pot_att_reb] == cs.com_polity[target], 
                                                cs.tax_rate * get_cur(cs.com_power[pot_att_reb]), 0.0))
    ###/!| The ifelse part does not seem to be in Gavrilets code, but it is more consistent like this
    ##Probability of attacking is from the code provided by Sergey Gavrilets, note: this differs from what is in the article
    p_attack = p_pot_att_reb_win * get_cur(cs.com_power[pot_att_reb]) / get_max(cs.com_power[pot_att_reb]) *
               exp(- 2.0 * cs.will_to_attack * (1.0 - p_pot_att_reb_win))
    
    return p_attack
end


###Returns a dictionary with polities and their initial communities targets (if there is an attack)
###pol_niegh_comm is a dictionary with indices (2d indices) of polities as keys and a
###set of neighboring communities (i.e. belonging to other polities)
###pol_facing_reb is a disctionary of polities with an array of subchiefs rebelling against them
function get_pol_target(cs::CommunitySystem, pol_neigh_comm, pol_facing_reb)
    ###Dictionary with attacking polities as keys and initial target community
    pol_attck_target = Dict{ComCoord, ComCoord}()
    for att_p in keys(pol_neigh_comm)
        ##Only attack if no facing rebellion
        if !haskey(pol_facing_reb, att_p)
            ###Find weakest (least powerful) polity neighbor; initialize with 
            t_pol = ComCoord(0, 0) ##current target polity
            t_pol_pow = -1.0 ##current polity power
            t_pol_com = ComCoord(0, 0) ##current target community of the target polity
            ##### Counting number of polities or communities visited in case powers or resources are the same
            ##### in order to select one of them randomly (by doing reservoir sampling)
            n_wkst_visited = 0.0 ## Number of weakest polities the iterator went through
            n_wlthst_visited = 0.0 ## Number of wealthiest communities visited (belonging to same polity)
            for a_n_c in pol_neigh_comm[att_p]
                ###Second part of the conjunction: polities do not attack polities facing rebellions
                if (get_cur(cs.com_power[cs.com_polity[a_n_c]]) < t_pol_pow) & (!haskey(pol_facing_reb, cs.com_polity[a_n_c]))
                    t_pol = cs.com_polity[a_n_c]
                    t_pol_pow = cs.com_power[t_pol]
                    t_pol_com = a_n_c
                    n_wkst_visited = 1.0
                ##Case the polity is already a target but subordinate community is wealthier
                elseif t_pol == cs.com_polity[a_n_c]
                    if get_cur(cs.com_resource[a_n_c]) > get_cur(cs.com_resource[t_pol_com]) ##case considered communities has more resources
                        t_pol_com = a_n_c
                        n_wlthst_visited = 1.0
                    elseif get_cur(cs.com_resource[a_n_c]) == get_cur(cs.com_resource[t_pol_com])
                        n_wlthst_visited += 1.0
                        if rand() <= (1.0 / n_wlthst_visited)
                            t_pol_com = a_n_c
                        end
                    end
                ## case: powers are the same, choose randomly using reservoir sampling (reservoir size = 1)
                elseif get_cur(cs.com_power[cs.com_polity[a_n_c]]) == t_pol_pow
                    n_wkst_visited += 1.0
                    if rand() <= (1.0 / n_wkst_visited)
                        t_pol = cs.com_polity[a_n_c]
                        t_pol_com = a_n_c
                    end
                end
            end
            ##Should polity attack? if yes, add wealthiest community of weakest neighboring polity as target
            if t_pol != ComCoord(0, 0) && rand() <= prob_attack_rebel(cs, p, target)
                pol_attck_target[p] = t_pol_com
            end
        end
    end
    return pol_attck_target
end

####Find sub chiefs that are rebelling
function pol_rebel(cs::CommunitySystem)
    pol_rebel = Dict{ComCoord, Vector[ComCoord]}()
    for sub_c in all_comm(sc)
        ###Check if indeed direct sub chief
        if sc.com_superior[sub_c] == sc.com_polity[sub_c]
            ##sub_chief decides weather to attack or not
            p_reb = prob_attack_rebel(cs, sub_c, sc.com_superior[sub_c])
            ###If sub chief decides to rebel add it to list of rebels against its superior
            if rand() <= p_reb
                if haskey(pol_rebel, sc.com_superior[sub_c])
                    push!(pol_rebel[sc.com_superior[sub_c]], sub_c)
                else
                    pol_rebel[sc.com_superior[sub_c]] = ComCoord[sub_c]
                end
            end
        end
    end
    return pol_rebel
end

function rebel!(cs::CommunitySystem, pol_rebel)
    ####/!\ war cost array should be reinitialized
    ##Reinitialize rebellion success indicator
    sc.rebel_success .= false
    ##Check if rebellions are successfull and calculate their cost
    for pol in keys(pol_rebel)
        ##Get power of chief against rebels
        tot_reb_pow = sum(get_cur(sc.com_power[r]) for r in pol_rebel[reb])
        chief_def_pow = get_cur(cs.com_power[pol]) - sc.tax_rate * tot_reb_pow
        for reb in pol_rebel[reb]
            p_reb_succ = prob_att_reb_win(cs, get_cur(cs.com_power[reb]), 
                                              chief_def_pow * cs.com_power[reb] / tot_reb_pow)
            ##Check if rebellion is successful
            if rand() <= p_reb_succ
                ##Set rebellion to sucessful
                sc.secede_success[reb] = true
            end
            ##Update costs for chief and rebel
            sc.com_war_cost[reb] = cs.cost_factor * ifelse(sc.secede_success[reb], 1.0 - p_reb_succ, p_reb_succ)
            sc.com_war_cost[pol] += chief_def_pow * cs.com_power[reb] *
                                    cs.cost_factor * ifelse(sc.secede_success[reb], 1.0 - p_reb_succ, p_reb_succ)
        end
    end
    ##Calculate costs and update superiors and polities
    reb_vassals = sc.com_list_cache
    for c in all_comm(cs)
        curr_com = c
        n_reb_vassals = 0
        if haskey(pol_rebel, cs.com_polity[curr_com])
            while sc.com_war_cost[curr_com] == 0.0 ##Ends when rebel or chief is reached and does not enter in loop if vassal line already visited
                n_reb_vassals += 1
                reb_vassals[n_reb_vassals] = curr_com
                curr_com = cs.com_superior[curr_com]
            end
        end
        ##Update war cost
        cs.com_war_cost[view(reb_vassals, 1:n_reb_vassals)] .= cs.com_war_cost[curr_com]
        ##Update polity of rebel and vassals if rebellion is successful
        if cs.secede_success[curr_com]
            cs.superior[curr_com] = curr_com
            cs.com_polity[curr_com] = curr_com
            cs.com_polity[view(reb_vassals, 1:n_reb_vassals)] = curr_com
        end
    end
end


function attack!(cs, pol_attck_target)
    for att in keys(pol_attck_target)
        target = pol_attck_target[att]
        ###Calculate probability that attacker wins
        p_attck_win = prob_att_reb_win(cs, )

    end
end



##Updates resources given war or rebellion costs
function apply_war_reb_cost!(cs::CommunitySystem)
    cs.com_war_cost .= (1.0 .- cs.com_war_cost) 
    mult_cur!.(cs.com_resource, cs.com_war_cost)
end

