# I/O:
# function to read the simulation parameters (total number of filaments, params...)
# function to write out configurations
# function to write out ops (radius, energy, total overlap, etc)

# US:
# type for us info
# function to run simple adaptive US sim

module ActinRingMC

using Random

const kb = 1.380649e-23

struct SystemParams
    ks::Float64
    kd::Float64
    T::Float64
    delta::Float64
    Xc::Float64
    EI::Float64
    Lf::Float64
    lf::Float64
    Nfil::Float64
    Nsca::Float64
end

struct SimulationParams
    steps::Int
    write_interval::Int
    ring_radius_move_freq::Float64
end

mutable struct Filament
    lf::Float64
    coors::Array{Int,2}
    current_coors::Array{Int,2}
    trial_coors::Array{Int,2}
    index::Int
end

# make ops as fields?
mutable struct System
    filaments::Vector{Filament}
    filament_length::Int
    ring_r::Float64
    parms::SystemParams
end

mutable struct Lattice
    occupancy::Dict{Vector{Int},Tuple{Int,Int}}
    current_occupancy::Dict{Vector{Int},Tuple{Int,Int}}
    trial_occupancy::Dict{Vector{Int},Tuple{Int,Int}}
    height::Int
end

function use_current_coors!(system, lattice)
    lattice.occupancy = lattice.current_occupancy
    for filament in system
        filament.coors = filament.current_coors
    end
end

function use_trial_coors!(system, lattice)
    lattice.occupancy = lattice.trial_occupancy
    for filament in system
        filament.coors = filament.trial_coors
    end
end

function wrap_pos!(lattice, pos::Vector{Int})
    if pos[2] > lattice.height
        pos[2] = -lattice.height + pos[2]
    elseif pos[2] < 0
        pos[2] = lattice.height + pos[2]
    end

    return nothing
end

function update_radius(system, lattice, height)
    lattice.height = height
    system.ring_r = system.parms.delta*(height - 1)
end

# only works with even number of sca filaments
function generate_starting_config(lattice, Nfil, Nsca, lf)
    filaments = []
    pos = [0, 0]
    for N_i in 1:Nfil
        for _ in 1:div(Nsca, 2)
            coors = zeros(2, lf)
            for j in 1:lf
                coors[:, j] = pos
                pos[2] += 1
                wrap_pos!(lattice, pos)
            end
            filament = Filament(lf, coors, coors, deepcopy(coors), N_i)
            filaments.push(filament)
            pos[2] += lf - 1
            wrap_pos!(lattice, pos)
        end
        pos[1] += 1
        pos[1] % 2 ? pos[2] = lf - 1: pos[2] = 0
        wrap_pos!(lattice, pos)
    end
end

function select_rand_filament(system)
    rand(system.filaments)
end

function generate_random_vector()
    [rand(0:5), rand(0:5)]
end

function calc_total_energy(system, lattice)
    ene = 0
    for filament in system.filaments
        ene += calc_filament_overlap_energy(system, lattice, filament)/2
        ene += calc_filament_bending_energy(system)
    end
    return ene
end

function calc_filament_bending_energy(system)
    system.parms.EI * system.parms.Lf / (2 * system.ring_r^2)
end

function calc_overlap_energy(system, l)
    L = system.delta * (l - 1)
    factor_t = L * kb * system.parms.T / system.parms.delta 
    ks = system.parms.ks
    kd = system.parms.kd
    Xc = system.parms.Xc
    log_t = log(1 + ks^2 * Xc / (kd * (ks + Xc)^2))
    factor_t * log_t
end

function calc_filament_overlap_energy(system, lattice, filament)
    l = 0
    for i in 1:system.lf
        pos = filament.coors[:, i]
        if pos in lattice.occupancy.keys
            return Inf
        else
            for dx in [-1, 1]
                adj_pos = [pos[1] + dx, pos[2]]
                if adj_pos in lattice.occupancy.keys
                    l += 1
                end
            end
        end
    end
    calc_overlap_energy(system, l)
end

function calc_filament_energy_diff!(system, lattice, filament)
    use_current_coors!(system, lattice)
    cur_overlap_ene = calc_filament_overlap_energy(system, lattice, filament)
    cur_bending_ene = calc_filament_bending_energy(system)

    use_trial_coors!(system, lattice)
    trial_overlap_ene = calc_filament_overlap_energy(system, lattice, filament)
    trial_bending_ene = calc_filament_bending_energy(system)

    trial_overlap_ene + trial_bending_ene - cur_overlap_ene - cur_bending_ene
end

function calc_system_energy_diff!(system, lattice)
    use_current_coors!(system, lattice)
    current_ene = calc_total_energy(system, lattice)

    use_trial_coors!(system, lattice)
    trial_ene = calc_total_energy(system, lattice)

    trial_ene - current_ene
end

function ring_and_system_connected(system, filament, lattice)
    filaments_in_path = [filament.index]
    pos = deepcopy(filament.trial_coors[:, 0])
    site_i = 0
    path_length = 0
    while site_i > filament.lf
        for dx in [-1, 1]
            adj_pos = [pos[1] + dx, pos[2]]
            if adj_pos in lattice.occupancy.keys
                filament_i, site_i = lattice.occupancy[adj_pos]
                filaments_in_path.push!(filament.index)
                next_filament = system.filaments[filament_i]
                if search_filament_for_path(system, lattice, next_filament,
                                            site_i, filaments_in_path,
                                            path_length, false)
                    return true
                end
            end
        end
        site_i + 1
        pos += [0, 1]
        wrap_pos!(lattice, pos)
        path_length += 1
    end
    return false
end

function search_filament_for_path(system, lattice, filament, site_i,
        filaments_in_path, path_length, ring_contig)

    initial_pos = deepcopy(filament[:, site_i])
    initial_path_length = path_length
    pos = deepcopy(initial_pos)
    dir = -1
    initial_site_i = site_i
    while true
        if site_i < 1 || site_i > filament.lf
            if dir == -1
                dir = 1
                path_length = initial_path_length + dir
                site_i = initial_site_i + dir
                pos = initial_pos + [0, dir]
                wrap_pos!(lattice, pos)
                continue
            else
                break
            end
        end
        for dx in [-1, 1]
            adj_pos = [pos[1] + dx, pos[2]]
            if adj_pos in lattice.occupancy.keys
                next_filament_i, next_site_i = lattice.occupancy[adj_pos]
                if (next_filament_i == filaments_in_path[1] &&
                     path_length > filament.lf)

                    ring_contig = true
                end
                if ring_contig && length(filaments_in_path) == system.parms.Nfil
                    return true
                elseif !(next_filament_i in filaments_in_path)
                    filaments_in_path.push!(filament.index)
                    next_filament = system.filaments[next_filament_i]
                    if search_filament_for_path(system, lattice,
                                                next_filament,
                                                next_site_i, filaments_in_path,
                                                path_length, ring_contig)
                        return true
                    end
                    filaments_in_path.pop!()
                end
            end
        end
        site_i + dir
        pos += [0, dir]
        wrap_pos!(lattice, pos)
        path_length += dir
    end
    return false
end

function translate_filament!(filament, lattice, move_vector)
    for i in 1:filament.length
        pos = filament.coors[:, i]
        delete!(lattice.occupancy, pos)
        pos += move_vector
        wrap_pos!(lattice, pos)
        lattice.occupancy[deepcopy(pos)] = filament.index
    end
end

function filament_accept_trial!(filament, lattice)
    filament.current_coors = deepcopy(filament.trial_coors)
    lattice.current_occupancy = deepcopy(lattice.trial_occupancy)
end

function filament_accept_current!(filament, lattice)
    filament.trial_coors = deepcopy(filament.current_coors)
    lattice.trial_occupancy = deepcopy(lattice.current_occupancy)
end

function system_accept_trial!(system, lattice)
    for filament in system.filaments
        filament.current_coors = deepcopy(filament.trial_coors)
    end
    lattice.current_occupancy = deepcopy(lattice.trial_occupancy)
end

function system_accept_current!(system, lattice)
    for filament in system.filaments
        filament.trial_coors = deepcopy(filament.current_coors)
    end
    lattice.trial_occupancy = deepcopy(lattice.current_occupancy)
end

function select_move(simparms)
    if simparms.ring_radius_move_freq > rand(Float64)
        return attempt_ring_radius_move!
    else
        return attempt_translation_move!
    end
end

function accept_move(delta_energy)
    p_accept = min(1, exp(delta_energy))
    accept = false
    if p_accept == 1 || p_accept > rand(Float64)
        accept = true;
    end

    return accept;
end

function attempt_translation_move!(system, lattice)
    use_trial_coors!(system, lattice)
    filament = select_rand_filament(system)
    move_vector = generate_random_vector()
    translate_filament!(filament, lattice, move_vector)
    # Should I be accepting the current or just propose again?
    if !ring_and_system_connected(system, filament, lattice)
        filament_accept_current!(filament, lattice)
        return nothing
    end
    delta_energy = calc_filament_energy_diff!(system, lattice, filament)
    if accept_move(delta_energy)
        filament_accept_trial!(filament, lattice)
    else
        filament_accept_current!(filament, lattice)
    end

    return nothing
end

function find_filaments_on_boundary(system, lattice)
    filaments_on_boundary = []
    for filament in system.filaments
        for i in 1:filament.lf
            if filament.coors[2, i] == lattice.height
                filaments_on_boundary.push((filament, i))
                break
            end
        end
    end
end

function translate_filaments_on_boundary!(filaments_on_boundary, dir)
    for (filament, i) in filaments_on_boundary
        for j in 1:i
            filament[2, j] += dir
        end
    end
end

function attempt_ring_radius_move!(system, lattice)
    dir = rand([-1, 1])
    use_trial_coors!(system, lattice)
    filaments_on_boundary = find_filaments_on_boundary(system, lattice)
    translate_filaments_on_boundary!(filaments_on_boundary, dir)
    update_radius(system, lattice, lattice.height + dir)
    if dir == 1
        # Should I be accepting the current or just propose again?
        while !ring_and_system_connected(system, system.filaments[0], lattice)
            system_accept_current!(system, lattice)
            return nothing
        end
    end
    delta_energy = calc_system_energy_diff!(system, lattice)
    if accept_move(delta_energy)
        system_accept_trial!(system, lattice)
    else
        system_accept_current!(system, lattice)
    end
end

function write_ops()
end

# filebase or already setup files
function run(system::System, lattice::Lattice, filebase, simparms)
    for _ in 1:simparms.steps
        attempt_move! = select_move(simparms)
        attempt_move!(system, lattice)
        if step % simparms.write_interval
            write_ops()
        end
    end
end

end
