module ActinRingsMC

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
    lf::Int
    Nfil::Int
    Nsca::Int
end

struct SimulationParams
    iters::Int
    steps::Int
    max_bias_diff::Float64
    write_interval::Int
    radius_move_freq::Float64
    filebase::String
end

mutable struct Filament
    lf::Int
    coors::Array{Int,2}
    current_coors::Array{Int,2}
    trial_coors::Array{Int,2}
    index::Int
end

mutable struct System
    parms::SystemParams
    filaments::Vector{Filament}
    radius::Float64
    energy::Float64
end

function System(parms::SystemParams, filaments::Vector{Filament}, radius::Float64)
    System(parms, filaments, radius, 0)
end

mutable struct Lattice
    occupancy::Dict{Vector{Int},Tuple{Int,Int}}
    current_occupancy::Dict{Vector{Int},Tuple{Int,Int}}
    trial_occupancy::Dict{Vector{Int},Tuple{Int,Int}}
    height::Int
    current_height::Int
    trial_height::Int
    max_height::Int
    min_height::Int
end

function Lattice(height::Int, max_height, min_height)
    occupancy::Dict{Vector{Int},Tuple{Int,Int}} = Dict()
    Lattice(occupancy, occupancy, copy(occupancy), height, height, height,
            max_height, min_height)
end

function calc_lf(Lf::Float64, delta::Float64)::Int
    div(Lf, delta)
end

function calc_Lf(lf::Int, delta::Float64)
    delta * lf
end

function calc_max_lattice_height(Nsca::Int, lf::Int)
    Nsca * lf - Nsca - 1
end

function calc_min_lattice_height(Nsca::Int, lf::Int)
    div(Nsca, 2)*lf - 1
end

function calc_radius(delta::Float64, height::Int)
    delta * (height + 1)/(2pi)
end

function update_radius!(system::System, lattice::Lattice, height::Int)
    lattice.height = height
    system.radius = calc_radius(system.parms.delta, height)
    return nothing
end

function use_current_coors!(system::System, lattice::Lattice)
    lattice.occupancy = lattice.current_occupancy
    update_radius!(system, lattice, lattice.current_height)
    for filament in system.filaments
        filament.coors = filament.current_coors
    end
    return nothing
end

function use_trial_coors!(system::System, lattice::Lattice)
    lattice.occupancy = lattice.trial_occupancy
    update_radius!(system, lattice, lattice.trial_height)
    for filament in system.filaments
        filament.coors = filament.trial_coors
    end
    return nothing
end

function accept_trial!(filament::Filament, lattice::Lattice)
    filament.current_coors = copy(filament.trial_coors)
    lattice.current_occupancy = copy(lattice.trial_occupancy)
end

function accept_trial!(system::System, lattice::Lattice)
    for filament in system.filaments
        filament.current_coors = copy(filament.trial_coors)
    end
    lattice.current_occupancy = copy(lattice.trial_occupancy)
    lattice.current_height = lattice.trial_height
end

function accept_current!(filament::Filament, lattice::Lattice)
    filament.trial_coors = copy(filament.current_coors)
    lattice.trial_occupancy = copy(lattice.current_occupancy)
end

function accept_current!(system::System, lattice::Lattice)
    for filament in system.filaments
        filament.trial_coors = copy(filament.current_coors)
    end
    lattice.trial_occupancy = copy(lattice.current_occupancy)
    lattice.trial_height = lattice.current_height
end

function wrap_pos!(lattice::Lattice, pos::Vector{Int})
    if pos[2] > lattice.height
        pos[2] = -(lattice.height + 1) + pos[2]
    elseif pos[2] < 0
        pos[2] = (lattice.height + 1) + pos[2]
    end
    return nothing
end

# only works with even number of sca filaments
function generate_starting_config(lattice::Lattice, Nfil::Int, Nsca::Int,
        lf::Int, overlap::Int)
    filaments::Vector{Filament} = []
    pos = [0, 0]
    Ni = 1
    while Ni <= Nfil
        for _ in 1:div(Nsca, 2)
            coors = zeros(2, lf)
            for j in 1:lf
                coors[:, j] = pos
                pos[2] += 1
                wrap_pos!(lattice, pos)
            end
            filament = Filament(lf, coors, coors, copy(coors), Ni)
            push!(filaments, filament)
            pos[2] += lf - 2*overlap
            wrap_pos!(lattice, pos)
            Ni += 1
            if Ni > Nfil
                break
            end
        end
        pos[1] += 1
        pos[1] % 2 == 0 ? pos[2] = 0 : pos[2] = lf - overlap
        wrap_pos!(lattice, pos)
    end
    return filaments
end

function update_lattice_occupancies!(filaments::Vector{Filament}, lattice::Lattice)
    lattice.current_occupancy = Dict()
    lattice.trial_occupancy = Dict()
    for filament in filaments
        for i in 1:filament.lf
            pos = filament.coors[:, i]
            lattice.current_occupancy[copy(pos)] = (filament.index, i)
            lattice.trial_occupancy[copy(pos)] = (filament.index, i)
        end
    end
    lattice.occupancy = lattice.current_occupancy
    return nothing
end

function select_rand_filament(system::System)
    rand(system.filaments)
end

function generate_random_vector()
    [0, rand([-1, 1])]
end

function calc_filament_bending_energy(system::System)
    system.parms.EI * system.parms.Lf / (2 * system.radius^2)
end

function calc_total_bending_energy(system::System)
    ene = 0
    for _ in 1:system.parms.Nfil
        ene += calc_filament_bending_energy(system)
    end
    return ene
end

function calc_overlap_energy(system::System, L::Float64)
    factor_t = L * kb * system.parms.T / system.parms.delta 
    ks = system.parms.ks
    kd = system.parms.kd
    Xc = system.parms.Xc
    log_t = log(1 + ks^2 * Xc / (kd * (ks + Xc)^2))
    -factor_t * log_t
end

function calc_overlap_energy(system::System, lattice::Lattice, filament::Filament)
    l = 0
    for i in 1:system.parms.lf
        pos = filament.coors[:, i]
        for dx in [-1, 1]
            adj_pos = [pos[1] + dx, pos[2]]
            if adj_pos in keys(lattice.occupancy)
                l += 1
            end
        end
    end
    L = system.parms.delta * l
    calc_overlap_energy(system, L)
end

function calc_overlap_energy(system::System, lattice::Lattice)
    ene = 0
    for filament in system.filaments
        ene += calc_overlap_energy(system, lattice, filament)/2
    end
    return ene
end

function calc_bias_energy(lattice::Lattice, biases::Vector{Float64})
    biases[lattice.height - lattice.min_height + 1]
end

function calc_total_energy(system::System, lattice::Lattice)
    ene = 0
    for filament in system.filaments
        ene += calc_overlap_energy(system, lattice, filament)/2
        ene += calc_filament_bending_energy(system)
    end
    return ene
end

function calc_total_energy(system::System, lattice::Lattice, biases::Vector{Float64})
    ene = calc_total_energy(system, lattice)
    ene += calc_bias_energy(system, biases)
    return ene
end

function calc_energy_diff!(system::System, lattice::Lattice, filament::Filament)
    use_current_coors!(system, lattice)
    cur_overlap_ene = calc_overlap_energy(system, lattice, filament)
    cur_bending_ene = calc_filament_bending_energy(system)

    use_trial_coors!(system, lattice)
    trial_overlap_ene = calc_overlap_energy(system, lattice, filament)
    trial_bending_ene = calc_filament_bending_energy(system)

    trial_overlap_ene + trial_bending_ene - cur_overlap_ene - cur_bending_ene
end

function calc_energy_diff!(system::System, lattice::Lattice, biases::Vector{Float64})
    use_current_coors!(system, lattice)
    current_ene = calc_total_energy(system, lattice)
    current_ene += calc_bias_energy(lattice, biases)

    use_trial_coors!(system, lattice)
    trial_ene = calc_total_energy(system, lattice)
    trial_ene += calc_bias_energy(lattice, biases)

    trial_ene - current_ene
end

function ring_and_system_connected(system::System, lattice::Lattice, filament::Filament)
    pos = filament.coors[:, 1]
    site_i = 1
    path::Vector{Int} = [filament.index]
    path_length = 0
    ring_contig = false
    connected_filaments::Set{Int} = Set(filament.index)
    while site_i <= filament.lf
        for dx in [-1, 1]
            adj_pos = [pos[1] + dx, pos[2]]
            if adj_pos in keys(lattice.occupancy)
                adj_filament_i, adj_site_i = lattice.occupancy[adj_pos]
                if adj_filament_i in connected_filaments
                    continue
                end
                union!(connected_filaments, adj_filament_i)
                adj_filament = system.filaments[adj_filament_i]
                ring_contig, connected_filaments = search_filament_for_path(
                    system, lattice, adj_filament, adj_site_i, path,
                    path_length, ring_contig, connected_filaments)
                if ring_contig && length(connected_filaments) == system.parms.Nfil
                    return true
                end
            end
        end
        site_i += 1
        pos += [0, 1]
        wrap_pos!(lattice, pos)
        path_length += 1
    end
    return false
end

# this seems to work but its really long and has too many arguments
function search_filament_for_path(system::System, lattice::Lattice,
        filament::Filament, site_i::Int, path::Vector{Int},
        path_length::Int, ring_contig::Bool, connected_filaments::Set{Int})

    searched_filaments::Set{Int} = Set()
    initial_pos = filament.coors[:, site_i]
    initial_path_length = path_length
    pos = copy(initial_pos)
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
            if adj_pos in keys(lattice.occupancy)
                next_filament_i, next_site_i = lattice.occupancy[adj_pos]
                if (next_filament_i == path[1] &&
                    abs(path_length) > filament.lf)
                    ring_contig = true
                end
                if ring_contig && length(connected_filaments) == system.parms.Nfil
                    return (ring_contig, connected_filaments) 
                end
                if (next_filament_i in searched_filaments ||
                    next_filament_i in path)
                    continue
                end
                union!(connected_filaments, next_filament_i)
                union!(searched_filaments, next_filament_i)
                push!(path, next_filament_i)
                next_filament = system.filaments[next_filament_i]
                ring_contig, connected_filaments = search_filament_for_path(
                    system, lattice, next_filament, next_site_i,
                    path, path_length, ring_contig, connected_filaments)
                if ring_contig && length(connected_filaments) == system.parms.Nfil
                    return (ring_contig, connected_filaments)
                end
                pop!(path)
            end
        end
        site_i += dir
        pos += [0, dir]
        wrap_pos!(lattice, pos)
        path_length += dir
    end
    return (ring_contig, connected_filaments)
end

function filaments_contiguous(system::System, lattice::Lattice)
    for filament in system.filaments
        prev_pos = filament.coors[:, 1]
        for i in 2:filament.lf
            pos = filament.coors[:, i]
            exp_pos = [prev_pos[1], prev_pos[2] + 1]
            wrap_pos!(lattice, exp_pos)
            if exp_pos != pos
                return false
            end
            prev_pos = pos
        end
    end
    return true
end

function accept_move(system::System, delta_energy::Float64)
    p_accept = min(1, exp(-delta_energy/kb/system.parms.T))
    accept = false
    if p_accept == 1 || p_accept > rand(Float64)
        accept = true;
    end

    return accept;
end

function translate_filament!(filament::Filament, lattice::Lattice, move_vector::Vector{Int})
    for i in 1:filament.lf
        pos = filament.coors[:, i]
        delete!(lattice.occupancy, pos)
    end
    for i in 1:filament.lf
        pos = filament.coors[:, i]
        pos += move_vector
        wrap_pos!(lattice, pos)
        if pos in keys(lattice.occupancy)
            return false
        end
        filament.coors[:, i] = pos
        lattice.occupancy[copy(pos)] = (filament.index, i)
    end
    return true
end

# Should I be accepting the current or just propose again?
function attempt_translation_move!(system::System, lattice::Lattice, biases::Vector{Float64})
    use_trial_coors!(system, lattice)
    filament = select_rand_filament(system)
    move_vector = generate_random_vector()
    if !translate_filament!(filament, lattice, move_vector)
        accept_current!(system, lattice)
        use_current_coors!(system, lattice)
        return false
    end
    if !ring_and_system_connected(system, lattice, filament)
        accept_current!(filament, lattice)
        use_current_coors!(system, lattice)
        return false
    end
    delta_energy = calc_energy_diff!(system, lattice, filament)
    accept = accept_move(system, delta_energy)
    if accept
        accept_trial!(filament, lattice)
    else
        accept_current!(filament, lattice)
    end
    use_current_coors!(system, lattice)
    return accept
end

function find_split_points(filaments::Vector{Filament}, lattice::Lattice, dir::Int)
    split_points::Vector{Int} = []
    for filament in filaments
        push!(split_points, 0)
        if filament.coors[2, 1] == 0 && dir == -1
            break
        end
        for i in 1:filament.lf
            split_points[end] += 1
            if filament.coors[2, i] == lattice.height
                break
            end
        end
    end
    return split_points
end

function translate_filaments_with_split_points!(filaments::Vector{Filament},
    split_points::Vector{Int}, lattice::Lattice, dir::Int)
    for (filament, split_point) in zip(filaments, split_points)
        for i in 1:split_point
            pos = filament.coors[:, i]
            delete!(lattice.occupancy, pos)
        end
    end

    # Shift segments of filaments below split point only
    for (filament, split_point) in zip(filaments, split_points)
        for i in 1:split_point
            filament.coors[2, i] += dir
            pos = filament.coors[:, i]
            if dir == -1
                wrap_pos!(lattice, pos)
            end
            if !(pos in keys(lattice.occupancy))
                lattice.occupancy[pos] = (filament.index, i)
            else
                return false
            end
        end
    end
    return true
end

function attempt_radius_move!(system::System, lattice::Lattice, biases::Vector{Float64})
    dir = rand([-1, 1])
    if dir == 1 && lattice.height == lattice.max_height
        dir = -1
    elseif dir == -1 && lattice.height == lattice.min_height
        dir = 1
    end
    use_trial_coors!(system, lattice)
    filaments = rand(system.filaments, rand(1:system.parms.Nfil))
    split_points = find_split_points(filaments, lattice, dir)
    if !translate_filaments_with_split_points!(filaments, split_points, lattice, dir)
        accept_current!(system, lattice)
        use_current_coors!(system, lattice)
        return false
    end
    lattice.trial_height += dir
    update_radius!(system, lattice, lattice.trial_height)
    if !filaments_contiguous(system, lattice)
        accept_current!(system, lattice)
        use_current_coors!(system, lattice)
        return false
    end
    if dir == 1
        # Should I be accepting the current or just propose again?
        if !ring_and_system_connected(system, lattice, system.filaments[1])
            accept_current!(system, lattice)
            use_current_coors!(system, lattice)
            return false
        end
    end
    delta_energy = calc_energy_diff!(system, lattice, biases)
    accept = accept_move(system, delta_energy)
    if accept
        accept_trial!(system, lattice)
        use_current_coors!(system, lattice)
    else
        update_radius!(system, lattice, lattice.height - dir)
        accept_current!(system, lattice)
        use_current_coors!(system, lattice)
    end
    return accept
end

function select_move(simparms::SimulationParams)
    if simparms.radius_move_freq > rand(Float64)
        return attempt_radius_move!
    else
        return attempt_translation_move!
    end
end

function prepare_ops_file(filename::String)
    file = open(filename, "w")
    header = "step energy height radius"
    println(file, header)
    return file
end

function prepare_us_file(filename::String, lattice::Lattice)
    file = open(filename, "w")
    header = ""
    for h in lattice.min_height:lattice.max_height
        header = "$(header)$h "
    end
    println(file, header)
    return file
end

function prepare_vtf_file(filename::String, system::System)
    file = open(filename, "w")
    atom_i = 0
    for filament in system.filaments
        end_atom_i = atom_i + filament.lf - 1
        println(file, "a $atom_i:$end_atom_i c $(filament.index) r 2.5")
        atom_i += filament.lf
    end
    println(file, "")
    return file
end

function write_vtf(system::System, file::IOStream)
    println(file, "t")
    for filament in system.filaments
        for i in 1:filament.lf
            # Widen the aspect ratio for easier viewing
            x = filament.coors[1, i]*10
            y = filament.coors[2, i]
            println(file, "$x $y 0")
        end
    end
    println(file, "")
end

function write_ops(system::System, lattice::Lattice, step::Int, file::IOStream)
    println(file, "$step $(system.energy) $(lattice.height) $(system.radius)")
end

function write_us_data(data::Vector, file::IOStream)
    for d in data
        print(file, "$d ")
    end
    println(file, "")
end

function update_counts(counts::Vector{Int}, lattice)
    counts[lattice.height - lattice.min_height + 1] += 1
    return nothing
end

# filebase or already setup files
function run!(system::System, lattice::Lattice, simparms::SimulationParams,
        counts::Vector{Int}, biases::Vector{Float64}, ops_file::IOStream,
        vtf_file::IOStream)
    for step in 1:simparms.steps
        attempt_move! = select_move(simparms)
        attempt_move!(system, lattice, biases)
        update_counts(counts, lattice)
        if step % simparms.write_interval == 0
            println(step)
            system.energy = calc_total_energy(system, lattice)
            write_ops(system, lattice, step, ops_file)
            write_vtf(system, vtf_file)
        end
    end
end

function run!(system::System, lattice::Lattice, simparms::SimulationParams)
    ops_file = prepare_ops_file("$(simparms.filebase).ops")
    vtf_file = prepare_vtf_file("$(simparms.filebase).vtf", system)
    counts::Vector{Int} = zeros(lattice.max_height - lattice.min_height + 1)
    biases::Vector{Float64} = zeros(lattice.max_height - lattice.min_height + 1)
    run!(system, lattice, simparms, counts, biases, ops_file, vtf_file)
    close(ops_file)
    close(vtf_file)
end

function update_biases!(counts::Vector{Int}, freqs::Vector{Float64},
        probs::Vector{Float64}, biases::Vector{Float64}, T::Float64,
        max_bias_diff::Float64)
    norm = sum(counts .* exp.(biases))
    max_bias_diff *= kb*T
    for i in 1:length(counts)
        if counts[i] == 0
            freqs[i] = 0
            probs[i] = 0
            bias_diff = -max_bias_diff
        else
            freqs[i] = 1 / counts[i]
            probs[i] = counts[i] / norm
            bias_diff = kb*T*log.(probs[i]) - biases[i]
            if bias_diff > max_bias_diff
                bias_diff = max_bias_diff
            elseif bias_diff < -max_bias_diff
                bias_diff = -max_bias_diff
            end
        end
        biases[i] += bias_diff
        counts[i] = 0
    end
end

function run_us!(system::System, lattice::Lattice, simparms::SimulationParams)
    counts::Vector{Int} = zeros(lattice.max_height - lattice.min_height + 1)
    freqs::Vector{Float64} = zeros(lattice.max_height - lattice.min_height + 1)
    probs::Vector{Float64} = zeros(lattice.max_height - lattice.min_height + 1)
    biases::Vector{Float64} = zeros(lattice.max_height - lattice.min_height + 1)
    freqs_file = prepare_us_file("$(simparms.filebase).freqs", lattice)
    biases_file = prepare_us_file("$(simparms.filebase).biases", lattice)
    for i in 1:simparms.iters
        println("Iter $i")
        iter_filebase = "$(simparms.filebase)_iter-$i"
        ops_file = prepare_ops_file("$iter_filebase.ops")
        vtf_file = prepare_vtf_file("$iter_filebase.vtf", system)
        run!(system, lattice, simparms, counts, biases, ops_file, vtf_file)
        update_biases!(counts, freqs, probs, biases, system.parms.T,
                       simparms.max_bias_diff)
        write_us_data(freqs, freqs_file)
        write_us_data(biases, biases_file)
        close(ops_file)
        close(vtf_file)
    end
    close(freqs_file)
    close(biases_file)
end

end
