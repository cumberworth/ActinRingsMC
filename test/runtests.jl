using Test

import ActinRingsMC as armc

function main()
    energy_test()
end

function energy_test()
    ks = 7.82e-9
    kd = 7.82e-9/5
    T = 300
    delta = 36e-9
    Xc = 12e-9
    EI = 3.6e-26
    lf = 84
    Nfil = 6
    Nsca = 4

    Lf = armc.calc_Lf(lf, delta)

    # Test energy functions
    max_height = armc.calc_max_lattice_height(Nsca, lf)
    min_height = armc.calc_min_lattice_height(Nsca, lf)
    overlap = 21
    start_height = max_height - Nsca*overlap
    start_radius = armc.calc_radius(delta, start_height)

    lattice = armc.Lattice(start_height, max_height, min_height)
    filaments = armc.generate_starting_config(lattice, Nfil, Nsca, lf, overlap)
    armc.update_occupancies!(filaments, lattice)
    sysparms = armc.SystemParams(ks, kd, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
    system = armc.System(sysparms, filaments, start_radius)

    @test isapprox(armc.overlap_energy(system, lattice), -5.729146253430745e-19)
    @test isapprox(armc.total_bending_energy(system), 1.617547027233688e-19)
    @test isapprox(armc.total_energy(system, lattice), -4.111599226197057e-19)
end

ks = 7.82e-9
kd = 7.82e-9/5
T = 300
delta = 36e-9
Xc = 12e-9
EI = 3.6e-26
lf = 84
Nfil = 8
Nsca = 2

Lf = armc.calc_Lf(lf, delta)

max_height = armc.calc_max_lattice_height(Nsca, lf)
min_height = armc.calc_min_lattice_height(Nsca, lf)
overlap = 21
start_height = max_height - Nsca*overlap
start_radius = armc.calc_radius(delta, start_height)

lattice = armc.Lattice(start_height, max_height, min_height)
filaments = armc.generate_starting_config(lattice, Nfil, Nsca, lf, overlap)
armc.update_occupancies!(filaments, lattice)
sysparms = armc.SystemParams(ks, kd, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
system = armc.System(sysparms, filaments, start_radius)

file = armc.prepare_vtf_file("test.vtf", system)
armc.write_vtf(system, file)
close(file)

iters = 10
steps = 1e4
max_bias_diff = 5
write_interval = 1e3
radius_move_freq = 0.5
filebase = "test"
analytical_biases = true

simparms = armc.SimulationParams(
    iters,
    steps,
    max_bias_diff,
    write_interval,
    radius_move_freq,
    filebase,
    analytical_biases
)

armc.run!(system, lattice, simparms)
#armc.run_us!(system, lattice, simparms)

main()
