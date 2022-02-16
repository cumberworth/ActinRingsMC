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
Nfil = 6
Nsca = 4

Lf = armc.calc_Lf(lf, delta)

max_height = armc.calc_max_lattice_height(Nsca, lf)
min_height = armc.calc_min_lattice_height(Nsca, lf)
overlap = 1
start_height = max_height - Nsca*(overlap - 1)
start_radius = armc.calc_radius(delta, start_height)

lattice = armc.Lattice(start_height, max_height, min_height)
filaments = armc.generate_starting_config(lattice, Nfil, Nsca, lf, overlap)
armc.update_occupancies!(filaments, lattice)
sysparms = armc.SystemParams(ks, kd, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
system = armc.System(sysparms, filaments, start_radius)

mkpath("outs")
file = armc.prepare_vtf_file("outs/test.vtf", system)
armc.write_vtf(system, file)
close(file)

iters = 2
steps = 1e4
max_bias_diff = 5
write_interval = 1e3
radius_move_freq = 0.5
filebase = "outs/test"
analytical_biases = true
read_biases = false
biases_filename = ""
restart_iter = 0
binwidth = 1

simparms = armc.SimulationParams(
    iters,
    steps,
    max_bias_diff,
    write_interval,
    radius_move_freq,
    filebase,
    analytical_biases,
    read_biases,
    biases_filename,
    restart_iter,
    binwidth
)

armc.run!(system, lattice, simparms)
armc.run_us!(system, lattice, simparms)

filebase = "outs/test_run-1"
read_biases = true
biases_filename = "outs/test.biases"
restart_iter = 2
simparms = armc.SimulationParams(
    iters,
    steps,
    max_bias_diff,
    write_interval,
    radius_move_freq,
    filebase,
    analytical_biases,
    read_biases,
    biases_filename,
    restart_iter,
    binwidth
)

armc.run_us!(system, lattice, simparms)

main()
rm("outs", recursive=true)
