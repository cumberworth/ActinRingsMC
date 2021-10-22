using Test

import ActinRingsMC as armc

ks = 7.82e-9
kd = 7.82e-9/5
T = 300
delta = 36e-9
Xc = 12e-9
EI = 3.6e-26
Lf = 3e-6
Nfil = 6
Nsca = 4

#lf = armc.calc_lf(Lf, delta)
lf = 84
Lf = armc.calc_Lf(lf, delta)

# Test energy functions
max_height = armc.calc_max_lattice_height(Nsca, lf)
min_height = div(max_height, 2)
overlap = 21
start_height = max_height - Nsca*overlap
start_radius = armc.calc_radius(delta, start_height)

lattice = armc.Lattice(start_height, max_height, min_height)
filaments = armc.generate_starting_config(lattice, Nfil, Nsca, lf, overlap)
armc.update_lattice_occupancies!(filaments, lattice)
sysparms = armc.SystemParams(ks, kd, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
system = armc.System(sysparms, filaments, start_radius)

@test isapprox(armc.calc_overlap_energy(system, lattice), -5.729146253430745e-19)
@test isapprox(armc.calc_total_bending_energy(system), 1.617547027233688e-19)
@test isapprox(armc.calc_total_energy(system, lattice), -4.111599226197057e-19)

#file = armc.prepare_vtf_file("test.vtf", system)
#armc.write_vtf(system, file)
#close(file)

#iters = 10
#steps = 1e5
#max_bias_diff = 10
#write_interval = 1e3
#radius_move_freq = 0.5
#filebase = "test"

#simparms = armc.SimulationParams(iters, steps, max_bias_diff, write_interval,
#                                 radius_move_freq, filebase)

#armc.run!(system, lattice, simparms)
#armc.run_us!(system, lattice, simparms)
