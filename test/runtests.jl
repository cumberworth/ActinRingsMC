using Test

import ActinRingsMC as armc

kd = 7.82e-9
ks = 7.82e-9/5
T = 300
delta = 36e-9
Xc = 12e-9
EI = 3.6e-26
Lf = 3e-6
Nfil = 6
Nsca = 4

iters = 10
steps = 1e4
max_bias_diff = 10
write_interval = 1e2
radius_move_freq = 1.0
filebase = "test"

lf = armc.calc_lf(Lf, delta)
Lf = armc.calc_Lf(lf, delta)
lattice_height = armc.calc_lattice_height(Nsca, lf)
max_radius = Lf*Nsca
max_height = lf*Nsca
min_height = div(max_height, 2)

lattice = armc.Lattice(lattice_height, max_height, min_height)
filaments = armc.generate_starting_config!(lattice, Nfil, Nsca, lf)
sysparms = armc.SystemParams(kd, ks, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
simparms = armc.SimulationParams(iters, steps, max_bias_diff, write_interval,
                                 radius_move_freq, filebase)
system = armc.System(sysparms, filaments, max_radius)
armc.run!(system, lattice, simparms)
#armc.run_us!(system, lattice, simparms)
