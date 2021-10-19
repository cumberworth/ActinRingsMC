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

steps = 1e4
write_interval = 1e2
radius_move_freq = 0.5
op_filename = "test.ops"

lf = armc.calc_lf(Lf, delta)
Lf = armc.calc_Lf(lf, delta)
lattice_height = armc.calc_lattice_height(Nsca, lf)
max_radius = Lf*Nsca

lattice = armc.Lattice(lattice_height)
filaments = armc.generate_starting_config!(lattice, Nfil, Nsca, lf)
sysparms = armc.SystemParams(kd, ks, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
simparms = armc.SimulationParams(steps, write_interval, radius_move_freq, op_filename)
system = armc.System(sysparms, filaments, max_radius)
op_file = armc.prepare_op_file(op_filename)
armc.run(system, lattice, op_file, simparms)
close(op_file)
