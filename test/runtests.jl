using Test

import ActinRingsMC as armc

kd = 7.82e-9
ks = 7.82e-9/5
T = 300
delta = 36e-9
Xc = 12e-9
EI = 3.6e-26
Lf = 2.988e-6
lf = 4
Nfil = 6
Nsca = 4

lattice = armc.Lattice(Nsca*lf - Nsca - 1)
filaments = armc.generate_starting_config!(lattice, Nfil, Nsca, lf)
sysparms = armc.SystemParams(kd, ks, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
system = armc.System(filaments, 2.988e-6*Nsca, sysparms)
