using ActinRingsMC

# System variables
ks = 7.82e-9
kd = 7.82e-9/5
T = 300
delta = 36e-9
Xc = 12e-9
EI = 3.6e-26
lf = 84
Nfil = 6
Nsca = 4

# Best to set lf directly and calculate the nearest Lf
Lf = calc_Lf(lf, delta)

# Calculate the lattice parameters
max_height = calc_max_lattice_height(Nsca, lf)
min_height = calc_min_lattice_height(Nsca, lf)

# This sets how many sites overlap between the filaments in the initial configuration
overlap = 1
start_height = max_height - Nsca*(overlap - 1)
start_radius = calc_radius(delta, start_height)
lattice = Lattice(start_height, max_height, min_height)

filaments = generate_starting_config(lattice, Nfil, Nsca, lf, overlap)

# This is required to make sure internal of lattice matches the starting configuration
update_occupancies!(filaments, lattice)

sysparms = SystemParams(ks, kd, T, delta, Xc, EI, Lf, lf, Nfil, Nsca)
system = System(sysparms, filaments, start_radius)

# Simulation parameters
iters = 2
steps = 1e4
max_bias_diff = 5
write_interval = 1e3
radius_move_freq = 0.5
filebase = "outs/test_run-0"
analytical_biases = true
read_biases = false
biases_filename = ""
restart_iter = 0
binwidth = 1

simparms = SimulationParams(
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

mkpath("outs")

# First run
run_us!(system, lattice, simparms)

# Update simulation parameters for continuation run
filebase = "outs/test_run-1"
read_biases = true
biases_filename = "outs/test_run-0.biases"
restart_iter = 2
simparms = SimulationParams(
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

# Continuation run
run_us!(system, lattice, simparms)
