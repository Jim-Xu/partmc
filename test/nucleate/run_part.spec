run_type particle               # particle-resolved run
output_prefix out/nucleate_part # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 1000                     # total number of particles
kernel sedi                     # coagulation kernel
nucleate sulf_acid              # nucleation parameterization
restart no                      # whether to restart from saved state (yes/no)

t_max 3600                      # total simulation time (s)
del_t 1                         # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
d_min 1e-10                     # minimum diameter (m)
d_max 1e-6                      # maximum diameter (m)

weight none                     # weighting function

gas_data gas_data.dat           # file containing gas data
gas_init gas_init.dat           # initial gas mixing ratios

aerosol_data aero_data.dat      # file containing aerosol data
aerosol_init aero_init_dist.dat # aerosol initial condition file

temp_profile temp.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.999              # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

rand_init 0                     # random initialization (0 to auto-generate)
do_coagulation no               # whether to do coagulation (yes/no)
allow_doubling no               # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
