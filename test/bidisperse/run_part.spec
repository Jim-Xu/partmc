run_type particle               # particle-resolved run
output_prefix out/bidisperse_part # prefix of output files
n_loop 1                        # number of Monte Carlo loops
n_part 10001                    # number of Monte Carlo particles
kernel sedi                     # coagulation kernel
restart no                      # whether to restart from saved state (yes/no)

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 10                     # output interval (0 disables) (s)
t_progress 10                   # progress printing interval (0 disables) (s)

n_bin 255                       # number of bins
r_min 1e-8                      # minimum radius (m)
r_max 1e0                       # maximum radius (m)

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
mix_rate 0                      # mixing rate between processes (0 to 1)
do_coagulation yes              # whether to do coagulation (yes/no)
allow_doubling no               # whether to allow doubling (yes/no)
allow_halving no                # whether to allow halving (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
record_removals no              # whether to record particle removals (yes/no)
