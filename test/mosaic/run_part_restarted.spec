run_type particle               # particle-resolved run
output_prefix out/mosaic_restarted # prefix of output files
n_loop 1                        # number of Monte Carlo loops
n_part 2                        # total number of particles
kernel golovin                  # coagulation kernel
restart yes                     # whether to restart from saved state (yes/no)
restart_file out/mosaic_0001_00000013.nc # saved state file to restart from

t_max 43200                     # total simulation time (s)
del_t 300                       # timestep (s)
t_output 3600                   # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
r_min 1e-8                      # minimum radius (m)
r_max 1e-3                      # maximum radius (m)

gas_data gas_data.dat           # file containing gas data

aerosol_data aero_data.dat      # file containing aerosol data

temp_profile temp.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rand_init 0                     # random initialization (0 to auto-generate)
mix_rate 0                      # mixing rate between processes (0 to 1)
do_coagulation no               # whether to do coagulation (yes/no)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
record_removals no              # whether to record particle removals (yes/no)
