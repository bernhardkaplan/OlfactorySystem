import os
import numpy as np

gor_min_min_sweep = 5e-5
gor_min_max_sweep = 3e-4
gor_max_min_sweep = 1e-3
gor_max_max_sweep = 5e-3

gkcag_min_min_sweep = 5e-4
gkcag_min_max_sweep = 5e-3
gkcag_max_min_sweep = 5e-3
gkcag_max_max_sweep = 5e-2

gcal_min_min_sweep = 5e-5
gcal_min_max_sweep = 5e-4
gcal_max_min_sweep = 5e-4
gcal_max_max_sweep = 5e-3

gleak_min_min_sweep = 5e-5
gleak_min_max_sweep = 2e-4
gleak_max_min_sweep = 3e-4
gleak_max_max_sweep = 1e-3

#gor_max_sweep = 1e-3
#gkcag_min_sweep = 1e-3
#gkcag_max_sweep = 1e-2
#gcal_min_sweep = 1e-4
#gcal_max_sweep = 1e-3
#gleak_min_sweep = 1e-4
#gleak_max_sweep = 1e-5

n_gor_steps = 3
gor_min_range = np.linspace(gor_min_min_sweep, gor_min_max_sweep, n_gor_steps)
gor_max_range = np.linspace(gor_max_min_sweep, gor_max_max_sweep, n_gor_steps)
n_gkcag_steps = 3
gkcag_min_range = np.linspace(gkcag_min_min_sweep, gkcag_min_max_sweep, n_gkcag_steps)
gkcag_max_range = np.linspace(gkcag_max_min_sweep, gkcag_max_max_sweep, n_gkcag_steps)
n_gcal_steps = 3
gcal_min_range = np.linspace(gcal_min_min_sweep, gcal_min_max_sweep, n_gcal_steps)
gcal_max_range = np.linspace(gcal_max_min_sweep, gcal_max_max_sweep, n_gcal_steps)
n_gleak_steps = 3
gleak_min_range = np.linspace(gleak_min_min_sweep, gleak_min_max_sweep, n_gleak_steps)
gleak_max_range = np.linspace(gleak_max_min_sweep, gleak_max_max_sweep, n_gleak_steps)

sim_cnt = 0

log_file = file('orn_sweep_params.log', 'w')
for gor_min in gor_min_range:
    for gor_max in gor_max_range:
        for gkcag_min in gkcag_min_range:
            for gkcag_max in gkcag_max_range:
                for gcal_min in gcal_min_range:
                    for gcal_max in gcal_max_range:
                        for gleak_min in gleak_min_range:
                            for gleak_max in gleak_max_range:
                                script_cmd = 'python HandTuneOrnParameters.py %d %f %f %f %f %f %f %f %f' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max)
                                log_info = '%d %f %f %f %f %f %f %f %f\n' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max)
                                log_file.write(log_info)
                                log_file.flush()
                                print script_cmd
                                os.system(script_cmd)
                                sim_cnt += 1
                #gor_min_max = np.array([7e-5, 1e-3])
#script_cmd = 'python HandTuneOrnParameters.py %d %f %f %f %f %f %f %f %f' % (sim_cnt, gor_min, gor_max, gkcag_min, gkcag_max, gcal_min, gcal_max, gleak_min, gleak_max)

