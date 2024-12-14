####################################################################################################
# Simple PK/PD modeller
####################################################################################################
# * Drug concentration in system reduces exponentially with time, with fixed rate constant
# * Drug added to system at constant rate over given time period after administration
#
# * Assumes that absorption is complete before next dose administered
#
# * Outputs total drug in system vs time
####################################################################################################

import math
import numpy as np
import matplotlib.pyplot as plt

####################################################################################################
# EDIT THIS SECTION TO CHANGE INPUTS
####################################################################################################

g_half_life = 21           # Drug half life (h)
g_absorption_time = 1      # Dose absorption time (h)

g_output_dt = 1            # Output time step (h)

def administer_doses():
 """Edit this function to change dosing schedule. Use dose()."""
 for i in range(30):
  dose(1, 24)

 for i in range(4):
  dose(1, 2*24)
  dose(1, 2*24)
  dose(1, 3*24)




####################################################################################################
# FUNCTIONS
####################################################################################################

def add_output_data(t, conc):
 global g_t_wp, g_conc_wp, g_next_op, g_count_data_points, g_num_data_points, g_filtered
  
 if g_count_data_points:
  g_num_data_points += 1
  return
 else:
  # Filter duplicates
  if g_next_op>=1 and g_t[g_next_op-1]==t and g_conc[g_next_op-1]==conc:
   g_filtered += 1
   return
  
  g_t[g_next_op] = t
  g_conc[g_next_op] = conc
 
  g_next_op += 1




def time_step(dose, dt, output_dt=None):
 """Record concentrations, for dose administered between most recent time waypoint and most recent 
 time waypoint+dt
 
 dose               Dose to administer
 dt                 Time over which to administer dose. May be 0.
 output_dt          Time interval at which to output concentrations (default: g_output_dt)
 
 Calculation:
  - t time
  - c(t) total drug in system at time t
  - k rate constant
  - r rate of addition of drug to system

       dc(t)/dt = -kc + r
   =>      c(t) = r/k + [c(0)-r/k]*exp(-kt)
 """
 global g_t_wp, g_conc_wp, g_k, g_output_dt
 
 if output_dt is None:
  output_dt = g_output_dt
 
 if dt==0:
  g_conc_wp += dose
  add_output_data(g_t_wp, g_conc_wp)
 else:
  rate = dose/dt

  # Intermediate t, conc
  t = math.ceil(g_t_wp/output_dt)*output_dt     # next output time >= g_t_wp
  t_end = g_t_wp+dt
  conc_0 = g_conc_wp
 
  while t<t_end:
   conc = rate/g_k + (conc_0-rate/g_k)*math.exp(-g_k*(t-g_t_wp))
   add_output_data(t,conc)
  
   t += output_dt
 
  g_t_wp = t_end
  g_conc_wp = rate/g_k + (conc_0-rate/g_k)*math.exp(-g_k*dt)
 
  add_output_data(g_t_wp, g_conc_wp)
 
 
 
 
def dose(dose, dose_interval, absorption_time=None, output_dt=None):
 """Simulate time between administering dose and waiting for next dose
 
 dose               Dose to administer
 dose_interval      Time between taking this dose and taking next dose (h)
 absorption_time    Time for dose to absorb (h). Must be < dose_interval. May be 0.
                    (Default: g_absorption_time)
 output_dt          Output time interval (defaults to global output_interval)
 """
 if absorption_time is None:
  absorption_time = g_absorption_time
  
 if absorption_time>dose_interval:
  raise Exception("Absorption time must be <= dose interval")
 
 time_step(dose, absorption_time, output_dt)
 time_step(0, dose_interval-absorption_time, output_dt)




####################################################################################################
# MAIN ENTRY POINT
####################################################################################################

g_k = math.log(2)/g_half_life       # Time constant
g_num_data_points = 0

# Dummy run to find number of data points

g_t_wp = 0                             # Waypoint time (hours)
g_conc_wp = 0                          # Waypoint drug concentration
g_count_data_points = True

administer_doses()

# Run it for real this time

g_t = np.zeros(g_num_data_points)
g_conc = np.zeros(g_num_data_points)

g_t_wp = 0
g_conc_wp = 0
g_next_op = 0                       # Next output index
g_filtered = 0                      # Number of duplicate outputs filtered
g_count_data_points = False

administer_doses()

# Output data on graph

num_data_points = g_num_data_points-g_filtered
plt.plot(g_t[0:num_data_points]/24, g_conc[0:num_data_points])

plt.xlabel("Time (days)")
plt.ylabel("Total drug in system")

plt.show()