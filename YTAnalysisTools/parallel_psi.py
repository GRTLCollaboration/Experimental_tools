# parallel_psi.py
# James Widdicombe
# Last Updated 17/10/2018
# Last Formatted Dec 2019
# Calculate psi evolution

# Load the modules
import yt
import numpy as np
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Timings
start_time = time.time()

# Enable Parallelism
yt.enable_parallelism()

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]
rcParams["font.size"] = 12

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../../outMatterSF_00*"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Program Parameters
center = ts[0].domain_right_edge / 2.0
cutoff = 60
adjusted_right = int(center[0]) * 2 - cutoff

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):

    # All Data
    ad = i.r[cutoff:adjusted_right, cutoff:adjusted_right, cutoff:adjusted_right]

    # L2H
    maxpsi = ad.max("phi")

    # L2M
    minpsi = ad.min("phi")

    array = [i.current_time, maxpsi, minpsi]

    sto.result = array
    sto.result_id = str(i)
if yt.is_root():
    timedata = []
    maxpsidata = []
    minpsidata = []

    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        maxpsidata.append(L[1][1])
        minpsidata.append(L[1][2])
    # L2H
    plt.figure(1)
    plt.plot(timedata, maxpsidata)
    plt.plot(timedata, minpsidata)
    plt.ylabel("$\\phi$ $[M_{pl}]$")
    plt.xlabel("Time $[1/m]$")
    plt.grid()
    plt.savefig("phi.png", bbox_inches="tight")
    plt.close()

    np.savetxt("time.out", timedata)
    np.savetxt("maxpsi.out", maxpsidata)
    np.savetxt("minpsi.out", minpsidata)
