# Last Formatted Dec 2019
# Load the modules
import yt
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams


# Enable Parallelism
yt.enable_parallelism()

# Timings
start_time = time.time()

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../../outMatterSF*.hdf5"  # Data file location

# Loading dataset
ts = yt.load(data_location)

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

# Other factors
total_box_size = float(ts[0].domain_right_edge[0])
center = total_box_size / 2.0
Quality = "cubic"

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):
    # Timings
    L_start = time.time()

    # Look at initial Stars
    # Line out from 0 to center (x), at the center of (y and z)

    e = i.ray((0, center, center), (center * 2, center, center))
    ex = e["x"]
    rhoe = e["psi"]

    if float(i.current_time) < 10.0:
        naming = "0000" + str(float(i.current_time))
    elif float(i.current_time) < 100.0 and float(i.current_time) > 9.9999999:
        naming = "000" + str(float(i.current_time))
    elif float(i.current_time) < 1000.0 and float(i.current_time) > 99.9999999:
        naming = "00" + str(float(i.current_time))
    else:
        naming = "0" + str(float(i.current_time))
    plt.figure(figsize=(20, 20))
    plt.scatter(ex, rhoe, alpha=alp, label="t =" + str(i.current_time), color="blue")
    plt.xlim([center - 35, center + 35])
    # plt.ylim(-0.001,0.001)
    plt.legend(loc="upper right")
    plt.xlabel(r"$x~[1/m]$")
    plt.ylabel(r"$\psi$")
    plt.savefig(naming + ".png", bbox_inches="tight")
    plt.close()
