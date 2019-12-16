# parallel_line.py
# James Widdicombe
# Last Updated 18/10/2018
# Last Formatted Dec 2019
# lineout

# Load the modules
import yt
import numpy as np
from scipy.interpolate import interp1d
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
data_location = "../../plt*.hdf5"  # Data file location

# Loading dataset
ts = yt.load(data_location)

outputdirectory = "pictures/"

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

    e = i.ray((0, center, center), (center, center, center))
    ex = e["x"]
    rhoe = e["rho"]
    rhoxe_max = float(rhoe.max())
    x_maxe = np.where(rhoe == rhoxe_max)[0][0]
    x_max_vale = ex[x_maxe]

    c = i.ray((0, center, center), (total_box_size, center, center))
    d = i.ray((x_max_vale, 0, center), (x_max_vale, total_box_size, center))

    # Name x coordinate and rho from lineout
    x = c["x"]
    rho = c["rho"]
    y = d["y"]
    rhoy = d["rho"]

    rhoy_max = float(rhoy.max())
    y_max = np.where(rhoy == rhoy_max)[0][0]
    y_max_val = y[y_max]

    # fx = interp1d(x,rho, fill_value="extrapolate")
    # fy = interp1d(y,rhoy, fill_value="extrapolate")

    fx_jam = interp1d(x - x_max_vale, rho, fill_value="extrapolate")
    fy_jam = interp1d(y - y_max_val, rhoy, fill_value="extrapolate")

    linda = np.linspace(-40, 40, num=100)

    eugene = []

    for q in linda:
        temp = fx_jam(q) - fy_jam(q)
        eugene.append(temp)
    if float(i.current_time) < 10.0:
        naming = "0000" + str(float(i.current_time))
    elif float(i.current_time) < 100.0 and float(i.current_time) > 9.9999999:
        naming = "000" + str(float(i.current_time))
    elif float(i.current_time) < 1000.0 and float(i.current_time) > 99.9999999:
        naming = "00" + str(float(i.current_time))
    else:
        naming = "0" + str(float(i.current_time))
    plt.figure(figsize=(20, 20))
    # plt.plot(x,fx(x), alpha = alp,label = "t ="+str(i.current_time),color='blue')
    # plt.scatter(x,rho, alpha = alp,label = "t ="+str(i.current_time),color='blue')
    # plt.plot(y-y_max_val,fy(y), alpha = alp,label = "t ="+str(i.current_time),color='red')
    # plt.plot(x_s,rho_e, alpha = alp,label = "t ="+str(i.current_time),color='green')
    # plt.plot(x-x_max_vale,fx(x), alpha = alp,label = "t ="+str(i.current_time),color='blue')
    # plt.plot(linda,fx_jam(linda), alpha = alp,label = "t ="+str(i.current_time),color='blue')
    # plt.plot(linda,fy_jam(linda), alpha = alp,label = "t ="+str(i.current_time),color='red')
    plt.plot(linda, eugene, alpha=alp, label="t =" + str(i.current_time), color="blue")
    # plt.xlim([center-35,center+35])
    plt.xlim([0 - 35, 0 + 35])
    plt.ylim(-0.001, 0.001)
    plt.legend(loc="upper right")
    plt.xlabel(r"$x~[1/m]$")
    plt.ylabel(r"$\rho~[M_{pl}^2 m^2]$")
    plt.savefig(outputdirectory + naming + ".png", bbox_inches="tight")
    plt.close()
