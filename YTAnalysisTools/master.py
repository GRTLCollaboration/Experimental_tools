# master.py
# James Widdicombe
# Last Updated 01/08/2018
# Last Formatted Dec 2019
# A master analysis script that include L2M, L2H and mass measurments
# for GRChombo HDF5 files, specifically for stars with density rho

# Load the modules
import yt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Enable Parallelism
yt.enable_parallelism()

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

# CHANGE ME
# Loading dataset (Load from one folder up)
ds = yt.load("../outMatterSF_*.3d.hdf5")

# Matplotlib Settings
rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]

# Define H2
def _H2(field, data):
    return data["Ham"] * data["Ham"]


# Define M2
def _M2(field, data):
    return (
        data["Mom1"] * data["Mom1"]
        + data["Mom2"] * data["Mom2"]
        + data["Mom3"] * data["Mom3"]
    )


# Define Modified Rho
def _RhoJam(field, data):
    return data["rho"] / (data["chi"] * np.sqrt(data["chi"]))


# Arrays for output data
time_data = []
L2H_data = []
L2M_data = []
RhoJam_Average_data = []
Mass_Total_data = []
star_size_initial_x_data = []
SizeyData = []
max_rho_pos_x_data = []
max_rho_pos_y_data = []
star_mass_center_x = []
star_mass_center_y = []
star_mass_center_average = []
star_size_x_precise = []
star_size_y_precise = []
star_size_z_precise = []
star_mass_ellipse = []

# Other factors
total_box_size = float(ds[0].domain_right_edge[0])
center = total_box_size / 2.0
Quality = "cubic"

# Plot Start Time
plot_start_time = time.time()
program_runtime_data = []
program_frames = []
counter = 0

for i in ds:

    # Timings
    i_start = time.time()
    counter += 1
    program_frames.append(counter)

    # Current Simulation Time
    time_data.append(i.current_time)

    # Add the M2 and L2 Fields
    i.add_field("H2", _H2, units="")
    i.add_field("M2", _M2, units="")

    # Add the RhoJam Field to the data
    i.add_field("RhoJam", _RhoJam, units="")

    # All Data
    ad = i.all_data()

    # L2H
    meanH2 = ad.mean("H2", weight="cell_volume")
    L2H = np.sqrt(meanH2)
    L2H_data.append(L2H)

    # L2M
    meanM2 = ad.mean("M2", weight="cell_volume")
    L2M = np.sqrt(meanM2)
    L2M_data.append(L2M)

    # Average RhoJam
    Rho_Jam_Average = ad.mean("RhoJam", weight="cell_volume")
    RhoJam_Average_data.append(Rho_Jam_Average)

    # Total RhoJam
    Mass_Total = Rho_Jam_Average * (total_box_size ** 3)
    Mass_Total_data.append(Mass_Total)

    program_runtime_data.append(time.time() - i_start)

    if yt.is_root():
        np.savetxt("program_runtime_data.out", program_runtime_data)
        np.savetxt("L2M.out", L2M_data)
        np.savetxt("L2H.out", L2H_data)
        np.savetxt("RhoJam_Average.out", RhoJam_Average_data)
        np.savetxt("Mass_Total.out", Mass_Total_data)
        np.savetxt("time_data.out", time_data)

    # Look at initial Stars
    # Line out from 0 to center (x), at the center of (y and z)
    c = i.ray((0, center, center), (center, center, center))

    # Name x coordinate and rho from lineout
    x = c["x"]
    rho = c["rho"]

    # Find maximum rho value within lineout
    rho_max = float(rho.max())

    # Find the coordinate of where the maximum value is
    x_max = np.where(rho == rho_max)[0][0]
    x_max_val = x[x_max]
    max_rho_pos_x_data.append(x_max_val)

    if yt.is_root():
        np.savetxt("max_rho_pos_x_data.out", max_rho_pos_x_data)

    # Create a function that goes to negative when we are below 95% of rho max
    rhofun = interp1d(x, rho - 0.05 * rho_max, kind=Quality, fill_value="extrapolate")

    # Solve the function for where it goes to zero using a best guess past and
    # before the the max value of x
    x_1 = fsolve(rhofun, float(x[x_max]) - 1)[0]
    x_2 = fsolve(rhofun, float(x[x_max]) + 1)[0]

    # Size of the 95% rho in center
    size_x = (x_2 - x_1) / 2.0
    star_size_initial_x_data.append(size_x)
    if yt.is_root():
        np.savetxt("star_size_inimax_rho_postial_x.out", star_size_initial_x_data)

    # If the stars have merged, we can look at the mass of the star
    if x_max_val > center - 3.5 and x_max_val < center + 3.5:

        # Measure in x plane
        c_precise = i.ray((0, center, center), (total_box_size, center, center))
        x_precise = c_precise["x"]
        rho_precise = c_precise["rho"]
        rho_max_precise = float(rho_precise.max())
        x_max_precise = np.where(rho_precise == rho_max_precise)[0][0]
        x_max_val_precise = x_precise[x_max_precise]
        rhofun_precise = interp1d(
            x_precise,
            rho_precise - 0.05 * rho_max_precise,
            kind=Quality,
            fill_value="extrapolate",
        )
        x_1_precise = fsolve(rhofun_precise, float(x_precise[x_max_precise]) - 1)[0]
        x_2_precise = fsolve(rhofun_precise, float(x_precise[x_max_precise]) + 1)[0]
        size_x_precise = (x_2_precise - x_1_precise) / 2.0
        star_size_x_precise.append(size_x_precise)

        # Measure in y plane
        cy_precise = i.ray(
            (x_max_val_precise, 0, center), (x_max_val_precise, total_box_size, center)
        )
        y_precise = cy_precise["y"]
        rhoy_precise = cy_precise["rho"]
        rho_max_precise = float(rhoy_precise.max())
        y_max_precise = np.where(rhoy_precise == rho_max_precise)[0][0]
        y_max_val_precise = y_precise[y_max_precise]
        rhofuny_precise = interp1d(
            y_precise,
            rhoy_precise - 0.05 * rho_max_precise,
            kind=Quality,
            fill_value="extrapolate",
        )
        y_1_precise = fsolve(rhofuny_precise, float(y_precise[y_max_precise]) - 1)[0]
        y_2_precise = fsolve(rhofuny_precise, float(y_precise[y_max_precise]) + 1)[0]
        size_y_precise = (y_2_precise - y_1_precise) / 2.0
        star_size_y_precise.append(size_y_precise)

        # Spherical Mass Estimates
        sp = i.sphere(
            [float(x_max_val_precise), float(y_max_val_precise), center],
            max(size_x_precise, 1.0),
        )
        RhoJamTotalx = (
            float(sp.mean("RhoJam", weight="cell_volume"))
            * (4.0 / 3.0)
            * np.pi
            * np.power((size_x_precise), 3)
        )
        star_mass_center_x.append(RhoJamTotalx)
        sp2 = i.sphere(
            [float(x_max_val_precise), float(y_max_val_precise), center],
            max(size_y_precise, 1.0),
        )
        RhoJamTotaly = (
            float(sp2.mean("RhoJam", weight="cell_volume"))
            * (4.0 / 3.0)
            * np.pi
            * np.power((size_y_precise), 3)
        )
        star_mass_center_y.append(RhoJamTotaly)

        # Elipsoidal Mass
        spp = i.sphere(
            [float(x_max_val_precise), float(y_max_val_precise), center],
            max(size_x_precise, size_y_precise, 1.0),
        )
        Rho5 = 0.05 * float(spp["RhoJam"].max())
        cr = spp.cut_region("obj['RhoJam'] > " + str(Rho5))
        total_vol = cr.sum("cell_volume")
        mass_precise = float(cr.mean("RhoJam", weight="cell_volume")) * total_vol
        star_mass_ellipse.append(mass_precise)

        # Plot Star Shape
        if yt.is_root():
            plt.figure(figsize=(20, 20))
            plt.scatter(
                x_precise, rho_precise, alpha=alp, label="t =" + str(i.current_time)
            )
            plt.scatter(
                y_precise, rhoy_precise, alpha=alp, label="t =" + str(i.current_time)
            )
            plt.xlim([x_1 - 10, x_2 + 10])
            plt.legend(loc="upper right")
            plt.xlabel(r"$x~[1/m]$")
            plt.ylabel(r"$\rho~[M_{pl}^2 m^2]$")
            plt.savefig(("rho%05d.png" % counter), bbox_inches="tight")
            plt.close()

    # Append zeros if there is no central star
    else:
        star_size_x_precise.append(0.0)
        star_size_y_precise.append(0.0)
        star_mass_center_x.append(0.0)
        star_mass_center_y.append(0.0)
        star_mass_ellipse.append(0.0)

    # Save data
    if yt.is_root():
        np.savetxt("star_size_x_precise.out", star_size_x_precise)
        np.savetxt("star_size_y_precise.out", star_size_y_precise)
        np.savetxt("star_mass_center_x.out", star_mass_center_x)
        np.savetxt("star_mass_center_y.out", star_mass_center_y)
        np.savetxt("star_mass_ellipse.out", star_mass_ellipse)

if yt.is_root():
    # Plots
    # Update rcParams so that we get all axis labelling
    rcParams.update({"figure.autolayout": True})
    rcParams["axes.formatter.limits"] = [-5, 5]

    # Timing Plot
    plt.figure(1)
    plt.title("Program Run Time")
    plt.plot(program_frames, program_runtime_data)
    plt.xlabel("Data Files")
    plt.ylabel("Time $[s]$")
    plt.savefig("pictures_program_time.png")
    plt.close()

    # L2H
    plt.figure(2)
    plt.plot(time_data, L2H_data)
    # plt.title('$\\mathcal{H}$ vs time')
    plt.ylabel("$\\mathcal{H}$")
    plt.xlabel("Time $[1/m]$")
    plt.savefig("H.png")
    plt.close()

    # L2M
    plt.figure(3)
    plt.plot(time_data, L2M_data)
    # plt.title('$\\mathcal{M}$ vs time')
    plt.ylabel("$\\mathcal{M}$")
    plt.xlabel("Time $[1/m]$")
    plt.savefig("M.png")
    plt.close()

    # RhoJam_Average
    plt.figure(4)
    plt.plot(time_data, RhoJam_Average_data)
    plt.title("$\\tilde{\\rho}$ vs time")
    plt.ylabel("$\\tilde{\\rho} \\, [M_{pl}^2m^2]$")
    plt.xlabel("Time $[1/m]$")
    plt.savefig("RhoJam_Average.png")
    plt.close()

    # Star Size
    plt.figure(5)
    plt.plot(time_data, star_size_x_precise, label="x size")
    plt.plot(time_data, star_size_y_precise, label="y size")
    plt.title("$r$ vs time")
    plt.ylabel("$r \\, [1/m]$")
    plt.xlabel("Time $[1/m]$")
    plt.legend()
    plt.savefig("Star_Size.png")
    plt.close()

    # Star Size Diagnostic
    plt.figure(6)
    plt.plot(time_data, star_size_x_precise, label="x size")
    plt.plot(time_data, star_size_y_precise, label="y size")
    plt.plot(time_data, star_size_initial_x_data, label="Inital x")
    plt.title("$r$ vs time")
    plt.ylabel("$r \\, [1/m]$")
    plt.xlabel("Time $[1/m]$")
    plt.legend()
    plt.savefig("Star_Size_Diagnostic.png")
    plt.close()

    # Total Mass Diagnostic
    plt.figure(7)
    plt.plot(time_data, Mass_Total_data, label="Total Mass")
    plt.plot(time_data, star_mass_center_y, label="Y Mass")
    plt.plot(time_data, star_mass_center_x, label="X Mass")
    plt.plot(time_data, star_mass_ellipse, label="Ellipsiodal Mass")
    plt.title("$M$ vs time")
    plt.ylabel("$M \\, [M_{pl}^2/m]$")
    plt.xlabel("Time $[1/m]$")
    plt.legend()
    plt.savefig("Mass_Total_Diag.png")
    plt.close()

    # Total Mass
    plt.figure(8)
    plt.plot(time_data, Mass_Total_data, label="Total Mass")
    plt.plot(time_data, star_mass_ellipse, label="Star Mass")
    # plt.title('$M$ vs time')
    plt.ylabel("$M \\, [M_{pl}^2/m]$")
    plt.xlabel("Time $[1/m]$")
    plt.legend()
    plt.savefig("Mass_Total.png")
    plt.close()

    # Max rho pos
    plt.figure(9)
    plt.plot(time_data, max_rho_pos_x_data)
    plt.title("$x,y$ vs time")
    plt.xlabel("Time $[1/m]$")
    plt.savefig("max_rho_pos.png")
    plt.close()
