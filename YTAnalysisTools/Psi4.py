# Psi4.py
# James Widdicombe
# Modified from a code written by Thomas Helfer
# Last Updated 14/08/2018
# Last Formatted Dec 2019
# Decomposition of Psi4 into spin weighted spherical harmonics
# l = 2,3,4

# Import the modules
import yt
import numpy as np
from numpy import pi
import time

# Enable YT Parallelism
yt.enable_parallelism()

# Script Parameters
extraction_radius = 120  # radius of extraction
data_location = "../outMatterSF_*.3d.hdf5"  # Data file location

# Loading dataset
ds = yt.load(data_location)

# Initialize Output arrays
# l = 2
Weyl4_l2_m0_data = []
# positive m
Weyl4_l2_m1_data = []
Weyl4_l2_m2_data = []
# negative m
Weyl4_l2_m1n_data = []
Weyl4_l2_m2n_data = []
# l = 3
Weyl4_l3_m0_data = []
# positive m
Weyl4_l3_m1_data = []
Weyl4_l3_m2_data = []
Weyl4_l3_m3_data = []
# negative m
Weyl4_l3_m1n_data = []
Weyl4_l3_m2n_data = []
Weyl4_l3_m3n_data = []
# l = 4
Weyl4_l4_m0_data = []
# positive m
Weyl4_l4_m1_data = []
Weyl4_l4_m2_data = []
Weyl4_l4_m3_data = []
Weyl4_l4_m4_data = []
# negative m
Weyl4_l4_m1n_data = []
Weyl4_l4_m2n_data = []
Weyl4_l4_m3n_data = []
Weyl4_l4_m4n_data = []

# Program Parameters
center = ds[0].domain_right_edge / 2.0

# Initialize Arrays
timedata = []
program_runtime_data = []

# Definitions for Quadrature scheme
N = 131

coord = np.loadtxt("PointDistFiles/lebedev/lebedev_%03d.txt" % N)
theta = coord[:, 1] * pi / 180
phi = coord[:, 0] * pi / 180 + pi
w = coord[:, 2]

phi_length = len(phi)

# Spinweighted spherical harmonics s = -2
# =========================================
# L = 2
# =========================================
# m = 0 zero
sY_l2_m0 = np.sqrt(15.0 / (32.0 * np.pi)) * np.sin(theta) ** 2
# positive m values : 1,2
sY_l2_m1 = (
    np.exp(1j * phi)
    * np.sqrt(5.0 / (np.pi))
    * np.cos(theta / 2) ** 3
    * np.sin(theta / 2)
)
sY_l2_m2 = (
    1.0 / 2.0 * np.exp(2 * 1j * phi) * np.sqrt(5.0 / (np.pi)) * np.cos(theta / 2) ** 4
)
# negative m values :-1,-2
sY_l2_m1n = (
    1.0
    / 2.0
    * np.exp(-1j * phi)
    * np.sqrt(5.0 / (np.pi))
    * np.sin(theta / 2) ** 2
    * np.sin(theta)
)
sY_l2_m2n = (
    1.0 / 2.0 * np.exp(-2 * 1j * phi) * np.sqrt(5.0 / (np.pi)) * np.sin(theta / 2) ** 4
)
# =========================================
# L = 3
# =========================================
# m = 0 zero
sY_l3_m0 = (
    1.0 / 4.0 * np.sqrt(105.0 / (2.0 * np.pi)) * np.cos(theta) * np.sin(theta) ** 2
)
# positive m values : 1, 2 ...
sY_l3_m1 = (
    1.0
    / 2.0
    * np.exp(1j * phi)
    * np.sqrt(35.0 / (2.0 * np.pi))
    * np.cos(theta / 2) ** 3
    * (-1 + 3 * np.cos(theta))
    * np.sin(theta / 2)
)
sY_l3_m2 = (
    1.0
    / 2.0
    * np.exp(2 * 1j * phi)
    * np.sqrt(7.0 / (np.pi))
    * np.cos(theta / 2) ** 4
    * (-2 + 3 * np.cos(theta))
)
sY_l3_m3 = (
    -np.exp(3 * 1j * phi)
    * np.sqrt(21.0 / (2 * np.pi))
    * np.cos(theta / 2) ** 5
    * np.sin(theta / 2)
)
# negative m values : -1, -2 ...
sY_l3_m1n = (
    1.0
    / 2.0
    * np.exp(-1j * phi)
    * np.sqrt(35.0 / (2 * np.pi))
    * np.cos(theta / 2)
    * (1 + 3 * np.cos(theta))
    * np.sin(theta / 2) ** 3
)
sY_l3_m2n = (
    1.0
    / 2.0
    * np.exp(-2 * 1j * phi)
    * np.sqrt(7.0 / (np.pi))
    * (2 + 3 * np.cos(theta))
    * np.sin(theta / 2) ** 4
)
sY_l3_m3n = (
    1.0
    / 2.0
    * np.exp(-3 * 1j * phi)
    * np.sqrt(21.0 / (2 * np.pi))
    * np.sin(theta / 2) ** 4
    * np.sin(theta)
)
# =========================================
# L = 4
# =========================================
sY_l4_m0 = (
    3.0
    / 16.0
    * np.sqrt(5.0 / (2.0 * np.pi))
    * (5.0 + 7.0 * np.cos(2.0 * theta))
    * np.sin(theta) ** 2
)
# positive m = 1,2,3,4
sY_l4_m1 = (
    3
    / (2 * np.sqrt(2 * np.pi))
    * np.exp(1j * phi)
    * np.cos(theta / 2) ** 3
    * (6 - 7 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.sin(theta / 2)
)
sY_l4_m2 = (
    3
    / (4 * np.sqrt(np.pi))
    * np.exp(2.0 * 1j * phi)
    * (9 - 14 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.cos(theta / 2) ** 4
)
sY_l4_m3 = (
    -3.0
    * np.sqrt(7.0 / (2.0 * np.pi))
    * np.cos(theta / 2) ** 5
    * np.exp(3.0 * 1j * phi)
    * (-1 + 2 * np.cos(theta))
    * np.sin(theta / 2)
)
sY_l4_m4 = (
    3.0
    * np.exp(4 * 1j * phi)
    * np.sqrt(7 / np.pi)
    * np.cos(theta / 2) ** 6
    * np.sin(theta / 2) ** 2
)
# negative m = -1,-2,-3,-4
sY_l4_m1n = (
    3.0
    / (2.0 * np.sqrt(2 * np.pi))
    * np.exp(-1j * phi)
    * np.cos(theta / 2)
    * (6 + 7 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.sin(theta / 2) ** 3
)
sY_l4_m2n = (
    3.0
    / (4.0 * np.sqrt(np.pi))
    * np.exp(-2 * 1j * phi)
    * (9 + 14 * np.cos(theta) + 7 * np.cos(2 * theta))
    * np.sin(theta / 2) ** 4
)
sY_l4_m3n = (
    3.0
    * np.sqrt(7.0 / (2.0 * np.pi))
    * np.exp(-3 * 1j * phi)
    * (1 + 2 * np.cos(theta))
    * np.sin(theta / 2) ** 5
    * np.cos(theta / 2)
)
sY_l4_m4n = (
    3.0
    / 4.0
    * np.exp(-4 * 1j * phi)
    * np.sqrt(7 / np.pi)
    * np.sin(theta / 2) ** 4
    * np.sin(theta) ** 2
)

# ==============================
#         Loop over all frames
# ==============================
for i in ds:

    # Timings
    i_start = time.time()

    # Initalising
    Sr = 0 + 1j * 0
    # l = 2
    Weyl4_l2_m0 = 0 + 1j * 0
    # positive m
    Weyl4_l2_m1 = 0 + 1j * 0
    Weyl4_l2_m2 = 0 + 1j * 0
    # negative m
    Weyl4_l2_m1n = 0 + 1j * 0
    Weyl4_l2_m2n = 0 + 1j * 0

    # l = 3
    Weyl4_l3_m0 = 0 + 1j * 0
    # positive m
    Weyl4_l3_m1 = 0 + 1j * 0
    Weyl4_l3_m2 = 0 + 1j * 0
    Weyl4_l3_m3 = 0 + 1j * 0

    # negative m
    Weyl4_l3_m1n = 0 + 1j * 0
    Weyl4_l3_m2n = 0 + 1j * 0
    Weyl4_l3_m3n = 0 + 1j * 0

    # l = 4
    Weyl4_l4_m0 = 0 + 1j * 0
    # positive m
    Weyl4_l4_m1 = 0 + 1j * 0
    Weyl4_l4_m2 = 0 + 1j * 0
    Weyl4_l4_m3 = 0 + 1j * 0
    Weyl4_l4_m4 = 0 + 1j * 0
    # negative m
    Weyl4_l4_m1n = 0 + 1j * 0
    Weyl4_l4_m2n = 0 + 1j * 0
    Weyl4_l4_m3n = 0 + 1j * 0
    Weyl4_l4_m4n = 0 + 1j * 0

    # k is a counter
    for k in range(phi_length):

        phi_var = phi[k]
        theta_var = theta[k]
        x1 = extraction_radius * np.cos(phi_var) * np.sin(theta_var) + float(center[0])
        y1 = extraction_radius * np.sin(phi_var) * np.sin(theta_var) + float(center[1])
        z1 = extraction_radius * np.cos(theta_var) + float(center[2])
        c = [x1, y1, z1]
        ptn = i.point(c)
        ReWeyl = float(ptn["ReWeyl4"][0])
        ImWeyl = float(ptn["ImWeyl4"][0])
        Weyl4 = ReWeyl + 1j * ImWeyl

        Weyl4_l2_m0 += (
            4 * pi * w[k] * np.conjugate(sY_l2_m0[k]) * Weyl4 * extraction_radius
        )
        # positive m
        Weyl4_l2_m1 += (
            4 * pi * w[k] * np.conjugate(sY_l2_m1[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l2_m2 += (
            4 * pi * w[k] * np.conjugate(sY_l2_m2[k]) * Weyl4 * extraction_radius
        )
        # negative m
        Weyl4_l2_m1n += (
            4 * pi * w[k] * np.conjugate(sY_l2_m1n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l2_m2n += (
            4 * pi * w[k] * np.conjugate(sY_l2_m2n[k]) * Weyl4 * extraction_radius
        )
        # l = 3
        Weyl4_l3_m0 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m0[k]) * Weyl4 * extraction_radius
        )
        # positive m
        Weyl4_l3_m1 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m1[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m2 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m2[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m3 += (
            4 * pi * w[k] * np.conjugate(sY_l3_m3[k]) * Weyl4 * extraction_radius
        )
        # negative m
        Weyl4_l3_m1n += (
            4 * pi * w[k] * np.conjugate(sY_l3_m1n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m2n += (
            4 * pi * w[k] * np.conjugate(sY_l3_m2n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l3_m3n += (
            4 * pi * w[k] * np.conjugate(sY_l3_m3n[k]) * Weyl4 * extraction_radius
        )
        # l = 4
        Weyl4_l4_m0 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m0[k]) * Weyl4 * extraction_radius
        )
        # positive m
        Weyl4_l4_m1 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m1[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m2 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m2[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m3 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m3[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m4 += (
            4 * pi * w[k] * np.conjugate(sY_l4_m4[k]) * Weyl4 * extraction_radius
        )
        # negative m
        Weyl4_l4_m1n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m1n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m2n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m2n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m3n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m3n[k]) * Weyl4 * extraction_radius
        )
        Weyl4_l4_m4n += (
            4 * pi * w[k] * np.conjugate(sY_l4_m4n[k]) * Weyl4 * extraction_radius
        )
    # Data Writeout
    # l = 2
    Weyl4_l2_m0_data.append(Weyl4_l2_m0)
    # positive m
    Weyl4_l2_m1_data.append(Weyl4_l2_m1)
    Weyl4_l2_m2_data.append(Weyl4_l2_m2)
    # negative m
    Weyl4_l2_m1n_data.append(Weyl4_l2_m1n)
    Weyl4_l2_m2n_data.append(Weyl4_l2_m2n)

    # l = 3
    Weyl4_l3_m0_data.append(Weyl4_l3_m0)
    # positive m
    Weyl4_l3_m1_data.append(Weyl4_l3_m1)
    Weyl4_l3_m2_data.append(Weyl4_l3_m2)
    Weyl4_l3_m3_data.append(Weyl4_l3_m3)
    # negative m
    Weyl4_l3_m1n_data.append(Weyl4_l3_m1n)
    Weyl4_l3_m2n_data.append(Weyl4_l3_m2n)
    Weyl4_l3_m3n_data.append(Weyl4_l3_m3n)

    # l = 4
    Weyl4_l4_m0_data.append(Weyl4_l4_m0)
    # positive m
    Weyl4_l4_m1_data.append(Weyl4_l4_m1)
    Weyl4_l4_m2_data.append(Weyl4_l4_m2)
    Weyl4_l4_m3_data.append(Weyl4_l4_m3)
    Weyl4_l4_m4_data.append(Weyl4_l4_m4)
    # negative m
    Weyl4_l4_m1n_data.append(Weyl4_l4_m1n)
    Weyl4_l4_m2n_data.append(Weyl4_l4_m2n)
    Weyl4_l4_m3n_data.append(Weyl4_l4_m3n)
    Weyl4_l4_m4n_data.append(Weyl4_l4_m4n)

    # time
    timedata.append(i.current_time)

    np.savetxt("time.out", timedata)
    # l = 2
    np.savetxt("Weyl4_l2_m0_data.out", Weyl4_l2_m0_data)
    np.savetxt("Weyl4_l2_m1_data.out", Weyl4_l2_m1_data)
    np.savetxt("Weyl4_l2_m2_data.out", Weyl4_l2_m2_data)
    np.savetxt("Weyl4_l2_m1n_data.out", Weyl4_l2_m1n_data)
    np.savetxt("Weyl4_l2_m2n_data.out", Weyl4_l2_m2n_data)
    # l = 3
    np.savetxt("Weyl4_l3_m0_data.out", Weyl4_l3_m0_data)
    np.savetxt("Weyl4_l3_m1_data.out", Weyl4_l3_m1_data)
    np.savetxt("Weyl4_l3_m2_data.out", Weyl4_l3_m2_data)
    np.savetxt("Weyl4_l3_m3_data.out", Weyl4_l3_m3_data)
    np.savetxt("Weyl4_l3_m1n_data.out", Weyl4_l3_m1n_data)
    np.savetxt("Weyl4_l3_m2n_data.out", Weyl4_l3_m2n_data)
    np.savetxt("Weyl4_l3_m3n_data.out", Weyl4_l3_m3n_data)
    # l = 4
    np.savetxt("Weyl4_l4_m0_data.out", Weyl4_l4_m0_data)
    np.savetxt("Weyl4_l4_m1_data.out", Weyl4_l4_m1_data)
    np.savetxt("Weyl4_l4_m2_data.out", Weyl4_l4_m2_data)
    np.savetxt("Weyl4_l4_m3_data.out", Weyl4_l4_m3_data)
    np.savetxt("Weyl4_l4_m4_data.out", Weyl4_l4_m4_data)
    np.savetxt("Weyl4_l4_m1n_data.out", Weyl4_l4_m1n_data)
    np.savetxt("Weyl4_l4_m2n_data.out", Weyl4_l4_m2n_data)
    np.savetxt("Weyl4_l4_m3n_data.out", Weyl4_l4_m3n_data)
    np.savetxt("Weyl4_l4_m4n_data.out", Weyl4_l4_m4n_data)

    program_runtime_data.append(time.time() - i_start)
    np.savetxt("run_time.out", program_runtime_data)
