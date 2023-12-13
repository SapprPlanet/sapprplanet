import time
import os
import sys
import math
import pandas
import numpy
from scipy import spatial

rootdir = 'maven.mag.calibrated'

# The given number of half rings, Nhalfring
Nhalfring = int(sys.argv[1])

# The given number of rings, Nring
Nring = Nhalfring * 2
print("Nring = " + str(Nring))

# A spherical surface is first split into several latitudinal rings (bands)
# of constant width, dB
dB = numpy.pi / Nring
print("dB = " + str(dB * 180 / numpy.pi) + "\N{DEGREE SIGN}")

b0, dL, dN, dI = [], [], [], []
Ncell = 0
# The ring number, i
for i in range(Nring):
    # The central latitude of the ring, b0[i]
    b0.append(dB * i + dB / 2)
    # print("b0[" + str(i) + "] = " + str(b0[i] * 180 / numpy.pi) + "\N{DEGREE SIGN}")

    # The longitudinal span of the cells in each ring, dL[i]
    dL.append(dB / numpy.sin(b0[i]))
    # print("dL[" + str(i) + "] = " + str(dL[i] * 180 / numpy.pi) + "\N{DEGREE SIGN}")

    # The number of cells in each ring, dN[i]
    dN.append(round(2 * numpy.pi / dL[i]))
    # The first index of cells in each ring, dI[i]
    if i > 0:
        dI.append(dI[i-1] + dN[i-1])
    else:
        dI.append(0)
    # print(str(dI[i]))
    Ncell += dN[i]

# The number of cells in the grid, Ncell
print("Ncell = " + str(Ncell))

# The cell area, A
A = 4.0 * numpy.pi / Ncell
print("A = " + str(A * 180 * 180 / numpy.pi / numpy.pi) +
      "\N{DEGREE SIGN}\N{DEGREE SIGN}")

ring = []
grid = numpy.empty([0, 3])

cell_filled = []
cell_B = numpy.empty([0, 3])
cell_XYZ = numpy.empty([0, 3])
cell_delta = []

k = 0
for i in range(Nring):
    theta = b0[i]
    for j in range(dN[i]):
        phi = j * 2 * numpy.pi / dN[i] + numpy.pi / dN[i]
        grid = numpy.append(
            grid, [[numpy.cos(phi) * numpy.sin(theta), numpy.sin(phi) * numpy.sin(theta), numpy.cos(theta)]], axis=0)
        ring.append(i)
        cell_filled.append(False)
        cell_B = numpy.append(cell_B, [[0., 0., 0.]], axis=0)
        cell_XYZ = numpy.append(cell_XYZ, [[0., 0., 0.]], axis=0)
        cell_delta.append(sys.float_info.max)
        k += 1
grid_tree = spatial.cKDTree(grid)

Total = 0
grid_track = numpy.empty([0, 3])
start = time.time()
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if file.endswith('.sts'):
            filename = os.path.join(subdir, file)

            # Find number of skipping rows (skip)
            stsfile = open(filename, 'r')
            skip = 1
            result = 0
            while result < 4:
                skip += 1
                line = stsfile.readline()
                if line.find('END_OBJECT') < 0:
                    result = 0
                else:
                    result += 1
            stsfile.close()

            print(filename)
            csv_reader = pandas.read_csv(
                filename, sep='\s+', skiprows=skip)
            space_track = numpy.empty([0, 3])
            for row in csv_reader.values:
                Total += 1

end = time.time()
print(end-start)
print("Total = " + str(Total) + "\n")

