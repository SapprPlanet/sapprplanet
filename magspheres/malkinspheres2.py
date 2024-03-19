import time
import os
import sys
import math
import pandas
import matplotlib.pyplot
import mpl_toolkits.mplot3d
import numpy
from scipy import spatial
from matplotlib import cm

rootdir = 'maven.mag.calibrated'
# rootdir = 'test'


def set_axes_equal(ax: matplotlib.pyplot.Axes):
    limits = numpy.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    origin = numpy.mean(limits, axis=1)
    radius = 0.5 * numpy.max(numpy.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)


def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])


def make_views(ax, angles, elevation=None, width=16, height=9, prefix='f_'):
    files = []
    ax.figure.set_size_inches(width, height)
    for i, angle in enumerate(angles):
        ax.view_init(elev=elevation, azim=angle)
        fname = '%s%03d.jpeg' % (prefix, i)
        ax.figure.savefig(fname)
        files.append(fname)
    return files


def make_movie(files, output, fps=60, bitrate=3500):
    output_name, output_ext = os.path.splitext(output)
    command = {'.mp4': 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc x264\
                         -x264encopts subq=6:partitions=all:8x8dct:me=umh:frameref=5:bframes=3:b_pyramid=normal:weight_b:bitrate=%d'
               % (",".join(files), fps, output_name, bitrate)}
    print(command[output_ext])
    output_ext = os.path.splitext(output)[1]
    os.system(command[output_ext])


def rotanimate(ax, angles, output):
    output_ext = os.path.splitext(output)[1]
    files = make_views(ax, angles)
    make_movie(files, output)
    for f in files:
        os.remove(f)


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

# bu, bl = [], []
# bu.append(numpy.pi / 2)
# for i in range(Nhalfring):
#     bl.append(numpy.arcsin(numpy.sin(bu[i]) - A / dL[i]))
#     bu.append(bl[i])
# bl[i] = 0.0
# bu[i+1] = bl[i]
# for i in range(1, Nhalfring+1):
#     bl.append(0 - bu[Nhalfring-i])
#     bu.append(0 - bl[Nhalfring-i-1])
# for i in range (Nring):
#  	  print(str(bu[i] * 180 / numpy.pi) + "\t" + str(bl[i] * 180 /numpy.pi))

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

if sys.argv[3] != "draw":
    Nfilled = 0
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
                    # Bx = float(row[7])
                    # By = float(row[8])
                    # Bz = float(row[9])

                    # Spacecraft position X-component in Planetocentric coordinate system
                    # posX = float(row[11])
                    # Spacecraft position Y-component in Planetocentric coordinate system
                    # posY = float(row[12])
                    # Spacecraft position Z-component in Planetocentric coordinate system
                    # posZ = float(row[13])
                    pos = numpy.array(
                        [float(row[11]), float(row[12]), float(row[13])])
                    # space_track = numpy.append(space_track, [pos], axis=0)

                    # Sr = math.sqrt(posX * posX + posY * posY + posZ * posZ)
                    # St = numpy.arctan2(numpy.sqrt(posX * posX + posY * posY), posZ)
                    # Sp = math.atan2(posY, posX)

                    # posX1 = math.sin(St) * math.cos(Sp)
                    # posY1 = math.sin(St) * math.sin(Sp)
                    # posZ1 = math.cos(St)

                    # positionXYZnorm = numpy.array(
                    #     (numpy.sin(St) * numpy.cos(Sp), numpy.sin(St) * numpy.sin(Sp), math.cos(St)))

                    # print(str(St) + " " + str(Sp))

                    # Find Nearest Sector
                    delta, indexNearest = grid_tree.query(pos)
                    if cell_filled[indexNearest] == False:
                        cell_filled[indexNearest] = True
                        Nfilled += 1
                        grid_track = numpy.append(
                            grid_track, [grid[indexNearest]], axis=0)
                    if delta < cell_delta[indexNearest]:
                        cell_B[indexNearest] = [
                            float(row[7]), float(row[8]), float(row[9])]
                        cell_XYZ[indexNearest] = [
                            float(row[11]), float(row[12]), float(row[13])]
                        cell_delta[indexNearest] = delta

                    # One Track
                    # x2.append(posX)
                    # y2.append(posY)
                    # z2.append(posZ)

                    # One Track on Sphere
                    # x3.append(math.sin(St) * math.cos(Sp))
                    # y3.append(math.sin(St) * math.sin(Sp))
                    # z3.append(math.cos(St))

                    # One Track on Grid
                    # x4.append(x[indexNearest])
                    # y4.append(y[indexNearest])
                    # z4.append(z[indexNearest])

                    # print(str(Bx) + "\t" + str(By) + "\t" + str(Bz) + "\t" + str(x) + "\t" + str(y) + "\t" + str(z))

    end = time.time()
    print(end-start)
    print("Nfilled = " + str(Nfilled) +
          " (" + str(int(100*Nfilled/Ncell)) + " %)")

    outF = open("out.txt", "w")
    outF.write("Ncell = " + str(Ncell) + "\n")
    outF.write("Nfilled = " + str(Nfilled) + "\n")
    for i in range(Ncell):
        outF.write(str(i+1) + "\t" +
                   str(cell_B[i][0]) + "\t" + str(cell_B[i][1]) + "\t" + str(cell_B[i][2]) + "\t" +
                   str(cell_XYZ[i][0]) + "\t" + str(cell_XYZ[i][1]) + "\t" + str(cell_XYZ[i][2]) + "\n")
    outF.close()
    print("File saved")

else:
    print("Drawing...")
    outF = open("testdata.txt", "w")
    csv_reader = pandas.read_csv(sys.argv[2], sep='\s+', skiprows=2)
    i = 0
    j = 0
    cell_R = []
    for row in csv_reader.values:
        if ring[i] > 14 and ring[i] < 165:
            cell_B[i] = [float(row[1]), float(row[2]), float(row[3])]
            cell_XYZ[i] = [float(row[4]), float(row[5]), float(row[6])]
            Btotal = numpy.sqrt(float(row[1])**2 +
                                float(row[2])**2 + float(row[3])**2)
            Sr = numpy.sqrt(float(row[4])**2 +
                            float(row[5])**2 + float(row[6])**2)
            St = numpy.arctan2(float(row[5]), float(row[4]))
            Sp = numpy.arctan2(numpy.sqrt(
                float(row[4])**2 + float(row[5])**2), float(row[6]))
            if numpy.linalg.norm(cell_XYZ[i]) > 0:
                cell_R.append(Sr)
                outF.write(str(j+1) + "\t" +
                           str(Btotal) + "\t" + str(Sr) + "\t" + str(St*180/numpy.pi) + "\t" +
                           str(Sp*180/numpy.pi) + "\n")
                j += 1
        i += 1
    print("\nmedian of cell_R : ", numpy.median(cell_R))
    print("\nmin of cell_R : ", numpy.min(cell_R))
    print("\nmax of cell_R : ", numpy.max(cell_R))
    print("\naverage of cell_R : ", numpy.average(cell_R))
    outF.close()

    # One Track on Grid
    # for i in range(Ncell):
    #    if cell_filled[i] == True:
    #        x4.append(x[i])
    #        y4.append(y[i])
    #        z4.append(z[i])

    fig = matplotlib.pyplot.figure()
    ax = fig.gca(projection='3d')

    # Generate and plot a unit sphere
    # u = numpy.linspace(0, 2*numpy.pi, 100)
    # v = numpy.linspace(0, numpy.pi, 100)
    # x = numpy.outer(numpy.cos(u), numpy.sin(v))
    # y = numpy.outer(numpy.sin(u), numpy.sin(v))
    # z = numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))
    # ax.plot_surface(x*3393.5, y*3393.5, z*3393.5, alpha=0.3, color="#ad6242")

    # Render
    fig.suptitle("Equal-area grid on a spherical surface\n(41252 cells)")
    ax.scatter(*grid.T, color="k", s=1)

    # fig.suptitle("Spacecraft trajectory\n(October 10, 2014)")
    # ax.scatter(*space_track.T, color="c", s=1)

    # fig.suptitle(
    #     "Vertical projection of the spacecraft trajectory onto the grid\n(2015)")
    # ax.scatter(*grid_track.T, color="r", s=1)

    fig.suptitle(
        "All-time minimum spacecraft position on the grid\n(October 10, 2014 - October 31, 2020)")
    cell_XYZ = cell_XYZ[~numpy.all(cell_XYZ == 0, axis=1)]
    # ax.scatter(*cell_XYZ.T, s=1, c="c")
    R = numpy.linalg.norm(cell_XYZ, axis=1)
    psm = ax.scatter(*cell_XYZ.T, s=1, c=R, cmap=cm.coolwarm)
    fig.colorbar(psm, ax=ax)

    ax.set_box_aspect([1, 1, 1])
    set_axes_equal(ax)

    # Animation
    # angles = numpy.linspace(0, 360, 600)[:-1]
    # rotanimate(ax, angles, 'movie.mp4')

    matplotlib.pyplot.show()
