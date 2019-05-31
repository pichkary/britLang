import argparse, sys, ast
from time import asctime
import numpy as np
from cHaversine import haversine

# parsing
parser = argparse.ArgumentParser(prog='Womble Script 0.1',
                                 description='Womble genetic or non-genetic data',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# parsing required arguments
required_arguments = parser.add_argument_group('required arguments')
data_input = required_arguments.add_mutually_exclusive_group(required=True)
data_input.add_argument('-n', '--non-genetic',
                        help='TSV (tab-separated values) or similar file containing non-genenetic data, '
                             'one sample per row. Either non-genetic data or PLINK-formatted genetic '
                             'data is required (see -p/--plink). May have non-data lines beginning with a #')
data_input.add_argument('-p', '--plink',
                        help='Genetic data in PLINK format (.bam, .fam, .bim) without the extension (if files'
                             ' are DATA.bam, DATA.fam, and DATA.bim, use "--plink DATA" ). Either non-genetic '
                             'data or PLINK-formatted genetic data is required (see -n/--non-genetic)')
required_arguments.add_argument('-c', '--coords', required=True,
                                help='TSV (tab-separated values) or similar file with the locations of each '
                                     'sample. One sample per row and longitude and latitude for each '
                                     'sample. May have non-data lines beginning with a #')

# parsing optional arguments
parser.add_argument('-o', '--output', type=str, default="womble_out",
                    help='Output prefix (default: "womble_out")')
parser.add_argument('--coord-sep', type=str,
                    help='Separator for the positions-file (default: whitespace)')
parser.add_argument('--non-gen-sep', type=str,
                    help='Separator for the non-genetic data (default: whitespace)')
parser.add_argument('--density', default=(1, 1), type=ast.literal_eval,
                    help='Density (per unit longitude and latitude) at which to compute wombling '
                         '(default: (1, 1) )')
parser.add_argument('--grid-bounds', type=ast.literal_eval,
                    help='Boundary of the grid in which wombling is computed, in longitude and latitude '
                         '(minLong, maxLong, minLat, maxLat, e.g. (120, 125, -6, 3). If this is not set, '
                         'the default boundaries are 1.5x DENSITY away from the most extreme coordinates.')
parser.add_argument('-i', '--identical', default=False, action='store_true',
                    help='If selected, the traits from samples at a single location are not averaged before '
                         'interpolation. Interpolation then treats a single location with two samples as '
                         'exerting twice as much influence. Only relevant when there are duplicate locations.')
parser.add_argument('-s', '--suppress-save', default=False, action='store_true',
                    help='By default, the script will save the grid\'s latitudes and longitudes, along with '
                         'the calculated rate-of-change, direction-of-change, and wighted direction-of-change. '
                         'This option suppresses these from being output.')
parser.add_argument('--percentiles', default=5., type=float,  # TODO percentiles
                    help='[WIP] Only plot the grid-squares of the top N-th percentile (default: 5, meaning that'\
                         'only the top 1 in 20 tiles, by absolute rate of change, will be shown).')
parser.add_argument('--no-plot', default=False, action='store_true',
                    help='Disable plotting and its output.')

args = parser.parse_args()

# load library for processing PLINK files
if args.plink:
    from plinkio import plinkfile

# load libraries for plotting
if not args.no_plot:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap


# importing data
sampleLoc = np.loadtxt(args.coords, delimiter=args.coord_sep)
if np.shape(sampleLoc)[1] != 2:
    raise IndexError("COORDS is not made up of 2 separate columns (check COORDS_SEP?)")


if args.grid_bounds:
    grid_bounds = args.grid_bounds
else:
    # If grid_bounds were not set from the console, use the data to create them.
    grid_bounds = (np.min(sampleLoc[:, 0]) - 1.5 * args.density[0] ** -1,  # Longitude
                   np.max(sampleLoc[:, 0]) + 1.5 * args.density[0] ** -1,
                   np.min(sampleLoc[:, 1]) - 1.5 * args.density[1] ** -1,  # Latitude
                   np.max(sampleLoc[:, 1]) + 1.5 * args.density[1] ** -1)

print("Wombling grid bounds: \n" 
      "  Min. longitude: {: 3.4f} \n"
      "  Max longitude:  {: 3.4f} \n"
      "  Min. latitude:  {: 3.4f} \n"
      "  Max latitude:   {: 3.4f} \n".format(*grid_bounds),
      file=sys.stderr)

if args.plink:
    sampleData = plinkfile.open(args.plink)
else:
    sampleData = np.transpose(np.loadtxt(args.non_genetic, delimiter=args.non_gen_sep))


# # TODO percentiles
# PERCENTILE = 25  # 95 == only show top 5 percentile of tiles

# Create grid for interpolation
longs = np.arange(grid_bounds[0], grid_bounds[1] + 1, 1 / args.density[0])
lats = np.arange(grid_bounds[2], grid_bounds[3] + 1, 1 / args.density[1])
grid_lats, grid_longs = np.meshgrid(lats, longs)
mid_lats, mid_longs = np.meshgrid((lats[1:] + lats[:-1])/2, (longs[1:] + longs[:-1])/2)

dist_matrix = np.zeros((len(longs),
                        len(lats),
                        len(np.unique(sampleLoc, axis=0))
                        ))

# For each point on ther interpolation-grid, find the distance to each sample
uniqueSampleLoc, uniqueSampleInverse, uniqueSampleIndices, uniqueSampleCount =\
    np.unique(sampleLoc, return_inverse=True, return_index=True, return_counts=True, axis=0)
if not args.identical:
    average_multiplier = np.copy(uniqueSampleCount)
else:
    average_multiplier = np.ones_like(uniqueSampleCount)

for i in range(len(longs)):
    for j in range(len(lats)):
        distances = []
        # weighted interpolation
        for location in range(uniqueSampleLoc.shape[0]):
            distances.append(haversine(tuple(uniqueSampleLoc[location, :]),
                                       tuple([longs[i], lats[j]])))

        distances = np.array(distances)
        dist_matrix[i, j] = (distances ** -2) * average_multiplier
        dist_matrix[i, j] /= np.sum(dist_matrix[i, j])

# Initialize matrices to store rate-of-change, along with weighted and unweighted direction-of-change
rateOfChange = np.zeros((len(longs) - 1, len(lats) - 1))
doubledDirectionOfChange = np.copy(rateOfChange)
doubledAngle_dX = np.copy(rateOfChange)
doubledAngle_dY = np.copy(rateOfChange)

# The wombling itself
# (for details, see Barbujani, Oden, and Sokal, (1989)
# "Detecting regions of abrupt change in maps of biological variables"
# in Syst. Zool. 38(4) 376-389
# tl;dr: find rate and direction of change via wombling at the midpoint of each grid-square,
# one variable (SNP, linguistic element, etc.) at a time
rowNum = 0
for row in sampleData:
    rowNum += 1
    if rowNum % 10000 == 1:
        print("Row: {}".format(rowNum), file=sys.stderr)

    if len(np.unique(row)) == 1:
        continue  # skip positions without variation

    # find the mean genotype per location
    mean_trait_value = []
    
    for location in range(uniqueSampleLoc.shape[0]):
        mean_trait_value.append(np.mean(np.array(row)[uniqueSampleIndices == location]))

    interpolated = dist_matrix @ np.array(mean_trait_value)

    # Calculate rates and directions of change per genotype at the midpoints of the grid-squares.
    dY = (interpolated[1:, 1:] - interpolated[:-1, 1:] +
          0.5 * (interpolated[:-1, 1:] - interpolated[1:, 1:] +
                 interpolated[1:, :-1] - interpolated[:-1, :-1])
          ).astype(np.float64, copy=False)
    dX = (interpolated[:-1, :-1] - interpolated[:-1, 1:] +
          0.5 * (interpolated[:-1, 1:] - interpolated[1:, 1:] +
                 interpolated[1:, :-1] - interpolated[:-1, :-1])
          ).astype(np.float64, copy=False)

    rateOfChange += np.sqrt(dX ** 2 + dY ** 2)

    with np.errstate(all='ignore'):
        midpoints_angles = np.arctan(np.divide(dY, dX)) + np.pi/2
        midpoints_angles[np.isnan(midpoints_angles)] = 0
        midpoints_angles *= 2
        doubledAngle_dX += np.cos(midpoints_angles) * rateOfChange
        doubledAngle_dY += np.sin(midpoints_angles) * rateOfChange

    doubledDirectionOfChange += midpoints_angles

finalDirectionOfChange = doubledDirectionOfChange / 2
with np.errstate(all='ignore'):
    finalWeightedDirectionOfChange = np.arctan(np.divide(doubledAngle_dY, doubledAngle_dX))/2

# End of wombling

# Reporting results
if not args.suppress_save:
    np.savetxt("{}_longitudes_{}_{}.txt".format(args.output, *args.density),
               (longs[1:] + longs[:-1])/2, fmt="%.10f")
    np.savetxt("{}_latitudes_{}_{}.txt".format(args.output, *args.density),
               (lats[1:] + lats[:-1])/2, fmt="%.10f")
    np.savetxt("{}_rateOfChange_{}_{}.txt".format(args.output, *args.density),
               rateOfChange, fmt="%.10f")
    np.savetxt("{}_weightedDirOfChange_{}_{}.txt".format(args.output, *args.density),
               finalWeightedDirectionOfChange, fmt="%.10f")
    np.savetxt("{}_dirOfChange_{}_{}.txt".format(args.output, *args.density),
               finalDirectionOfChange, fmt="%.10f")

# Plotting rates of change
if not args.no_plot:
    m = Basemap(projection='cyl',
                llcrnrlat=grid_bounds[2],
                urcrnrlat=grid_bounds[3],
                llcrnrlon=grid_bounds[0],
                urcrnrlon=grid_bounds[1],
                resolution='i')
    m.drawcoastlines(linewidth=0.5)

    plt.tight_layout(pad=0, h_pad=0, w_pad=0)

    plt.imshow(np.fliplr(rateOfChange).T,
               extent=[longs.min(), longs.max(),
                       lats.min(), lats.max()],
               cmap="Reds", alpha=0.9)

    plt.scatter(sampleLoc[:, 0], sampleLoc[:, 1], c='green', alpha=0.3)

    plt.quiver(mid_longs, mid_lats,  # Longitude, Latitude
               np.cos(finalWeightedDirectionOfChange),
               np.sin(finalWeightedDirectionOfChange),
               headwidth=0, headlength=0, color="Black",
               headaxislength=0, width=0.0005, pivot='mid')

    plt.savefig(fname="{}_{}_{}.png".format(args.output, *args.density), dpi=300, format="png")

    plt.clf()

print(asctime(), file=sys.stderr)
