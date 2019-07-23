import argparse
import ast
import sys
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
                             'data is required (see -p/--plink). Lines beginning with a number sign '
                             '(#) are ignored.')
data_input.add_argument('-p', '--plink',
                        help='Genetic data in PLINK format (.bam, .fam, .bim) without the extension (if files'
                             ' are DATA.bam, DATA.fam, and DATA.bim, use "--plink DATA" ). Either non-genetic '
                             'data or PLINK-formatted genetic data is required (see -n/--non-genetic)')
required_arguments.add_argument('-c', '--coords', required=True,
                                help='TSV (tab-separated values) or similar file with the locations of each '
                                     'sample. One row contains one sample\'s longitude and latitude. '
                                     'Lines beginning with a number sign (#) are ignored.')

# parsing optional arguments
parser.add_argument('-o', '--output', type=str, default="womble_out",
                    help='Output prefix (default: "womble_out")')
parser.add_argument('--coord-sep', type=str,
                    help='Separator for the positions-file (default: whitespace)')
parser.add_argument('--non-gen-sep', type=str,
                    help='Separator for the non-genetic data (default: whitespace)')
parser.add_argument('-d', '--density', default=(1, 1), type=ast.literal_eval,
                    help='Density (per unit longitude and latitude) at which to compute wombling '
                         '(default: (1, 1) )')
parser.add_argument('--grid-bounds', type=ast.literal_eval,
                    help='Boundary of the grid in which wombling is computed, in longitude and latitude '
                         '(minLong, maxLong, minLat, maxLat, e.g. "(120, 125, -6, 3)"). If this is not set, '
                         'the default boundaries are set using the most extreme coordinates in COORDS.')
parser.add_argument('-i', '--identical', default=False, action='store_true',
                    help='If selected, the traits from samples at a single location are not averaged before '
                         'interpolation. Interpolation then treats a single location with two samples as '
                         'exerting twice as much influence. Only relevant when there are duplicate locations.')
parser.add_argument('-s', '--suppress-save', default=False, action='store_true',
                    help='By default, the script will save the grid\'s latitudes and longitudes, along with '
                         'the calculated rate-of-change, direction-of-change, and wighted direction-of-change. '
                         'This option suppresses these from being output.')
parser.add_argument('--percentile', default=5., type=float,
                    help='Only plot the grid-squares of the top N-th percentile (default: 5. The top 5th '
                         'percentile will select only the top 1 in 20 tiles by absolute rate of change).')
parser.add_argument('--no-plot', default=False, action='store_true',
                    help='Disable plotting and its output.')
parser.add_argument('--no-plot-save', default=False, action='store_true',
                    help='If plotting, only display (do not save) the plot.')
parser.add_argument('--alpha', default=0.5, type=float,
                    help='Alpha value (opacity) for rate-of-change heat-map (Default: 0.5)')

args = parser.parse_args()

# load library for processing PLINK files
if args.plink:
    from plinkio import plinkfile

# load libraries for plotting
if not args.no_plot:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from scipy.ndimage import label

# importing data
sampleLoc = np.loadtxt(args.coords, delimiter=args.coord_sep)
if np.shape(sampleLoc)[1] != 2:
    raise IndexError('{} has {} separate columns, but ought to be 2 '
                     '(check COORDS_SEP?)'.format(args.coords, np.shape(sampleLoc)[1]))

if args.plink:
    sampleData = plinkfile.open(args.plink)
    print('Using data file (PLINK):  {}'.format(args.plink), file=sys.stderr)
    printEveryNthLine = 10000
else:
    sampleData = np.transpose(np.loadtxt(args.non_genetic, delimiter=args.non_gen_sep))
    if np.shape(sampleData)[1] != np.shape(sampleLoc)[0]:
        raise IndexError('Make sure that {} ({} rows) has one row per position '
                         'in {} ({} positions)!'.format(args.non_genetic, np.shape(sampleData)[1],
                                                        args.coords, np.shape(sampleLoc)[0]))
    print('Using data file:  {}'.format(args.non_genetic), file=sys.stderr)
    printEveryNthLine = int(np.shape(sampleData)[0]/10)

if args.density:
    if not (type(args.density) is tuple):
        raise TypeError('DENSITY must be a tuple of size 2 (default: "(1, 1)", input: {} )'.format(args.density))
    elif len(args.density) != 2:
        raise IndexError('DENSITY must be a tuple of size 2 (default: "(1, 1)", input: {} )'.format(args.density))
    else:
        print('Using density:  {}'.format(args.density), file=sys.stderr)

if args.grid_bounds:
    # If grid_bounds are set, make sure that they are properly formatted
    grid_bounds = args.grid_bounds
    if not (type(grid_bounds) is tuple):
        raise TypeError('GRID_BOUNDS must be a tuple with 4 elements (e.g. "(120, 125, -6, 3)"), '
                        'but input was: {}, of type {}'.format(args.grid_bounds,
                                                               type(args.grid_bounds)))
    elif len(grid_bounds) != 4:
        raise IndexError('GRID_BOUNDS must be a tuple with 4 elements (e.g. "(120, 125, -6, 3)"), '
                         'but input was: {}, of length {}'.format(args.grid_bounds,
                                                                  len(args.grid_bounds)))
    if args.grid_bounds[0] >= args.grid_bounds[1] or args.grid_bounds[2] >= args.grid_bounds[3]:
        print('WARNING: Are you sure that GRID_BOUNDS is formatted correctly? \n'
              '         (minLong, maxLong, minLat, maxLat, e.g. "(120, 125, -6, 3)")', file=sys.stderr)
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

# Create grid for interpolation
longs = np.arange(grid_bounds[0], grid_bounds[1] + 1/ args.density[0], 1 / args.density[0])
lats = np.arange(grid_bounds[2], grid_bounds[3] + 1/ args.density[1], 1 / args.density[1])
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
        dist_matrix[i, j] = (distances ** -3) * average_multiplier
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
    if rowNum % printEveryNthLine == 1:
        print("Row: {}".format(rowNum), file=sys.stderr)

    if len(np.unique(row)) == 1:
        continue  # skip positions without variation

    # find the mean genotype per location
    mean_trait_value = []

    for location in range(uniqueSampleLoc.shape[0]):
        mean_trait_value.append(np.mean(np.array(row)[uniqueSampleIndices == location]))

    interpolated = dist_matrix @ np.array(mean_trait_value)

    # Calculate rates and directions of change per variable at the midpoints of the grid-squares.
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
        doubledAngle_dY += np.sin(midpoints_angles) * rateOfChange
        doubledAngle_dX += np.cos(midpoints_angles) * rateOfChange

    doubledDirectionOfChange += midpoints_angles
print("All rows processed", file=sys.stderr)

finalDirectionOfChange = doubledDirectionOfChange / 2
with np.errstate(all='ignore'):
    finalWeightedDirectionOfChange = np.arctan(np.divide(doubledAngle_dY, doubledAngle_dX))/2
    finalWeightedDirectionOfChange[doubledAngle_dX < 0] += np.pi/2

print("Wombling completed.", file=sys.stderr)
# End of wombling

# Reporting results
if not args.suppress_save:
    print('Saving results (rate of change, direction of change) and info needed for plotting \n'
          '  (longitudes and latitudes of the grid).', file=sys.stderr)
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
    print("Drawing plots...", file=sys.stderr)
    # Using a percentile, select regions with the top N-th percentile rates-of-change. These are the core regions.
    #    The 2nd N-th percentile is to be used as a "linked"; they are only to be plotted if they are connected to
    #    regions in the top N-th percentile. Plot quivers for only the core regions and connected regions in the
    #    2nd N-th percentile.
    # 1. Find 1st and 2nd N-th percentile cutoffs
    gridSize = len(rateOfChange.ravel().tolist())
    if args.percentile > 50:
        perc = 50
    elif int(args.percentile * gridSize / 100) < 1:
        perc = 100 / gridSize
    else:
        perc = args.percentile
    cutoff_1, cutoff_2 = np.array([value[0] for value in
                                   sorted(zip((rateOfChange.ravel().tolist()),
                                              range(gridSize)),
                                          reverse=True)])[[int(perc * gridSize/100),
                                                           int(perc * 2 * gridSize/100)]]

    # 2. Use rate-of-change matrix to find regions within the top first 2 N-th percentiles.
    #    These are connected and labeled.
    grid_labeled, num_features = label(rateOfChange > cutoff_2,
                                       structure=[[1, 1, 1],
                                                  [1, 1, 1],
                                                  [1, 1, 1]])
    for label_ in np.unique(grid_labeled):
        # For any label which does not encompass a region in the first N-th percentile,
        if label_ not in np.unique(grid_labeled[rateOfChange > cutoff_1]):
            grid_labeled[grid_labeled == label_] = 0

    m = Basemap(projection='cyl',
                llcrnrlat=grid_bounds[2],
                urcrnrlat=grid_bounds[3],
                llcrnrlon=grid_bounds[0],
                urcrnrlon=grid_bounds[1],
                resolution='i')
    # m.drawcoastlines(linewidth=0.5)
    m.drawmapboundary(fill_color='#bfffff', zorder=0)
    m.fillcontinents(color='#ffffff', lake_color='#bfffff', zorder=1)

    plt.tight_layout(pad=0, h_pad=0, w_pad=0)

    # Heat-map of rates-of-change, with midpoints of grid-squares being the points for which
    # rate-of-change was calculated.
    plt.imshow(np.fliplr(rateOfChange).T,
               extent=[longs.min(), longs.max(),
                       lats.min(), lats.max()],
               cmap="Reds", alpha=args.alpha, zorder=2)

    plt.quiver(mid_longs[grid_labeled > 0], mid_lats[grid_labeled > 0],  # Longitude, Latitude
               np.cos(finalWeightedDirectionOfChange[grid_labeled > 0]),
               np.sin(finalWeightedDirectionOfChange[grid_labeled > 0]),
               headwidth=0, headlength=0, color="Black",
               headaxislength=0, width=0.005, pivot='mid', zorder=3)

    if args.no_plot_save:
        plt.show()
    else:
        plt.savefig(fname="{}_{}_{}.png".format(args.output, *args.density), dpi=300, format="png")

    plt.clf()

print(asctime(), "\nDone\n", file=sys.stderr)
