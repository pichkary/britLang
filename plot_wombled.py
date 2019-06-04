import argparse, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from skimage.morphology import label

parser = argparse.ArgumentParser(prog='Womble plotter',
                                 description='Plot saved data from womble.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

required_arguments = parser.add_argument_group('required arguments')
required_arguments.add_argument('-r', '--rates', required=True, type=str,
                                help='Rates of change generated by womble.py for each point on the grid.')
required_arguments.add_argument('-d', '--directions', required=True, type=str,
                                help='Directions of change generated by womble.py for every point of the grid.'
                                     ' (weighted direction of change is preferred.')
required_arguments.add_argument('-g', '--longitudes', required=True, type=str,
                                help='Longitudes of the grid generated by womble.py')
required_arguments.add_argument('-t', '--latitudes', required=True, type=str,
                                help='Latitudes of the grid generated by womble.py')
parser.add_argument('-p', '--percentile', required=False, default=5., type=float,
                    help='Only plot the grid-squares of the top N-th percentile (default: 5, meaning that '
                    'the top 1 in 20 tiles (top by absolute rate of change) will be shown).')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.0,
                    help='The opacity of the heatmap of rates of change (Default: 0; 0 <= alpha <= 1).')

args = parser.parse_args()

lats = np.loadtxt(args.latitudes)
longs = np.loadtxt(args.longitudes)
rateOfChange = np.loadtxt(args.rates)
directionOfChange= np.loadtxt(args.directions)
mid_lats, mid_longs = np.meshgrid(lats, longs)

# How this ought to work:
# Using a percentile, select regions with the top N-th percentile rates-of-change. These are the core regions.
#    The 2nd N-th percentile is to be used as a "linked"; they are only to be plotted if they are connected to
#    regions in the top N-th percentile. Plot quivers for only the core regions and connected regions in the
#    2nd N-th percentile.
# 1. Find 1st and 2nd N-th percentile cutoffs
cutoff_1, cutoff_2 = np.array([value[0] for value in
                               sorted(zip((rateOfChange.ravel().tolist()),
                                          range(len(rateOfChange.ravel().tolist()))),
                                      reverse=True)])[[int(args.percentile),
                                                       int(args.percentile*2)]]

# 2. Use rate-of-change matrix to find regions within the top first 2 N-th percentiles.

grid_labeled = label(rateOfChange > cutoff_2, connectivity=2)
for label_ in np.unique(grid_labeled):
    if label_ not in np.unique(grid_labeled[rateOfChange > cutoff_1]):
        grid_labeled[grid_labeled == label_] = 0

m = Basemap(projection='cyl',
            llcrnrlat=min(lats),
            urcrnrlat=max(lats),
            llcrnrlon=min(longs),
            urcrnrlon=max(longs),
            resolution='i')
m.drawcoastlines(linewidth=0.5)

plt.tight_layout(pad=0, h_pad=0, w_pad=0)

plt.imshow(np.fliplr(rateOfChange).T,
           extent=[longs.min(), longs.max(),
                   lats.min(), lats.max()],
           cmap="Reds", alpha=args.alpha)

plt.quiver(mid_longs[grid_labeled>0], mid_lats[grid_labeled>0],  # Longitude, Latitude
           np.cos(directionOfChange[grid_labeled>0]),
           np.sin(directionOfChange[grid_labeled>0]),
           headwidth=0, headlength=0, color="Black",
           headaxislength=0, width=0.0005, pivot='mid')

plt.show()
