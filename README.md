# britLang
A tool to detect genetic and non-genetic barriers in geography, written in Python 3.7

Required python packages:
* numpy 
* matplotlib
* cHaversine
* basemap  (from matplotlib_toolkits)
* plinkio  (when using genetic data)
* sci-kit learn (for skimage labeling)

### Usage
  python3 womble.py [-h] (-n NON_GENETIC | -p PLINK) -c COORDS
                    [-o OUTPUT] [--coord-sep COORD_SEP]
                    [--non-gen-sep NON_GEN_SEP] [--density DENSITY]
                    [--grid-bounds GRID_BOUNDS] [-i] [-s]
                    [--percentiles PERCENTILES] [--no-plot]

required arguments:
  -n NON_GENETIC, --non-genetic NON_GENETIC
                        TSV (tab-separated values) or similar file containing
                        non-genenetic data, one sample per row. Either non-
                        genetic data or PLINK-formatted genetic data is
                        required (see -p/--plink). May have non-data lines
                        beginning with a #
  -p PLINK, --plink PLINK
                        Genetic data in PLINK format (.bam, .fam, .bim)
                        without the extension (if files are DATA.bam,
                        DATA.fam, and DATA.bim, use "--plink DATA" ). Either
                        non-genetic data or PLINK-formatted genetic data is
                        required (see -n/--non-genetic)
  -c COORDS, --coords COORDS
                        TSV (tab-separated values) or similar file with the
                        locations of each sample. One sample per row and
                        longitude and latitude for each sample. May have non-
                        data lines beginning with a #

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output prefix (default: "womble_out")
  --coord-sep COORD_SEP
                        Separator for the positions-file (default: whitespace)
  --non-gen-sep NON_GEN_SEP
                        Separator for the non-genetic data (default:
                        whitespace)
  --density DENSITY     Density (per unit longitude and latitude) at which to
                        compute wombling (default: "(1, 1)" )
  --grid-bounds GRID_BOUNDS
                        Boundary of the grid in which wombling is computed, in
                        longitude and latitude (minLong, maxLong, minLat,
                        maxLat, e.g. (120, 125, -6, 3). If this is not set,
                        the default boundaries are 1.5x DENSITY away from the
                        most extreme coordinates.
  -i, --identical       If selected, the traits from samples at a single
                        location are not averaged before interpolation.
                        Interpolation then treats a single location with two
                        samples as exerting twice as much influence. Only
                        relevant when there are duplicate locations.
  -s, --suppress-save   By default, the script will save the grid's latitudes
                        and longitudes, along with the calculated rate-of-
                        change, direction-of-change, and wighted direction-of-
                        change. This option suppresses these from being
                        output.
  --percentiles PERCENTILES
                        [WIP] Only plot the grid-squares of the top N-th
                        percentile (default: 5, meaning thatonly the top 1 in
                        20 tiles, by absolute rate of change, will be shown).
  --no-plot             Disable plotting and its output.

