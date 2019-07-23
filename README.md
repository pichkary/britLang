# britLang
A tool to detect genetic and non-genetic barriers in geography, written in Python 3.7

Required python packages:
  * numpy 
  * matplotlib
  * cHaversine
  * basemap  (from matplotlib_toolkits)
  * plinkio  (when using genetic data)
  * scipy (for scipy.ndimage.label)

### Usage
    python3 womble.py [-h] (-n NON_GENETIC | -p PLINK) -c COORDS
                      [-o OUTPUT] [--coord-sep COORD_SEP]
                      [--non-gen-sep NON_GEN_SEP] [-d DENSITY]
                      [--grid-bounds GRID_BOUNDS] [-i] [-s]
                      [--percentile PERCENTILE] [--no-plot]
                      [--no-plot-save] [--alpha ALPHA]

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
      -d DENSITY, --density DENSITY
                            Density (per unit longitude and latitude) at which to
                            compute wombling (default: (1, 1) )
      --grid-bounds GRID_BOUNDS
                            Boundary of the grid in which wombling is computed, in
                            longitude and latitude (minLong, maxLong, minLat,
                            maxLat, e.g. "(120, 125, -6, 3)"). If this is not set,
                            the default boundaries are set using the most extreme
                            coordinates in COORDS.
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
      --percentile PERCENTILE
                            Only plot the grid-squares of the top N-th percentile
                            (default: 5. The top 5th percentile will select only
                            the top 1 in 20 tiles by absolute rate of change).
      --no-plot             Disable plotting and its output.
      --no-plot-save        If plotting, only display (do not save) the plot.
      --alpha ALPHA         Alpha value (opacity) for rate-of-change heat-map
                            (Default: 0.5)

To plot saved data (Note: either [prefix from womble.py & density] or [files containing rates, directions, longitudes, & latitudes] are required):

    usage: Womble plotter [-h] (([-p PREFIX] [--density DENSITY]) | ([-r RATES]
                          [-d DIRECTIONS] [-g LONGITUDES] [-t LATITUDES]))
                          [-c PERCENTILE] [-a ALPHA]
    
    required arguments: 
      -p PREFIX, --prefix PREFIX
                            Prefix for rates-of-change, directions, longitudes,
                            and latitudes files produced by womble.py. Uses
                            weighted direction-of-change by default. Individual
                            files can be overridden by calling other arguments
                            (e.g. -d direction-of-change.txt).
      --density DENSITY     Density (per unit longitude and latitude) used to
                            generate grid. Required if using the "prefix" option
                            (tuple of length 2, e.g.: "(2, 2)" ).
      -r RATES, --rates RATES
                            Rates-of-change file generated by womble.py for each
                            point on the grid.
      -d DIRECTIONS, --directions DIRECTIONS
                            Directions-of-change file generated by womble.py for
                            every point of the grid. (weighted direction of change
                            is preferred.
      -g LONGITUDES, --longitudes LONGITUDES
                            File with longitudes of the grid generated by
                            womble.py
      -t LATITUDES, --latitudes LATITUDES
                            File with latitudes of the grid generated by womble.py

    optional arguments:
      -h, --help            show this help message and exit
      -c PERCENTILE, --percentile PERCENTILE
                            Only plot the grid-squares of the top N-th percentile
                            (default: 5, meaning that the top 1 in 20 tiles (top
                            by absolute rate of change) will be shown).
      -a ALPHA, --alpha ALPHA
                            The opacity of the heat-map of rates of change
                            (Default: 0; 0 <= alpha <= 1).
    
