# Manhattan Maker

`manhattan_maker` is a python tool to create Manhattan plots from GWAS output data.

It is heavily based on `manhattan_generator` (https://github.com/pgxcentre/manhattan_generator), the only differences
being that some options have been removed, and a functionality has been added to color the significant variants by the
group they belong to.


## Dependencies

Here are the dependencies of the tool:

- [Python](http://python.org/) version 2.7 or 3.4 or latest
- [numpy](http://www.numpy.org/) version 1.8.0 or latest
- [matplotlib](http://matplotlib.org/) version 1.3.1 or latest
- [pandas](http://pandas.pydata.org/) version 0.17 or latest


## Installation

We recommend installing the tool in a Python virtual environment.

`manhattan_generator` should work on Windows and MacOS, even though it hasn't
been fully tested for full compatibility.


## Basic usage

```console
$ manhattan_maker --help
usage: manhattan_maker [-h] [--in-file FILE] [--col-chr COL] [--col-name COL] [--col-pos COL] [--col-pvalue COL]
                       [--col-groups COL] [-o NAME] [-f FORMAT] [--web] [--dpi INT] [--use-groups]
                       [--exclude-chr STRING] [--only-chr CHR] [--no-negative-values] [--max-ylim FLOAT]
                       [--min-ylim FLOAT] [--no-y-padding] [--graph-title TITLE] [--graph-xlabel TEXT]
                       [--graph-ylabel TEXT] [--graph-width WIDTH] [--graph-height HEIGHT] [--point-size SIZE]
                       [--significant-point-size SIZE] [--abline POS1,POS2,...] [--significant-threshold FLOAT]
                       [--axis-text-size INT] [--chr-text-size INT] [--label-text-size INT]
                       [--chromosome-box-color COLOR] [--even-chromosome-color COLOR] [--odd-chromosome-color COLOR]
                       [--significant-color COLOR] [--group-colors COLOR] [--plot-margins PLOT_MARGINS]


This script produces nice Manhattan plots for GWAS results.

optional arguments:
  -h, --help            show this help message and exit

Input Options:
  Options for the input file(s) (name of the file, type of graph, etc.)

  --in-file FILE        The input FILE.

Column Options:
  The name of the different columns in the input file(s).

  --col-chr COL         The name of the column containing the chromosomes [Default: chr].
  --col-name COL        The name of the column containing the marker names [Default: name].
  --col-pos COL         The name of the column containing the marker positions [Default: pos].
  --col-pvalue COL      The name of the column containing the marker p values [Default: p_value]
  --col-groups COL      The name of the column containing the group names [Default: groups].

Graph Output Options:
  Options for the output file (name of the file, type of graph, etc.).

  -o NAME, --output NAME
                        The NAME of the output file [Default: manhattan].
  -f FORMAT, --format FORMAT
                        The FORMAT of the plot (ps, pdf, png, eps) [Default: png].
  --web                 Always write a PNG file for web display, and return the path of the PNG file.
  --dpi INT             The quality of the output (in dpi) [Default: 600].

Graph Options:
  Options for the graph type (two-point, multipoint, etc.).

  --use-groups          Use group information to color the variants differently.
  --exclude-chr STRING  Exclude those chromosomes (list of chromosomes, separated by a coma) [Default: None].
  --only-chr CHR        Print only the results for a single chromosome. The xaxis will hence show the positions/cm instead of the chromosome number.

Graph Presentation Options:
  Options for the graph presentation (title, axis label, etc.).

  --no-negative-values  Do not plot negative values.
  --max-ylim FLOAT      The maximal Y value to plot [Default: maximum of max(LOD) and 1+significant-threshold].
  --min-ylim FLOAT      The minimal Y value to plot [Default: -2.0].
  --no-y-padding        Do not add Y padding to the Y limit
  --graph-title TITLE   The TITLE of the graph [Default: empty].
  --graph-xlabel TEXT   The TEXT for the x label. [Default: Chromosome].
  --graph-ylabel TEXT   The TEXT for the y label. [Default: LOD].
  --graph-width WIDTH   The WIDTH of the graph, in inches [Default: 14].
  --graph-height HEIGHT
                        The HEIGHT of the graph, in inches [Default: 7].
  --point-size SIZE     The SIZE of each points [Default: 2.1].
  --significant-point-size SIZE
                        The SIZE of each significant points [Default: 4.5].
  --abline POS1,POS2,...
                        The y value where to create a horizontal line, separated by a comma [Default: 3,-2].
  --significant-threshold FLOAT
                        The significant threshold for linkage or association [Default: 3.0]
  --axis-text-size INT  The axis font size [Default: 12]
  --chr-text-size INT   The chromosome font size [Default: 12]
  --label-text-size INT
                        The label font size [Default: 12]

Graph Colors Options:
  Options for the graph colors.

  --chromosome-box-color COLOR
                        The COLOR for the box surrounding even chromosome numbers [Default: #E5E5E5].
  --even-chromosome-color COLOR
                        The COLOR for the box surrounding even chromosome numbers [Default: #1874CD].
  --odd-chromosome-color COLOR
                        The COLOR for the box surrounding odd chromosome numbers [Default: #4D4D4D].
  --significant-color COLOR
                        The COLOR for points representing significant linkage [Default: #FF0000].
  --group-colors COLOR  The COLOR for the groups, separated by a comma [Default: None].
  --plot-margins PLOT_MARGINS
                        The margins for the plot passed to plt.pyplot.margins [Default: None].

```

