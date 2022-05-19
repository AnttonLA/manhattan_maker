"""
    manhattan_maker.py
    ~~~~~~~~~~~~~~~~~~~

    This module is used to create Manhattan plots.
    ~~~~~~~~~~~~~~~~~~~

    Heavily based on the manhattan plot generator from: https://github.com/pgxcentre/manhattan_generator

"""

import argparse
import sys
import os

import numpy as np
import pandas as pd

pd.options.mode.use_inf_as_na = True


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.

    :param msg: the message to print to the user before exiting.
    :type msg: string

    """

    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        """Creates a string representation of the message."""
        return self.message


def main():
    """The main method of the program."""
    # Getting and checking the options
    args = parse_args(args=sys.argv[1:])
    check_args(args)

    # Reading the input file for two point linkage
    file_data = read_input_file(args.in_file, args.use_groups_flag, args)

    # Creating the plots
    create_manhattan_plot(file_data, args)


def safe_main():
    """The main method of the program, with exception handling."""
    try:
        main()

    except KeyboardInterrupt:
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    except ProgramError as e:
        print(e, file=sys.stderr)
        sys.exit(1)


def read_input_file(in_filename, use_groups, options):
    """Reads the input file and returns a pandas DataFrame.

    This function reads any kind of input file, as long as the file is
    tab-separated and that it contains columns with the following headers:

    ======================  ===============================================
            Header                           Description
    ======================  ===============================================
    ``chr``                 The name of the chromosome
    ``snp``                 The name of the marker
    ``pos``                 The physical position (bps)
    ``p_value``             The confidence value (*p value*)
    ``group``               The group the variant belongs to (OPTIONAL)
    ======================  ===============================================

    Note:
        If there is a problem while reading the input file(s),
        a :py:class:`ProgramError` will be raised, and the program will be
        terminated.

    :param in_filename: the input file name
    :type in_filename: string
    :param use_groups: whether to use the different colors for the different groups
    :type use_groups: boolean
    :param options: the command line options
    :type options: :py:class:`argparse.Namespace`
    :returns: a pandas dataframe. It will contain the following names:  'chr', 'pos', 'snp', 'conf',
              and optionally 'group'
    :rtype: pandas.DataFrame

    """
    # Reading the data
    csv_iterator = pd.read_csv(in_filename, sep="\t", chunksize=1e6, low_memory=False,
                               na_values=["NA", "na", "NaN", "None"])

    if not use_groups:
        data = next(csv_iterator).dropna()
        for chunk in csv_iterator:
            data = data.append(chunk.dropna(), ignore_index=True)
    else:  # If we are using groups, we can't drop NA values
        data = next(csv_iterator)
        for chunk in csv_iterator:
            data = data.append(chunk, ignore_index=True)

    # Checking we have the required columns
    required_cols = {options.col_chr,
                     options.col_name,
                     options.col_pos,
                     options.col_pvalue}
    if use_groups:
        required_cols.add(options.col_groups)

    same_col = required_cols & set(data.columns)
    if same_col != required_cols:
        raise ProgramError("{}: missing column(s) {}".format(
            in_filename,
            ", ".join(required_cols - same_col),
        ))

    # Renaming the columns
    data = data.rename(columns={
        options.col_chr: "chrom",
        options.col_pos: "pos",
        options.col_name: "snp",
        options.col_pvalue: "conf",
    })
    if use_groups:
        data = data.rename(columns={options.col_groups: "group"})
        data = data[["chrom", "pos", "snp", "conf", "group"]]

        color_list = options.group_colors.split(',')

        # If data.group has NaN values
        if data.group.isnull().any():
            data.group.fillna("<no data>", inplace=True)     # Filling the empty group column entries
            color_list.append("no_data_color")                       # Adding black color to the list

        if len(color_list) < len(data.group.unique()):
            msg = "The number of group colors is smaller than the number of groups in the file." + \
                  "Please input more colors.\n"
            raise ProgramError(msg)
    else:
        data = data[["chrom", "pos", "snp", "conf"]]

    # Encoding the chromosomes
    data["chrom"] = [encode_chr(chrom) for chrom in data.chrom]

    # Keeping only the required chromosome
    if options.only_chr is not None:
        data = data.loc[data.chrom == options.only_chr, :]
    else:
        if options.exclude_chr:
            data = data[~data.chrom.isin(options.exclude_chr)]

    # Modify p-value to get negative log10
    data["conf"] = -1 * np.log10(data.conf)

    # Ordering and returning
    return data.sort_values(by=["chrom", "pos"])


def encode_chr(chromosome):
    """Encode a chromosome in integer format.

    This function encodes sex chromosomes, pseudo-autosomal regions and
    mitochondrial chromosomes in 23, 24, 25 and 26, respectively. If the
    chromosome is none of the above, the function returns the integer
    representation of the chromosome, if possible.

    Note:
        If the chromosome is invalid, a :py:class:`ProgramError` will be raised, and the program terminated.

    Warning!
    -------

        No check is done whether the chromosome is higher than 26 and below 1.
        As long as the chromosome is an integer or equal to ``X``, ``Y``,
        ``XY`` or ``MT``, no :py:class:`ProgramError` is raised.

    :param chromosome: the chromosome to encode in integer.
    :type chromosome: string
    :returns: the chromosome encoded in integer instead of string.
    :rtype: int
    """
    try:
        return int(chromosome)

    except ValueError:
        chromosome = chromosome.upper()
        if chromosome == 'X':
            return 23
        elif chromosome == 'Y':
            return 24
        elif chromosome == 'XY':
            return 25
        elif chromosome == 'MT':
            return 26

        msg = "%(chromosome)s: not a valid chromosome" % locals()
        raise ProgramError(msg)


def parse_args(args):
    """Parses the command line options and arguments. The input options are the following:

    ============================  =======  ====================================
         Options                   Type                   Description
    ============================  =======  ====================================
    ``--in-file``                 File     The input *file* for two-point
                                           linkage
    ``--output``                  String   The name of the ouput *file*
    ``--format``                  String   The format of the plot (ps, pdf
                                           png)
    ``--dpi``                     Int      The quality of the output (in dpi)
    ``--use-groups``              Boolean  Whether or not to separate the
                                           significant variants into differently
                                           colored groups
    ``--exclude-chr``             String   List of chromosomes to exclude
    ``--only-chr``                String   Only plot the specified chromosome
    ``--no-negative-values``      Boolean  Do not plot negative values
    ``--max-ylim``                Float    The maximal Y *value* to plot
    ``--min-ylim``                Float    The minimal Y *value* to plot
    ``--graph-title``             String   The *title* of the graph
    ``--graph-xlabel``            String   The *text* for the x label
    ``--graph-ylabel``            String   The *text* for the y label
    ``--graph-width``             Float      The *width* of the graph, in
                                           inches
    ``--graph-height``            Float      The *height* of the graph, in
                                           inches
    ``--point-size``              Float    The *size* of each points.
    ``--significant-point-size``  Float    The *size* of each significant
                                           points
    ``--abline``                  String   The y *value* where to create a
                                           horizontal line, separated by a comma
    ``--significant-threshold``   Float    The significant threshold for
                                           linkage
    ``--chromosome-box-color``    String   The *color* for the box surrounding
                                           even chromosome numbers
    ``--even-chromosome-color``   String   The *color* for the box surrounding
                                           even chromosome numbers
    ``--odd-chromosome-color``    String   The *color* for the box surrounding
                                           odd chromosome numbers
    ``--significant-color``       String   The *color* for points representing
                                           significant linkage
    ``--group-colors``            String   The *color* for each of the groups
    ============================  =======  ====================================

    Note:
        No option check is done here (except for the one automatically done by :py:mod:`argparse`. Those need to be
        done elsewhere (see :py:func:`check_args`).

    :returns: An object created by the :py:mod:`argparse` module. It contains the values of the different options.
    :rtype: argparse.Namespace
    """

    # Creating the parser object
    parser = argparse.ArgumentParser(
        description="This script produces nice Manhattan plots for either "
                    "linkage or GWAS results.",
    )

    parser.add_argument("--mode", default='client')
    parser.add_argument("--port", default=38475)

    # --- The input options ---
    input_options_group = parser.add_argument_group(
        "Input Options",
        "Options for the input file(s) (name of the file, type of graph, etc.)",
    )
    # The input file (for two point)
    input_options_group.add_argument(
        "--in-file", type=str, metavar="FILE",
        help="The input FILE.",
    )

    # --- The column options ---
    column_options_group = parser.add_argument_group(
        "Column Options",
        "The name of the different columns in the input file(s).",
    )
    # The chromosome column
    column_options_group.add_argument(
        "--col-chr", type=str, metavar="COL", default="chr",
        help="The name of the column containing the chromosomes [Default: %(default)s].",
    )
    # The marker name column
    column_options_group.add_argument(
        "--col-name", type=str, metavar="COL", default="name",
        help="The name of the column containing the marker names [Default: %(default)s].",
    )
    # The marker position column
    column_options_group.add_argument(
        "--col-pos", type=str, metavar="COL", default="pos",
        help="The name of the column containing the marker positions [Default: %(default)s].",
    )
    # The marker p value
    column_options_group.add_argument(
        "--col-pvalue", type=str, metavar="COL", default="p_value",
        help="The name of the column containing the marker p values [Default: %(default)s]",
    )
    # The groups column
    column_options_group.add_argument(
        "--col-groups", type=str, metavar="COL", default="groups",
        help="The name of the column containing the group names [Default: %(default)s].",
    )

    # --- The output options ---
    output_options_group = parser.add_argument_group(
        "Graph Output Options",
        "Options for the output file (name of the file, type of graph, etc.).",
    )
    # The output file name
    output_options_group.add_argument(
        "-o", "--output", dest="outFile_name", type=str, default="manhattan",
        metavar="NAME",
        help="The NAME of the output file [Default: %(default)s].",
    )
    # The type of the graph (png, ps or pdf)
    format_choices = ["ps", "pdf", "png", "eps"]
    output_options_group.add_argument(
        "-f", "--format", dest="graph_format", type=str, default="png",
        metavar="FORMAT", choices=format_choices,
        help="The FORMAT of the plot ({}) "
             "[Default: %(default)s].".format(", ".join(format_choices)),
    )
    output_options_group.add_argument(
        "--web", action="store_true",
        help="Always write a PNG file for web display, and return the path of the PNG file.",
    )
    output_options_group.add_argument(
        "--dpi", type=int, default=600, metavar="INT",
        help="The quality of the output (in dpi) [Default: %(default)d].",
    )

    # --- The graph type options ---
    graph_options_group = parser.add_argument_group(
        "Graph Options",
        "Options for the graph type (two-point, multipoint, etc.).",
    )
    # Using p values instead of LOD score
    graph_options_group.add_argument(
        "--use-groups", dest='use_groups_flag', action="store_true",
        help="Use group information to color the variants differently.",
    )
    # Exclude some chromosomes
    graph_options_group.add_argument(
        "--exclude-chr", metavar="STRING",
        help="Exclude those chromosomes (list of chromosomes, separated by a "
             "coma) [Default: None].",
    )
    graph_options_group.add_argument(
        "--only-chr", metavar="CHR",
        help="Print only the results for a single chromosome. The xaxis will "
             "hence show the positions/cm instead of the chromosome number.",
    )

    # --- The graph presentation options ---
    graph_aesthetics_options_group = parser.add_argument_group(
        "Graph Presentation Options",
        "Options for the graph presentation (title, axis label, etc.).",
    )
    # print negative values
    graph_aesthetics_options_group.add_argument(
        "--no-negative-values", action="store_true",
        help="Do not plot negative values.",
    )
    # The maximal y limit of the graph
    graph_aesthetics_options_group.add_argument(
        "--max-ylim", type=float, metavar="FLOAT",
        help="The maximal Y value to plot [Default: maximum of max(LOD) "
             "and 1+significant-threshold].",
    )
    # The minimal y limit of the graph
    graph_aesthetics_options_group.add_argument(
        "--min-ylim", type=float, default=-2.0, metavar="FLOAT",
        help="The minimal Y value to plot [Default: %(default).1f].",
    )
    # Do we want padding?
    graph_aesthetics_options_group.add_argument(
        "--no-y-padding", action="store_true",
        help="Do not add Y padding to the Y limit",
    )
    # The graph's title
    graph_aesthetics_options_group.add_argument(
        "--graph-title", type=str, dest='graph_title', default="",
        metavar="TITLE",
        help="The TITLE of the graph [Default: empty].",
    )
    # The graph's x label
    graph_aesthetics_options_group.add_argument(
        "--graph-xlabel", dest='graph_x_label', type=str, default="Chromosome",
        metavar="TEXT",
        help="The TEXT for the x label. [Default: %(default)s].",
    )
    # The graph's y label
    graph_aesthetics_options_group.add_argument(
        "--graph-ylabel", dest='graph_y_label', type=str, default="LOD",
        metavar="TEXT",
        help="The TEXT for the y label. [Default: %(default)s].",
    )
    # The graph width
    graph_aesthetics_options_group.add_argument(
        "--graph-width", type=float, default=14.0, metavar="WIDTH",
        help="The WIDTH of the graph, in inches [Default: %(default)d].",
    )
    # The graph height
    graph_aesthetics_options_group.add_argument(
        "--graph-height", type=float, default=7.0, metavar="HEIGHT",
        help="The HEIGHT of the graph, in inches [Default: %(default)d].",
    )
    # The size of each point
    graph_aesthetics_options_group.add_argument(
        "--point-size", type=float, default=2.1, metavar="SIZE",
        help="The SIZE of each points [Default: %(default).1f].",
    )
    # The size of each significant point
    graph_aesthetics_options_group.add_argument(
        "--significant-point-size", type=float, default=4.5, metavar="SIZE",
        help="The SIZE of each significant points [Default: %(default).1f].",
    )
    # The ablines positions
    graph_aesthetics_options_group.add_argument(
        "--abline", type=str, default="3,-2", metavar="POS1,POS2,...",
        help="The y value where to create a horizontal line, separated by a "
             "comma [Default: %(default)s].",
    )
    # The significant threshold
    graph_aesthetics_options_group.add_argument(
        "--significant-threshold", type=float, default=3.0, metavar="FLOAT",
        help="The significant threshold for linkage or association "
             "[Default: %(default).1f]",
    )
    # The size of the text
    graph_aesthetics_options_group.add_argument(
        "--axis-text-size", type=int, default=12, metavar="INT",
        help="The axis font size [Default: %(default)d]",
    )
    # The size of the text for the chromosome names
    graph_aesthetics_options_group.add_argument(
        "--chr-text-size", type=int, default=12, metavar="INT",
        help="The chromosome font size [Default: %(default)d]",
    )
    # The size of the text for the labels
    graph_aesthetics_options_group.add_argument(
        "--label-text-size", type=int, default=12, metavar="INT",
        help="The label font size [Default: %(default)d]",
    )

    # --- The graph color options ---
    color_options_group = parser.add_argument_group(
        "Graph Colors Options",
        "Options for the graph colors.",
    )
    color_options_group.add_argument(
        "--chromosome-box-color", type=str, default="#E5E5E5", metavar="COLOR",
        help="The COLOR for the box surrounding even chromosome numbers [Default: %(default)s].",
    )
    color_options_group.add_argument(
        "--even-chromosome-color", type=str, default="#1874CD",
        metavar="COLOR",
        help="The COLOR for the box surrounding even chromosome numbers [Default: %(default)s].",
    )
    color_options_group.add_argument(
        "--odd-chromosome-color", type=str, default="#4D4D4D", metavar="COLOR",
        help="The COLOR for the box surrounding odd chromosome numbers [Default: %(default)s].",
    )
    color_options_group.add_argument(
        "--significant-color", type=str, default="#FF0000", metavar="COLOR",
        help="The COLOR for points representing significant linkage [Default: %(default)s].",
    )
    color_options_group.add_argument(
        "--group-colors", type=str, default=None, metavar="COLOR",
        help="The COLOR for the groups, separated by a comma [Default: %(default)s].",
    )
    color_options_group.add_argument(
        "--plot-margins", type=float, default=None,
        help="The margins for the plot passed to plt.pyplot.margins [Default: %(default)s].",
    )

    return parser.parse_args(args)


def check_args(args):
    """Checks the arguments and options.

    Note:
        If there is a problem with an option, an exception is raised using the
        :py:class:`ProgramError` class, a message is printed to the
        :class:`sys.stderr` and the program exists with code 1.

    :param args: a :py:class:`Namespace` object containing the options of the program.
    :type args: argparse.Namespace
    :return: None.
    """

    if (args.max_ylim is not None) and (args.min_ylim is not None):
        if args.max_ylim <= args.min_ylim:
            msg = "Y max limit (%f) is <= Y min limit (%f)" % (args.max_ylim, args.min_ylim)
            raise ProgramError(msg)

    # The type of graph
    if args.in_file is None:
        msg = "No input file was specified."
        raise ProgramError(msg)

    # Check for input file (two-point)
    if args.in_file is not None:
        if not os.path.isfile(args.in_file):
            msg = "%s: no such file or directory" % args.in_file
            raise ProgramError(msg)

    try:
        args.abline = [float(i) for i in args.abline.split(',')]
    except ValueError:
        msg = "%s: not a valid LOD score (must be float)" % args.abline
        raise ProgramError(msg)

    # Checking that only one of --exclude-chr or --only-chr is used
    if args.exclude_chr is not None and args.only_chr is not None:
        msg = "Use only one of '--exclude-chr' or '--only-chr'"
        raise ProgramError(msg)

    # Checking if there are some chromosome to exclude
    if args.exclude_chr is None:
        args.exclude_chr = set()
    else:
        args.exclude_chr = {encode_chr(i) for i in args.exclude_chr.split(",")}

    if args.only_chr is not None:
        args.only_chr = encode_chr(args.only_chr)

    # Check if '--use_groups' was used on its own. If so, give default values to --group-colors
    if args.use_groups_flag & (args.group_colors is None):  # TODO: this is awkwardly located here. It may need to move.
        args.group_colors = "#FF0000,#0000FF"
    elif not args.use_groups_flag & (args.group_colors is not None):
        msg = "--group-colors option requires --use-groups to be used"
        raise ProgramError(msg)


def create_manhattan_plot(plotting_data, args):
    """Creates the manhattan plot from marker data.

    Creates manhattan plots from the given data. Results are shown in a manhattan plot using points (different color for
    each of the chromosomes).

    :param plotting_data: a pandas dataframe containing the data to plot.
    :type plotting_data: pandas.DataFrame
    :param args: a :py:class:`Namespace` object containing the options of the program.
    :type args: argparse.Namespace
    :return: None.
    """

    import matplotlib as mpl

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ProgramError("Could not import matplotlib.")

    plt.ioff()

    # The available chromosomes
    available_chrom = [sorted(plotting_data.chrom.unique())]
    if len(available_chrom) == 1:
        available_chrom = available_chrom[0]

    # Creating the figure
    figure = plt.figure(figsize=(args.graph_width, args.graph_height), frameon=True)

    # Getting the maximum and minimum of the confidence value
    conf_min = [0.0]
    conf_max = []

    tmp_conf_vals = plotting_data.conf
    tmp_conf_vals.dropna(inplace=True)
    tmp_conf_vals.dropna(inplace=True)
    conf_min.append(tmp_conf_vals.min())  # TODO: Add a flag that allows for the plot to NOT start at zero
    conf_max.append(tmp_conf_vals.max())

    conf_min = min(conf_min)
    conf_max = max(conf_max)
    if args.max_ylim is not None:
        conf_max = args.max_ylim
    if args.min_ylim is not None:
        conf_min = args.min_ylim
    if args.no_negative_values:
        conf_min = 0.0

    # The chromosome spacing
    chrom_spacing = 25000000

    # Creating the ax and modify it
    ax = figure.add_subplot(111)
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("left")
    ax.set_ylabel("LOD")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)

    ax.spines["left"].set_linewidth(1.5)  # Tweak to axis line width, as requested by Björn
    ax.spines["bottom"].set_linewidth(1.5)
    ax.yaxis.set_tick_params(width=1.5, length=5)  # Tweak to tick width and length

    if args.only_chr is None:
        ax.set_xticks([])
        ax.set_xticklabels([])
    ax.set_ylabel(r'$-\log_{10}$ $\it{P}$', fontsize=args.label_text_size, fontname='arial')

    ax.set_xlabel(args.graph_x_label, fontsize=args.label_text_size, fontname='arial')
    # ANTTON - Fontsize used to be 16 in the original
    ax.set_title(args.graph_title, fontsize=args.label_text_size, weight="bold", fontname='arial')

    if args.use_groups_flag:  # Create dict used to assign colors to groups
        group_colors_dict = {'<no data>': args.significant_color}
        color_list = args.group_colors.split(',')  # Splitting the group colors into a list

    # Now plotting for each of the chromosomes
    starting_pos = 0
    ticks = []
    for i, chrom in enumerate(available_chrom):
        chrom_data = None
        max_pos = []

        chrom_data = plotting_data[plotting_data.chrom == chrom]
        max_pos.append(chrom_data.pos.max())

        max_pos = max(max_pos)

        # The color of the points
        color = args.even_chromosome_color
        if i % 2 == 0:
            color = args.odd_chromosome_color

        # The box
        xmin = starting_pos - (chrom_spacing / 2)
        xmax = max_pos + starting_pos + (chrom_spacing / 2)
        if i % 2 == 1:
            ax.axvspan(xmin=xmin, xmax=xmax, color=args.chromosome_box_color)

        # The chromosome label
        ticks.append((xmin + xmax) / 2)

        # Plotting
        if args.in_file is not None:
            ax.plot(chrom_data.pos + starting_pos, chrom_data.conf,
                    marker="o", ms=args.point_size, mfc=color, mec=color,
                    ls="None")

        # Plotting the abline
        for abline_position in args.abline:
            ax.axhline(y=abline_position, color="black", ls="--", lw=0.6)
        if conf_min < 0:
            ax.axhline(y=0, color="black", ls="-", lw=1.0)

        # Plotting the significant markers
        if not args.use_groups_flag:
            sig_mask = chrom_data.conf >= args.significant_threshold
            ax.plot(chrom_data.pos[sig_mask] + starting_pos,
                    chrom_data.conf[sig_mask], marker="o", ls="None",
                    ms=args.significant_point_size, mfc=args.significant_color,
                    mec=args.significant_color)
        else:
            for group in chrom_data.group.unique():
                if group not in group_colors_dict.keys():
                    group_colors_dict[group] = color_list.pop(0)
                sig_mask = chrom_data.conf >= args.significant_threshold
                sig_mask = sig_mask & (chrom_data.group == group)
                ax.plot(chrom_data.pos[sig_mask] + starting_pos,
                        chrom_data.conf[sig_mask], marker="o", ls="None",
                        ms=args.significant_point_size, mfc=group_colors_dict[group],
                        mec=group_colors_dict[group], label=group)

        # Changing the starting point for the next chromosome
        starting_pos = max_pos + starting_pos + chrom_spacing

    # Setting the limits
    padding = 0.39
    if args.no_y_padding:
        padding = 0
    ax.set_ylim(conf_min - padding, conf_max + padding)
    ax.set_xlim(0 - chrom_spacing, starting_pos + chrom_spacing)

    # Putting the xticklabels
    if args.only_chr is None:
        ax.set_xticks(ticks)
        ax.set_xticklabels(available_chrom)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(args.axis_text_size)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(args.chr_text_size)

    # Saving or plotting the figure
    mpl.rcParams['axes.linewidth'] = 0.1  # ADDED BY LUDVIG
    mpl.rcParams['savefig.dpi'] = args.dpi
    mpl.rcParams['ps.papersize'] = "auto"
    mpl.rcParams['savefig.orientation'] = "landscape"

    plt.margins(x=None, y=0.0, tight=True)  # ADDED BY LUDVIG
    if args.use_groups_flag:
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))  # Get rid of duplicate labels
        if '<no data>' in by_label.keys():  # Remove the '<no data>' label
            by_label.pop('<no data>')
        plt.legend(by_label.values(), by_label.keys(),
                   loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, ncol=1)

    plt.savefig(args.outFile_name + "." + args.graph_format, bbox_inches="tight")

    if args.graph_format != "png":  # If a format other than png is used, a png file is also created
        plt.savefig(args.outFile_name + ".png", bbox_inches="tight")

    if args.web:
        print(args.outFile_name + ".png")
    else:
        # There is some two-point data and annotation is asked, se we show
        # the figure
        plt.show()


if __name__ == "__main__":
    safe_main()
