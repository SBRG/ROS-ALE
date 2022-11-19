"""
Functions for ROS TALE Manuscript
"""
__author__ = 'Kevin Rychel'
__email__ = 'krychel@eng.ucsd.edu'


############
## Set Up ##
############

# import general packages
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# iModulon functions; see PyModulon documentation
# github.com/SBRG/pymodulon
from pymodulon.io import *
from pymodulon.plotting import *

# import some more specific functions
from matplotlib.colors import to_rgba

# modify matplotlib settings for output
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
sns.set_style('ticks')
matplotlib.rcParams['text.color'] = '#000000'
matplotlib.rcParams['axes.labelcolor'] = '#000000'
matplotlib.rcParams['xtick.color'] = '#000000'
matplotlib.rcParams['ytick.color'] = '#000000'

######################
## Helper Functions ##
######################

def get_ev_strains(ros_meta):
    """
    Get a list of all the evolved strains
    
    Parameter
    ---------
    ros_meta: metadata for the samples
    
    Returns
    -------
    strains: list of strains
    """
    
    strains = ros_meta.lineage.unique()
    strains = [i for i in strains if not('0_0' in i)]
    strains = [i for i in strains if not('mid' in i)]
    
    return strains

# input: ros_meta
#        lineages, a list of lineages, e.g. ['0_0', '1_0']
#        pq, a list of paraquat concentrations, e.g. [0, 250]
# output: a list of samples in the desired conditions
def get_samples(ros_meta, lineages, pq):
    """
    Gets sample names for the given condition
    
    Parameters
    ----------
    ros_meta: metadata for the samples
    lineages: a list of strains, e.g. ['0_0', '1_0']
    pq: a list of paraquat concentrations, e.g. [0, 250]
    
    Returns
    -------
    samples: a Pandas series of sample names
    
    """
    
    lin_bool = ros_meta.lineage.isin(lineages)
    pq_bool = ros_meta.pq.isin(pq)
    series = ros_meta.index
    return series[lin_bool & pq_bool]

def get_pq_dicts():
    """
    Returns dictionaries for PQ plot settings
    
    Returns
    -------
    pq_color: colors for each PQ value
    pq_marker: marker shapes for each PQ value
    """
    pq_marker = {0:'o',
                 250: 's',
                 750: '^'}
    pq_color = {0: 'tab:blue',
                125: '#6D6CBA',
                250: 'tab:purple', 
                750:'tab:red'}
    pq_color = {k:to_rgba(v) for k, v in pq_color.items()}
    
    return pq_color, pq_marker

def get_mutations(all_mutations, gene):
    """
    Gets the set of mutations affecting a gene

    Parameters
    ----------
    all_mutations: table of mutations
    gene: name of gene, e.g. 'emrE'

    Returns
    -------
    mutations: Pandas DataFrame with details for relevant mutations
    """
    muts = []
    for i, row in all_mutations.iterrows():
        mut_genes = row.Gene.split(', ')
        if gene in mut_genes:
            muts += [i]
    return all_mutations.loc[muts]

####################
## iModulon Plots ##
####################

def im_bar_swarm(ica, ros_meta, k, strain_dict = None, ax = None,
              orientation = 'vertical', dodge = True, legend = False,
              legendsize = 7, labelsize = 7, ticksize = 7, 
              errwidth = 1, pointsize = 1.5, x_rot = 90):
    """
    Generates bar/swarm plots for iModulon activities
    Relies on seaborn's bar and swarm functions
    
    Parameters
    ----------
    ica, ros_meta: main object & metadata
    k: iModulon name
    strain_dict: dictionary with labels pointing to lists of strains
                 default: start vs. evolved
    ax: axes to plot to (will create one if not provided)
    orientation: 'horizontal' or 'vertical'
    dodge: Whether to separate each PQ level into different bars
    legend: Whether to draw the legend
    legendsize: fontsize for the legend
    labelsize: fontsize for the axis labels
    ticksize: fontsize for the tick labels on the axes
    err_width: width of errorbars in barplot
    pointsize: point size for swarm plot
    x_rot: rotation of x axis labels
    
    Returns
    -------
    ax: axes with plot
    """
    
    # organize groups of strains
    if strain_dict is None:
        strain_dict = {'Start':['0_0'],
                       'Evolved':get_ev_strains(ros_meta)}
    sample_dict = {label: get_samples(ros_meta, strains, [0, 250, 750]) 
                   for label, strains in strain_dict.items()}

    # make data table for seaborn
    data = pd.DataFrame(columns = ['x_label', 'y', 'Lineage', 'PQ'])
    for label, samples in sample_dict.items():
        for s in samples:
            data.loc[s] = [label, ica.A.loc[k, s], 
                           ros_meta.lineage[s],
                           ros_meta.pq[s]]
    
    # create axes if necessary
    if ax is None:
        fig, ax = plt.subplots(figsize = (0.75, 1.5),
                               dpi = 140)
    
    # define labels & columns based on orientation
    if not orientation in ['horizontal', 'vertical']:
        print('assuming vertical',
              '(else enter orientation = "horizontal")')
        orientation = 'vertical'
    if orientation == 'horizontal':
        x, y = 'y', 'x_label'
        x_label, y_label = k + ' iModulon Activity', ''
        ticksize_adjust = 'x'
    elif orientation == 'vertical':
        x, y = 'x_label', 'y'
        y_label, x_label = k + ' iModulon Activity', ''
        ticksize_adjust = 'y'
        
    # get PQ settings
    pq_color, _ = get_pq_dicts()
    
    # draw plots based on dodge setting
    if dodge:
        ax = sns.barplot(x = x, y = y, data = data, 
                         hue = 'PQ', palette = pq_color,
                         ax = ax, alpha = 0.7)
        ax = sns.swarmplot(x = x, y = y, data = data, 
                           hue = 'PQ', palette = pq_color,
                           s = pointsize, ax = ax, dodge = True)
    else:
        ax = sns.barplot(x = x, y = y, data = data, 
                         color = 'tab:gray',
                         ax = ax, alpha = 0.5)
        ax = sns.swarmplot(x = x, y = y, data = data, 
                           hue = 'PQ', palette = pq_color,
                           s = pointsize, ax = ax, dodge = False)
        
    # deal with the legend
    l = ax.legend(bbox_to_anchor = (1,1),loc = 'upper left', 
                  fontsize = legendsize, title = 'PQ (μM)')
    plt.setp(l.get_title(), fontsize=legendsize)
    if legend == False:
        l.remove()
    
    # ticks and labels
    ax.tick_params(labelsize = labelsize)
    ax.tick_params(ticksize_adjust, labelsize = ticksize)
    ax.set_ylabel(y_label, fontsize = labelsize)
    ax.set_xlabel(x_label, fontsize = labelsize)
    if x_rot != 0:
        ax.tick_params('x', rotation = x_rot)
    
    return ax

####################
## Mutation Plots ##
####################

def mut_color_table(all_mutations, ros_meta, genes, label_df,
                    ax, fontsize = 7, skip_muts = []):
    """
    Generates colored mutation heatmap for a set of genes.
    The labels need to be added to the plot in Adobe Illustrator,
        so they are output as a separate table here.

    Parameters
    ----------
    all_mutations: table of mutations
    ros_meta: metadata for this study
    genes: list of gene names
    label_df: table matching mutation names to their shorter labels and types
    ax: axes to plot to
    fontsize: fontsize applied to all labels
    skip_muts: mutation names to ignore

    Returns
    -------
    ax: axes with plot
    labeled_table: table with labels
    """

    # get the set of mutations
    muts = []
    for g in genes:
        g_muts = get_mutations(all_mutations, g)
        muts += g_muts.index.to_list()

    # use this to color code the heatmap
    # a default of 7 or any other unused value below 10 is fine
    label_colors = {'SNP': 0, #'tab:blue',
                  'upstream': 2, #'tab:green', 
                  'amp': 4, #'tab:purple', 
                  'frameshift': 1, #'tab:orange',
                  'del': 3,
                  'deletion': 3,
                  'nonsense SNP': 5,
                  'insertion': 8}

    # prepare tables
    strains = get_ev_strains(ros_meta)
    labeled_table = pd.DataFrame('', index = genes, columns = strains)
    t_hmap_val = pd.DataFrame(7, index = genes, columns = strains)

    # fill in tables
    for g in genes:
        g_muts = get_mutations(all_mutations, g)
        for m in g_muts.index:
            if m in skip_muts:
                continue
            for s in strains:
                if g_muts.loc[m, s]:
                    # account for multiple mutations in same gene
                    if len(labeled_table.loc[g, s]) > 0:
                        labeled_table.loc[g, s] = ', '.join([labeled_table.loc[g, s], 
                            label_df.label[m]])
                    else:
                        labeled_table.loc[g, s] = label_df.label[m]
                    # heatmap will only be colored by last mut type observed
                    # be careful of this when labeling in Illustrator
                    t_hmap_val.loc[g, s] = label_colors[label_df.type[m]]

    # make the heatmap
    sns.heatmap(t_hmap_val, cmap = sns.color_palette('tab10'), 
                vmin = 0, vmax = 10, annot = False, ax = ax,
                cbar = False, xticklabels = 1, yticklabels = 1)

    # change fontsizes
    ax.set_yticklabels(ax.get_yticklabels(), fontsize = fontsize)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize = fontsize)

    return ax, labeled_table

###########################
## Genome Coverage Plots ##
###########################

def read_gff(gff_file, left, right, scale):
    """
    Reads GFF file for genome_plot()
    
    Parameters
    ----------
    gff_file: path to gff file to read
    left, right: integer bounds for genome location
    scale: whether to normalize to the max observed coverage
    
    Returns
    -------
    plus, minus: pandas Series of coverage by nucleotide sequence
                 for each strand
    """
    
    # sometimes the first line is a comment which pandas can't handle
    skiprows = 0
    with open(gff_file, "r") as infile:
        if infile.read(1) == "#":
            skiprows = 1
    table = pd.read_table(gff_file, header=None,
        usecols=[0, 2, 3, 4, 5, 6], comment="#", skiprows=skiprows,
        names=["chromosome", "name", "leftpos", "rightpos", 
               "reads", "strand"])
    table = table[(table.rightpos >= left) & (table.leftpos <= right)]
    
    if (table.leftpos != table.rightpos).any():
        print('Need to account for each line not being a single point')
        
    table = table[["leftpos", "reads", "strand"]]
    table_plus = table[table.strand == "+"].set_index("leftpos")
    table_minus = table[table.strand == "-"].set_index("leftpos")
    
    # fill missing values with 0
    filler = pd.Series([range(left, right + 1)], 
                       [range(left, right + 1)])
    table_plus["filler"] = 1
    table_minus["filler"] = 1
    table_plus.fillna(0)
    table_minus.fillna(0)
    
    # extract only the series we need
    plus = table_plus.reads
    minus = table_minus.reads.abs()  # in case stored negative
    if scale:
        plus *= 100. / plus.max()
        minus *= 100. / minus.max()
    
    # downsample to 2000 pts
    collapse_factor = int((right - left) / 2000)
    if collapse_factor > 1:
        plus = plus.groupby(lambda x: x // collapse_factor).mean()
        plus.index *= collapse_factor
        minus = minus.groupby(lambda x: x // collapse_factor).mean()
        minus.index *= collapse_factor

    return plus, minus

def genome_plot(gff_file, left, right, ymax = None,
                ax = None, fontsize = 7, ylabel = 'Coverage',
                scale = False, legend = True):
    
    """
    Generates filled plots for genome coverage from DNAseq
    
    Parameters
    ----------
    gff_file: path to gff file to read
    left, right: integer bounds for genome location
    ymax: if provided, set y axis bounds to [0, ymax]
    ax: axes to plot to (will create if not provided)
    fontsize: fontsize for all text
    ylabel: ylabel to use
    scale: whether to normalize to the max observed coverage
    legend: whether to include a legend
    
    Returns
    -------
    ax: axes with plot
    """
    
    # read gff file
    final_plus, final_minus = read_gff(gff_file, left, right, scale)
    
    # generate axes if necessary
    if ax is None:
        fig, ax = plt.subplots(figsize = (4, 1.5), dpi = 140)
    
    # draw plots
    ax.fill_between(final_minus.index, 0, final_minus.values, 
                    color = 'tab:orange', label = 'minus', 
                    alpha = 0.75)
    ax.fill_between(final_plus.index, 0, final_plus.values,
                    color = 'blue', label = 'plus', alpha = 0.5)
    
    # clean up appearance
    ax.set_xlim(left, right)
    ax.set_xlabel('Genome', fontsize = fontsize)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    ax.tick_params(labelsize = fontsize)
    ax.ticklabel_format(axis="x", style="sci", scilimits=(1, 3))
    if legend:
        ax.legend(fontsize = fontsize, title = 'Strand', 
              title_fontsize = fontsize)
    if ymax is None:
        _, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax)
    
    return ax

#######################
## Metabolomic Plots ##
#######################

def met_bar_swarm(final_rates, ros_meta, k, strain_dict = None, ax = None,
              orientation = 'vertical', dodge = True, legend = False,
              legendsize = 7, labelsize = 7, ticksize = 7, 
              errwidth = 1, pointsize = 1.5, x_rot = 90):
    """
    Generates bar/swarm plots for iModulon activities
    Relies on seaborn's bar and swarm functions
    
    Parameters
    ----------
    final_rates: table of metabolomic data
    ros_meta: main metadata
    k: column in final_rates to plot
    strain_dict: dictionary with labels pointing to lists of strains
                 default: start vs. evolved
    ax: axes to plot to (will create one if not provided)
    orientation: 'horizontal' or 'vertical'
    dodge: Whether to separate each PQ level into different bars
    legend: Whether to draw the legend
    legendsize: fontsize for the legend
    labelsize: fontsize for the axis labels
    ticksize: fontsize for the tick labels on the axes
    err_width: width of errorbars in barplot
    pointsize: point size for swarm plot
    x_rot: rotation of x axis labels
    
    Returns
    -------
    ax: axes with plot
    """
    
    # organize groups of strains
    if strain_dict is None:
        strain_dict = {'Start':['0_0'],
                       'Evolved':get_ev_strains(ros_meta)}
    sample_dict = {label: final_rates.index[final_rates.my_strain.isin(strains1)] 
                   for label, strains1 in strain_dict.items()}

    # make data table for seaborn
    data = pd.DataFrame(columns = ['x_label', 'y', 'Lineage', 'PQ'])
    for label, samples in sample_dict.items():
        for s in samples:
            data.loc[s] = [label, final_rates.loc[s, k], 
                           final_rates.my_strain[s],
                           final_rates.pq[s]]
    if (data.y <= 0).all():
        data.y = -data.y
    
    # create axes if necessary
    if ax is None:
        fig, ax = plt.subplots(figsize = (0.75, 1.5),
                               dpi = 140)
    
    # define labels & columns based on orientation
    if not orientation in ['horizontal', 'vertical']:
        print('assuming vertical',
              '(else enter orientation = "horizontal")')
        orientation = 'vertical'
    if orientation == 'horizontal':
        x, y = 'y', 'x_label'
        x_label, y_label = k.replace('_', ' '), ''
        ticksize_adjust = 'x'
    elif orientation == 'vertical':
        x, y = 'x_label', 'y'
        y_label, x_label = k.replace('_', ' '), ''
        ticksize_adjust = 'y'
        
    # get PQ settings
    pq_color, _ = get_pq_dicts()
    
    # draw plots based on dodge setting
    if dodge:
        ax = sns.barplot(x = x, y = y, data = data, 
                         hue = 'PQ', palette = pq_color,
                         ax = ax, alpha = 0.7)
        ax = sns.swarmplot(x = x, y = y, data = data, 
                           hue = 'PQ', palette = pq_color,
                           s = pointsize, ax = ax, dodge = True)
    else:
        ax = sns.barplot(x = x, y = y, data = data, 
                         color = 'tab:gray',
                         ax = ax, alpha = 0.5)
        ax = sns.swarmplot(x = x, y = y, data = data, 
                           hue = 'PQ', palette = pq_color,
                           s = pointsize, ax = ax, dodge = False)
        
    # deal with the legend
    l = ax.legend(bbox_to_anchor = (1,1),loc = 'upper left', 
                  fontsize = legendsize, title = 'PQ (μM)')
    plt.setp(l.get_title(), fontsize=legendsize)
    if legend == False:
        l.remove()
    
    # ticks and labels
    ax.tick_params(labelsize = labelsize)
    ax.tick_params(ticksize_adjust, labelsize = ticksize)
    ax.set_ylabel(y_label, fontsize = labelsize)
    ax.set_xlabel(x_label, fontsize = labelsize)
    if x_rot != 0:
        ax.tick_params('x', rotation = x_rot)
    
    return ax