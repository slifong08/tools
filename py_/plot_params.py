#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import seaborn as sns
import subprocess


"""
### SNS PALETTE EXAMPLE ###
#colors = [ "amber", "dusty purple", "windows blue"]
#PAL = sns.xkcd_palette(colors)
#sns.palplot(PAL)

### custom color palette
yg = '#f7fcb9'
stemg = '#addd8e'
kellyg = '#31a354'

y = mcolors.to_rgba(yg)
lg = mcolors.to_rgba(stemg)
g = mcolors.to_rgba(kellyg)


colors = [g, lg, y]  # kelly (core) -> light (derived) - > yellow (simple)
cmap_name = 'yg_'
CM = LinearSegmentedColormap.from_list(cmap_name, colors, N=3)
print("COLORMAP - CM", CM)
"""


# In[2]:


def summary_labels():
    summary_labels = {
                        "1111":'ConsAct',
                        "1000":'HuAct',
                        "1100":'HuDNA',
                        "1010":'HuEnv',
                        "1001":'SpSpAct',
                        "1101":None,
                        "1011":None,
                        "1110":None, 
                    }
    return summary_labels


# In[3]:


def fonts():
    font_fam = matplotlib.rcParams['font.family'] = "sans-serif"
    font = matplotlib.rcParams['font.sans-serif'] = "Arial"
    font_size = matplotlib.rcParams['font.size'] = 18
    return font_fam, font, font_size


# In[4]:


def annotate_bar(graph, ax):

    if ax == "x":

        for p in graph.patches:

            graph.annotate(round(p.get_height(), 1), 
                            (p.get_x() + p.get_width() / 2.0, 
                            p.get_height()), 
                            ha = 'center', 
                            va = 'center', 
                            xytext = (0, 5),
                            textcoords = 'offset points')
    else:
        for p in graph.patches:
            print(p)

            graph.annotate(round(p.get_width(), 1), 
                            (p.get_y() + p.get_height() / 2.0, 
                             p.get_width()), 
                             ha = 'center', 
                             va = 'center', 
                             xytext = (0, 5),
                             textcoords = 'offset points')


# In[5]:


def plot_pie(data, label, size, title):

    if type(size) is list: 
        sizes = size
        labels = label
    else:
        pie_data = data.groupby(label)[size].count().reset_index()
        print(pie_data,"\n\nsum", pie_data.sum())
        
        labels = pie_data[label]

        sizes = pie_data[size]
    #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, 
            #explode=explode, 
            labels=labels, autopct='%1.1f%%',
            startangle=90)
    ax1.set(title = title)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    #plt.show()


# In[6]:


def plot_OR_bar(x, y, data, ci_upper, ci_lower, ytick_multiple, ylim, title):
    fig, ax = plt.subplots(figsize=(6,6))
    
    # get difference between y and upper/lower CI for plotting. 
    l = np.array(data[y]-data[ci_lower])
    u = np.array(data[ci_upper]- data[y])
             
    yerr = [l, u]

    sns.barplot(
            x=x, y=y, data=data,
            hue = hue,
            linewidth=2.5, 
            #facecolor=(1, 1, 1, 0),
            edgecolor=".2",
            yerr=yerr
            )
   
    plt.axhline(0, color = "grey", linewidth = 2.5)  # plot a line at zero


    #ax.set_xticklabels(["Complex\nEnhancer", "Simple\nEnhancer", "Complex\nEnhancer\nCore v. Derived"])

    ax.set_xlabel("")

    #ax.get_yaxis().ticker.LogLocator(base=2)

    # set the y ticks
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(2**x))
    ax.yaxis.set_major_formatter(ticks)
    ax.yaxis.set_major_locator(MultipleLocator(ytick_multiple))
    ax.set(
        xlabel = "",
        ylabel="Odds ratio, log2-scaled",  # y label
        title = title,
        ylim=ylim)
    return fig, ax





