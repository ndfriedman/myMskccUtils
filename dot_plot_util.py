#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

import matplotlib
matplotlib.use('Agg') #the order of this import crap is cirtical for avoiding a breaking error
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#--------------------

#TAKEN From: http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap

#-------------------------

def dot_plot_with_words(l):
	fig, ax = plt.subplots()
	for entry in l:
		sense, specif, name = entry
		ax.text(sense, specif, name, ha='center', size=20)
	plt.savefig('noahTest')
	sys.exit()
	plt.clf()

def scatter_plot_with_legend(l, title):

	def make_scatter_plot_legend(legendInfo):
		scatInfo = [x[1] for x in legendInfo]
		names = [x[0] for x in legendInfo]
		plt.legend(scatInfo,
           names,
           scatterpoints=1,
           loc='lower left',
           ncol=3,
           fontsize=8)

	colors = cm.rainbow(np.linspace(0, 1, len(l)))
	cntr = 0
	legendL = []
	for entry in l:
		sense, specif, name = entry
		curColor = colors[cntr]
		scat = plt.scatter(sense, specif, s=100, c=curColor)
		cntr += 1
		legendL.append((name, scat))

	make_scatter_plot_legend(legendL)
	plt.xlabel('Sensitivity')
	plt.ylabel('Specificity')
	plt.xlim([0,1])
	#plt.ylim([.8,1])
	plt.title(title)

	filename = 'noahTest'
	print 'saving figure to ', filename
	plt.savefig(filename)
	sys.exit()

def scatter_plot_line(vals, title):
	senses = [x[0] for x in vals]
	specifs = [x[1] for x in vals]

	for i in range(len(senses)):
		plt.scatter(senses[i], specifs[i], s=50) #overlay colored dots with heatmapesque colors

	plt.title(title)
	plt.xlabel('Specificity')
	plt.ylabel('Sensitivity')
	plt.plot(senses, specifs)


	filename = 'signatureConfidenceCapture' + title
	print 'saving figure to ', filename
	plt.savefig(filename)
	plt.clf()

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')



	args = parser.parse_args()

if __name__ == '__main__':
    main()















