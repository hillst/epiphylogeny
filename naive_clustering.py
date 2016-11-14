#!/usr/bin/env python
from sklearn.preprocessing import Imputer
from sklearn.decomposition import RandomizedPCA
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import matplotlib as mlp
mlp.use('pdf')
from pylab import * 
import sys
import numpy as np


def dimensionality_reduction(data, n=100):
    """
    simplifies the dimension of the data to speed up clustering

    
    """ 


    #gotta do these.
    imp = Imputer(missing_values='NaN', strategy='mean', axis=1)
    data = imp.fit_transform(data.transpose()).transpose()

    pca = RandomizedPCA(n_components=n)
    data = pca.fit_transform(data)
    return data, pca

def evaluate_dimr(data):
    """
    Prints the cummulative explained variance of the princpal components
    """
    _, pca = dimensionality_reduction(data)
    running_total = 0
    for i, component in enumerate(sorted(pca.explained_variance_ratio_, reverse=True)):
        running_total += component
        print i, (running_total * 100) * "="


def h_clustering(data, dr=False):
    """
    performs some arbitrary heirarchical clustering algorithm

    in the presense of missing data (NaN), every NaN is assumed to be 1 distance unit away from anything.

    0, NaN = 1
    NaN, NaN = 1
    1, NaN = 1
    0, 1 = 1

    Then everything is normalized by teh size of the vector. Sites with many missing data will be farther away from everything *uniformly*. 

    If we are using pca, we should be using euclidean distance.
    """
    if dr == True:
        metric = 'euclidean'
    else:
        metric = 'hamming'
    dists = pdist(data, metric ) #how does hammign distance handle mising data?
    links = linkage(dists, 'ward')
    return links

def build_newick(clustering_tree):
    """
    Converts our heirarchical clustering into a newick tree
    Get this from dendogram probably.
    """
    ()
    pass

def plot_dendogram(links):
    figure(figsize=(19, 10))
    dendrogram(
        links,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    savefig("dendogram.pdf")


def permute_data(data):
    """
    applies random noise to the data to simulate a randomized algorithm to check stability of the data
    """
    pass

def read_csv(filename):
    names, site, data = [],[],[]
    for line in open(filename):
        if "RRBS" in line or "#" in line:
            names = line.strip().split(",")
        else:
            line_arr = line.strip().split(",")
            site = line_arr[0]
            row = []
            for d in line_arr[1:]:
                if d == "-":
                    row.append(np.NaN)
                else:
                    row.append(float(d))
            data.append(row)
    data = np.asarray(data) 
    #we need to transpose for the imputer to work correctly, then transpose back for the rest.
    #data = imp.fit_transform(data).transpose()
    print data
    return names, data #.transpose()

def main():
    names, data = read_csv(sys.argv[1])
    #evaluate_dr(data)
    print data
    data, pca = dimensionality_reduction(data, n=5)
    print data
    model = h_clustering(data, dr=True)
    print model
    plot_dendogram(model)
    #print model.children_
    #print len(names) 
    #print sorted(dict(enumerate(model.children_, model.n_leaves_)).keys())
    
    
    

if __name__ == "__main__":
    main()
    
