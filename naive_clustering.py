#!/usr/bin/env python
from sklearn.preprocessing import Imputer
from sklearn.decomposition import RandomizedPCA
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet, to_tree
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
    data = data.transpose() 
    data = imp.fit_transform(data)
    data = data.transpose()

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
    lens = set()
    for i, line in enumerate(open(filename)):
        if "RRBS" in line or "#" in line:
            names = line.strip().split(",")[1:]
        else:
            line_arr = line.strip().split(",")
            site = line_arr[0]
            row = []
            for d in line_arr[1:]:
                if d == "-" or d == "2":
                    row.append(np.NaN)
                else:
                    row.append(float(d))
            if len(row) == 103:
                print row
                print i
            if row.count(np.NaN) < 93:
                lens.add(len(row))
                data.append(row)
    data = np.asarray(data) 
    print lens
    #we need to transpose for the imputer to work correctly, then transpose back for the rest.
    #data = imp.fit_transform(data).transpose()
    return names, data #.transpose()

def get_newick(node, newick, parentdist, leaf_names):
    """
    http://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format

    tree = hierarchy.to_tree(Z,False)
    get_newick(tree, "", tree.dist, leaf_names)
    """
    if node.is_leaf():
        if node.id < len(leaf_names):
        #return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
            name = leaf_names[node.id]
        else:
            name = node.id
        return "%s:%.2f%s" % (name , parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def main():
    names, data = read_csv(sys.argv[1])
    #evaluate_dr(data)
    data, pca = dimensionality_reduction(data, n=25)
    print "clustering"
    import random
    print len(data)
    data = np.asarray(random.sample(data, 10000))
    model = h_clustering(data, dr=True)
    tree = to_tree(model, False)
    print get_newick(tree, "", tree.dist, names)

    #plot_dendogram(model)
    #print model.children_
    #print len(names) 
    #print sorted(dict(enumerate(model.children_, model.n_leaves_)).keys())
    
    
    

if __name__ == "__main__":
    main()
    
