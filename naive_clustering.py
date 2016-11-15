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


def build_colorful_tree(newick, filename=""):
    """
    Note that these will fail if we dont have all the pre-reqs and it is not triival to get them all.
    This stuff is NOT general purpose.
    """
    from ete3 import Tree, TreeStyle, CircleFace, TextFace
    tree = Tree(newick)

    #setup colors and treestyle
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -180 # 0 degrees = 3 o'clock
    ts.force_topology = True
    ts.arc_span = 360

    face = CircleFace(30, "MediumSeaGreen")
    face.margin_top = 1000
    ts.legend.add_face(face, column=0)
    face = TextFace("Normal B-cell", fsize=64)
    face.margin_right = 100
    face.margin_top = 1000
    ts.legend.add_face(face, column=1)

    ts.legend.add_face(CircleFace(30, "SeaGreen"), column=0)
    face = TextFace("Normal B CD19pcell", fsize=64)
    face.margin_right = 100
    ts.legend.add_face(face, column=1)

    ts.legend.add_face(CircleFace(30, "ForestGreen"), column=0)
    face = TextFace("Normal B CD19pCD27pcell", fsize=64)
    face.margin_right = 100
    ts.legend.add_face(face, column=1)

    ts.legend.add_face(CircleFace(30, "Green"), column=0)
    face = TextFace("Normal B CD19pCD27mcell", fsize=64)
    face.margin_right = 100
    ts.legend.add_face(face, column=1)

    ts.legend.add_face(CircleFace(30, "RoyalBlue"), column=0)
    face = TextFace("CLL all-batches", fsize=64)
    face.margin_right = 100
    ts.legend.add_face(face, column=1)

    #draw tree
    from ete3 import NodeStyle
    styles= {}
    styles["normal_B"] = NodeStyle( bgcolor="MediumSeaGreen", hz_line_color="Black", vt_line_color="Black")
    styles["NormalBCD19pcell"] = NodeStyle( bgcolor="SeaGreen", hz_line_color="Black", vt_line_color="Black")
    styles["NormalBCD19pCD27pcell"] = NodeStyle( bgcolor="ForestGreen", hz_line_color="Black", vt_line_color="Black")
    styles["NormalBCD19pCD27mcell"] = NodeStyle( bgcolor="Green", hz_line_color="Black", vt_line_color="Black")
    styles["CLL"] = NodeStyle( bgcolor="RoyalBlue", hz_line_color="Black", vt_line_color="Black" )

    for node in tree.traverse("postorder"):
        #print node.set_style()
        if len(node.get_leaf_names()) == 1:
            name = node.get_leaf_names()[0]
            if "normal_B" in name:
                node.set_style(styles["normal_B"])
            elif "NormalBCD19pcell" in name:
                node.set_style(styles["NormalBCD19pcell"])

            elif "NormalBCD19pCD27pcell" in name:
                node.set_style(styles["NormalBCD19pCD27pcell"])

            elif "NormalBCD19pCD27mcell" in name:
                node.set_style(styles["NormalBCD19pCD27mcell"])
            else:
                node.set_style(styles["CLL"])
    #lol
    tree.render(filename, w=10,dpi=600, units='in',tree_style=ts)
        
        

def main():
    names, data = read_csv(sys.argv[1])
    #evaluate_dr(data)
    data, pca = dimensionality_reduction(data, n=25)
    print "clustering"
    import random
    print len(data)
    #data = np.asarray(random.sample(data, 10000))
    data = np.asarray(data)
    model = h_clustering(data, dr=True)
    tree = to_tree(model, False)
    nw =  get_newick(tree, "", tree.dist, names)
    print nw
    build_colorful_tree(nw) 

    #plot_dendogram(model)
    #print model.children_
    #print len(names) 
    #print sorted(dict(enumerate(model.children_, model.n_leaves_)).keys())
    
    
    

if __name__ == "__main__":
    main()
    
