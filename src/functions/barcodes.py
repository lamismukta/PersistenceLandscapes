import numpy as np
from ripser import Rips
import gudhi as gd
import urllib

# Import point cloud data for a protein
def get_data(pdbid,chain):
    url = 'https://knotprot.cent.uw.edu.pl/chains/{pdbid}/{chain}/chain.xyz.txt'
    link = url.format(pdbid = pdbid,chain = chain)
    file = urllib.request.urlopen(link)
    g = np.genfromtxt(file, delimiter = '')
    data = g[0:,1:]
    return data

# Plot a ripser diagram from data
def ripser_plot(data, title):
    rips = Rips(maxdim=2)
    diagrams = rips.fit_transform(data)
    rips.plot(diagrams, title = title)

# Plot ph diagram for a protein with ripser
# Formerly rips_pd
def import_and_plot(pdbid,chain):
    data = get_data(pdbid, chain)
    plot = ripser_plot(data, title = 'Persistence Diagram for {pdbid}_{chain}'.format(pdbid=pdbid,chain=chain))
    return plot

# Interpolate k points in a point cloud
# Formerly interp
def interpolate(data,k):
    l = len(data)
    v = np.empty((l+(l-1)*k,3))
    for j in range(0,k+1):
        for i in range(0,l-1):
            v[(k+1)*i+j] = data[i] + (j/(k+1))*(data[i+1]-data[i])
    v[(k+1)*(l-1)]=data[l-1]
    return v

# Gudhi barcode from a point cloud
# Formerly gd_barcode
# Reference https://github.com/GUDHI/TDA-tutorial/blob/master/Tuto-GUDHI-persistence-diagrams.ipynb
def gudhi_barcode(data, max_edge_length=15):
    skeleton = gd.RipsComplex(data, max_edge_length=max_edge_length)
    rips_simplex_tree = skeleton.create_simplex_tree(max_dimension = 3)
    return rips_simplex_tree.persistence()

# Bottleneck distance between barcodes
# Formerly bn
def bottleneck(barcode_1, barcode_2):
    diag1 = [item[1] for item in barcode_1]
    diag2 = [item[1] for item in barcode_2]
    return gd.bottleneck_distance(diag1, diag2)

# Bottleneck separating classes
# Formerly my_bn
def bottleneck_by_dimension(barcode_1,barcode_2):
    diag1_0 = [item[1] for item in barcode_1 if item[0]==0]
    diag2_0 = [item[1] for item in barcode_2 if item[0]==0]
    diag1_1 = [item[1] for item in barcode_1 if item[0]==1]
    diag2_1 = [item[1] for item in barcode_2 if item[0]==1]
    diag1_2 = [item[1] for item in barcode_1 if item[0]==2]
    diag2_2 = [item[1] for item in barcode_2 if item[0]==2]
    d0 = gd.bottleneck_distance(diag1_0, diag2_0)
    d1 = gd.bottleneck_distance(diag1_1, diag2_1)
    d2 = gd.bottleneck_distance(diag1_2, diag2_2)
    return max(d0,d1,d2)
