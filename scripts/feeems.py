
# Running FEEMS (Fast Estimation of Effective Migration Surfaces)

# Juliette Archambeau - September 2022


# Docs 
######
# Marcus et al. (2019): https://elifesciences.org/articles/61927
# Associated tutorial: https://github.com/NovembreLab/feems/blob/main/docsrc/notebooks/getting-started.ipynb
# Building the outer polygon: https://www.keene.edu/campus/maps/tool/



# In this script, I run FEEMS on the 9,817 polymorphic SNPs
# Snps were genotyped with the Illumina Infinium assay and the Affymetrix assay
# The genomic files in plink format were generated in the R script 'createinputFilesFEEMs.R'

# Input files are generated in the R script 'GeneratingInputFilesFEEMs.R'


# Import the required packages and feems
########################################
# base
import os
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink


# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"



# Working directory
###################
os.chdir('/home/juliette/Documents/H2Pinpin/H2Pinpin/')



# Import the input files
########################
(bim, fam, G) = read_plink("{}/GenotypeMatrixPlink".format("data/FEEMs/")) # genomic matrix (5,165 SNPs x 523 genotypes)
G = np.array(G).T
print("n_samples={}, n_snps={}".format(G.shape[0], G.shape[1]))

imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)


#coord = np.loadtxt("{}/pop_coord.coord".format("data/FEEMs/"))  # population coordinates
coord = np.loadtxt("{}/clon_coord.coord".format("data/FEEMs/"))  # clones coordinates
coord[:5,]

outer = np.loadtxt("{}/OuterPolygonMaritimePinePopulations.outer".format("data/FEEMs/")) # coordinates of the outer polygon
outer[:5,]

grid_path = "{}/grid_100.shp".format("/home/juliette/Documents/H2Pinpin/H2Pinpin/data/FEEMs") # path to discrete global grid
# I have to build my own discrete global grid with a higher resolution (using the dggridR package in R

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)

edges[:5,]
grid[:5,]



# Map of the sample coordinates, graph and observed nodes
#########################################################

sp_graph = SpatialGraph(G, coord, grid, edges, scale_snps=True)

projection = ccrs.EquidistantConic(central_longitude=0.1647949, central_latitude=40.4470470)                                   
                                             
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False) # necessary that version networkx==2.4.0 is installed
v.draw_obs_nodes(use_ids=False)

                                      
plt.savefig('/home/juliette/Documents/H2Pinpin/H2Pinpin/figs/FEEMs/grid100.png')                                          
                                             
# Fit feems with a smooth regularization
########################################
sp_graph.fit(lamb = 20.0)        


fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10,cbar_loc="lower right")
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar() # I would like to move the legend to the lower right but I do not know how to do it!

plt.savefig('/home/juliette/Documents/H2Pinpin/H2Pinpin/figs/FEEMs/fitted_feems_lamb20.png') 



# Fit feems with a less smooth regularization
#############################################
sp_graph.fit(lamb = 2.0)

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10,cbar_loc="lower right")
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar() # I would like to move the legend to the lower right but I do not know how to do it!


plt.savefig('/home/juliette/Documents/H2Pinpin/H2Pinpin/figs/FEEMs/fitted_feems_lamb2.png') 

             
