import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

# taken from https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb
os.makedirs('tmp')
os.system('wget -O tmp/pbmc8k_filtered_gene_bc_matrices.tar.gz http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz')
os.system('cd tmp; tar xfz pbmc8k_filtered_gene_bc_matrices.tar.gz')

input_dir = 'tmp/filtered_gene_bc_matrices/GRCh38/'
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()

# run scrublet and save results to compare it with R version
output_dir = "tmp/filtered_gene_bc_matrices/out/"
os.makedirs(output_dir)

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
# 16.5 secs

np.savetxt(output_dir+"doublet_scores.csv", doublet_scores, delimiter=",")                                                        
np.savetxt(output_dir+"doublet_scores_sim.csv", scrub.doublet_scores_sim_, delimiter=",")                                 
np.savetxt(output_dir+"pca.csv", scrub.manifold_obs_, delimiter=",")
np.savetxt(output_dir+"pca.sim.csv", scrub.manifold_sim_, delimiter=",")
np.savetxt(output_dir+"doublet_parents_.csv", scrub.doublet_parents_, delimiter=",") # to make result comparable I'll use same doublets parent in R
np.savetxt(output_dir+"_gene_filter.csv", scrub._gene_filter, delimiter=",")
np.savetxt(output_dir+"predicted_doublets.csv", scrub.predicted_doublets_, delimiter=",")


# run exact (non approx) KNN to make results not stochastic
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30,
                                                          use_approx_neighbors=False)
#  57 secs                                                     
np.savetxt(output_dir+"doublet_scores_eknn.csv", doublet_scores, delimiter=",")                                                        
np.savetxt(output_dir+"doublet_scores_sim_eknn.csv", scrub.doublet_scores_sim_, delimiter=",")  

# run "run.scrublet.R" and compare results
