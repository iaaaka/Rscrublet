# run run.scrublet.py first
# install external dependences
# install.packages(c('Matrix','RcppAnnoy','RSpectra')) 
# load R version of scrublet (is it not wrapped into package yet)
source('scrublet.R')

# load data as sparse matrix
m = t(readMM('tmp/filtered_gene_bc_matrices/GRCh38/matrix.mtx'))
rownames(m) = read.table('tmp/filtered_gene_bc_matrices/GRCh38/barcodes.tsv')[,1]
colnames(m) = read.table('tmp/filtered_gene_bc_matrices/GRCh38/genes.tsv')[,1]

# run doublets calls #########
scrr = scrub_doublets(E_obs = m,
										 expected_doublet_rate=0.06,
										 min_counts=2, 
										 min_cells=3, 
										 min_gene_variability_pctl=85, 
										 n_prin_comps=30)
# bit longer than python version: 31 vs 16 secs
# set threshold and count stat
scrr=call_doublets(scrr)
# plot score hists
plot_doublet_histogram(scrr)
# comp to python
dsoa = as.numeric(readLines('tmp/filtered_gene_bc_matrices/out/doublet_scores.csv'))
dssa = as.numeric(readLines('tmp/filtered_gene_bc_matrices/out/doublet_scores_sim.csv'))

plot(dsoa,scrr$doublet_scores_obs,bty='n',xlab='Python',ylab='R')
plot(dssa,scrr$doublet_scores_sim,bty='n',xlab='Python',ylab='R') # since I used independent simulated doublets parents scores for simulated doublets are different
# compare calls
pyd = as.numeric(readLines('tmp/filtered_gene_bc_matrices/out/predicted_doublets.csv'))
table(python=pyd,R=scrr$doublet_scores_obs>scrr$threshold)

# now use same simulated doublets as in python #######
doublets_parents = as.matrix(read.csv('tmp/filtered_gene_bc_matrices/out/doublet_parents_.csv',header = F))+1 # same doublet parents as was usedin python version to compare with

scrrp = scrub_doublets(E_obs = m,
											expected_doublet_rate=0.06,
											min_counts=2, 
											min_cells=3, 
											min_gene_variability_pctl=85, 
											n_prin_comps=30,
											doublets_parents=doublets_parents)
scrrp=call_doublets(scrrp)
plot_doublet_histogram(scrrp)

# comp to python
plot(dsoa,scrrp$doublet_scores_obs,bty='n',xlab='Python',ylab='R') # looks bit more similar
plot(dssa,scrrp$doublet_scores_sim,bty='n',xlab='Python',ylab='R') # now scores for simulated doublets also correlates
table(python=pyd,R=scrrp$doublet_scores_obs>scrrp$threshold)

# there are still some differences because of stochasticity in approximate KNN, lets try exact KNN
# current realization of exact KNN is bit slow, it takes about 5 mins
scrrpe = scrub_doublets(E_obs = m,
										 use_approx_neighbors = FALSE,
										 expected_doublet_rate=0.06,
										 min_counts=2, 
										 min_cells=3, 
										 min_gene_variability_pctl=85, 
										 n_prin_comps=30,
										 doublets_parents=doublets_parents)
scrrpe=call_doublets(scrrpe)
plot_doublet_histogram(scrrpe)

dsoe = as.numeric(readLines('tmp/filtered_gene_bc_matrices/out/doublet_scores_eknn.csv'))
dsse = as.numeric(readLines('tmp/filtered_gene_bc_matrices/out/doublet_scores_sim_eknn.csv'))
plot(dsoe,scrrpe$doublet_scores_obs,bty='n',xlab='Python',ylab='R')
plot(dsse,scrrpe$doublet_scores_sim,bty='n',xlab='Python',ylab='R') 
table(dsoe-scrrpe$doublet_scores_obs)
table(dsse-scrrpe$doublet_scores_sim)
# so, now scores are identical

# TODO
# 1) wrap into r-package (S4? R5? seurat?)
# 2) check whether is can be accelerate to match python performance
# 3) Add functionalities:
#   a) other embeddings
#   b) accelerate exact knn (do we need it?)
#   c) smth else?
# 4) check
#   a) test with other examples/parameters
#   b) check fit_percentile in get_vscores - looks like there is a mistake in python version
# 5) justify usage of svn instead of pca (looks correct...)


