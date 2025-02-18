import numpy as np
import conga
import pandas as pd
import matplotlib.pyplot as plt

# call the Anndata file and conga on
gex_datafile = "/home/<user>/analysis/tcr/inf.pb_mtx/"
gex_datatype = '10x_mtx' # other possibilities right now: ['10x_mtx', 'h5ad'] (h5ad from scanpy)
bcr_datafile = '/home/<user>/analysis/tcr/conga_infant_metadata_pb.csv'
metadata_file = "/home/<user>/analysis/tcr/infant_metadata_conga_all.csv"
organism = 'human_ig'
clones_file = 'infant_bcell.pb_clones.tsv'
outfile_prefix = 'inf_bcell.pb'

# this creates the BCRdist 'clones file'
conga.tcrdist.make_10x_clones_file.make_10x_clones_file_batch( bcr_datafile, organism, clones_file )

# this is necessary for the BCR data (see https://github.com/phbradley/conga/tree/master?tab=readme-ov-file#human-melanoma-b-cell-dataset)
oldfile = clones_file
clones_file = clones_file[:-4]+'_condensed.tsv'
input_distfile = clones_file[:-4]+'_AB.dist'
conga.preprocess.condense_clones_file_and_barcode_mapping_file_by_tcrdist(
            oldfile, clones_file, 50, organism,
            output_distfile=input_distfile)

# this command will create another file with the kernel PCs for subsequent reading by conga
conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file( clones_file, organism )

adata = conga.preprocess.read_dataset(gex_datafile, gex_datatype, clones_file )

# store the organism info in adata
adata.uns['organism'] = organism

adata = conga.preprocess.filter_and_scale( adata )

adata = conga.preprocess.reduce_to_single_cell_per_clone( adata )
adata = conga.preprocess.cluster_and_tsne_and_umap( adata )

plt.figure(figsize=(12,6))
plt.subplot(121)
xy = adata.obsm['X_gex_2d']
clusters = np.array(adata.obs['clusters_gex'])
cmap = plt.get_cmap('tab20')
colors = [ cmap.colors[x] for x in clusters]
plt.scatter( xy[:,0], xy[:,1], c=colors)
plt.title('GEX UMAP colored by GEX clusters')

plt.subplot(122)
xy = adata.obsm['X_tcr_2d']
clusters = np.array(adata.obs['clusters_tcr'])
cmap = plt.get_cmap('tab20')
colors = [ cmap.colors[x] for x in clusters]
plt.scatter( xy[:,0], xy[:,1], c=colors)
plt.title('BCR UMAP colored by BCR clusters')

plt.savefig("UMAP_"+outfile_prefix+"_GEX.BCR_clusters.pdf", dpi=500)
plt.cla()

# these are the nbrhood sizes, as a fraction of the entire dataset:
nbr_fracs = [0.01, 0.1]

# we use this nbrhood size for computing the nndists
nbr_frac_for_nndists = 0.01

all_nbrs, nndists_gex, nndists_tcr = conga.preprocess.calc_nbrs(
    adata, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)

# stash these in obs array, they are used in a few places...                                                                                                                
adata.obs['nndists_gex'] = nndists_gex
adata.obs['nndists_tcr'] = nndists_tcr

conga.preprocess.setup_tcr_cluster_names(adata) #stores in adata.uns

results = conga.correlations.run_graph_vs_graph(adata, all_nbrs)

conga_scores = adata.obs['conga_scores']

good_mask = ( conga_scores <= 1.0 )
adata.obs['good_score_mask'] = good_mask

print(f'found {np.sum(good_mask)} conga hits')

# write the results to a tsv file
clusters_gex = np.array(adata.obs['clusters_gex'])
clusters_tcr = np.array(adata.obs['clusters_tcr'])

indices = results['clone_index']
results['gex_cluster'] = clusters_gex[indices]
results['bcr_cluster'] = clusters_tcr[indices]
for tag in 'va ja cdr3a vb jb cdr3b'.split():
    results[tag] = list(adata.obs[tag][indices])
tsvfile = outfile_prefix+'_graph_vs_graph_hits.tsv'
print('saving graph-vs-graph results to file:',tsvfile)

results.to_csv(tsvfile, sep='\t', index=False)

#put the conga hits on top
colors = np.sqrt(np.maximum(-1*np.log10(conga_scores),0.0))
reorder = np.argsort(colors)

plt.figure(figsize=(12,6))
plt.subplot(121)
xy = adata.obsm['X_gex_2d']
plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder], vmin=0, vmax=np.sqrt(5))
plt.title('GEX UMAP colored by conga score')

plt.subplot(122)
xy = adata.obsm['X_tcr_2d']
plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder], vmin=0, vmax=np.sqrt(5))
plt.title('BCR UMAP colored by conga score')
plt.savefig("UMAP_"+outfile_prefix+"_GEX.BCR_congascore.pdf", dpi=500)

nbrs_gex, nbrs_bcr = all_nbrs[0.1]

min_cluster_size = 5

# Generate graph-vs-graph results plot
# If batch keys have been specified and database matching results are stored in the object,
#  they will appear in the plot by default.
conga.plotting.make_graph_vs_graph_logos(
    adata, outfile_prefix, min_cluster_size,
    nbrs_gex,
    nbrs_bcr,
)

# run the graph-vs-features analysis and store results
#  in adata.uns[‘conga_results’]
conga.correlations.run_graph_vs_features(
    adata, all_nbrs, outfile_prefix = outfile_prefix)

# generate plots visualizing the results
conga.plotting.make_graph_vs_features_plots(
    adata, all_nbrs, outfile_prefix,
    clustermap_max_type_features=25,
)

# run an alternative graph-vs-features analysis
#  using the Yosef lab’s HotSpot algoritmh
conga.correlations.find_hotspots_wrapper(
    adata, all_nbrs, nbr_fracs = [0.1],
    outfile_prefix = outfile_prefix,
)

# generate plots visualizing the results
conga.plotting.make_hotspot_plots(
    adata, all_nbrs,
    outfile_prefix = outfile_prefix,
    clustermap_max_type_features=25,
)

# Run BCR clumping analysis
conga.tcr_clumping.assess_tcr_clumping(
    adata, outfile_prefix = outfile_prefix,
)

# make plots to visualize the results
conga.plotting.make_tcr_clumping_plots(
    adata, nbrs_gex, nbrs_bcr, outfile_prefix, 
    min_cluster_size_for_logos=min_cluster_size,
)

# calc bcr sequence features of good cluster pairs
good_bicluster_bcr_scores = conga.correlations.calc_good_cluster_tcr_features(
    adata, good_mask, clusters_gex, clusters_tcr, conga.tcr_scoring.all_tcr_scorenames, verbose=False,
    min_count=min_cluster_size)

# run rank_genes on most common biclusters
rank_genes_uns_tag = 'rank_genes_good_biclusters'
conga.correlations.run_rank_genes_on_good_biclusters(
    adata, good_mask, clusters_gex, clusters_tcr, min_count=min_cluster_size, key_added= rank_genes_uns_tag)

gex_header_bcr_score_names = ['imhc', 'cdr3len', 'mait', 'nndists_tcr']

logo_pngfile = outfile_prefix+'_bicluster_logos.pdf'

conga.plotting.make_logo_plots(
    adata, nbrs_gex, nbrs_bcr, min_cluster_size, logo_pngfile,
    good_bicluster_tcr_scores=good_bicluster_bcr_scores,
    rank_genes_uns_tag = rank_genes_uns_tag,
    gex_header_tcr_score_names = gex_header_bcr_score_names )

pval_threshold = 0.05
results = []
for nbr_frac in nbr_fracs:
    nbrs_gex, nbrs_bcr = all_nbrs[nbr_frac]
    print('finding biased GEX features for nbrhoods with size', nbr_frac, nbrs_gex.shape)
    results.append( conga.correlations.tcr_nbrhood_rank_genes_fast( adata, nbrs_bcr, pval_threshold, verbose=False))
    results[-1]['nbr_frac'] = nbr_frac

# save the results to a file
tsvfile = outfile_prefix+'_bcr_nbr_graph_vs_gex_features.tsv'
print('making:', tsvfile)
results_df = pd.concat(results, ignore_index=True)
results_df.to_csv(tsvfile, index=False, sep='\t')

pngfile = outfile_prefix+'_bcr_nbr_graph_vs_gex_features.pdf'
print('making:', pngfile)
conga.plotting.plot_ranked_strings_on_cells(
    adata, results_df, 'X_tcr_2d', 'clone_index', 'mwu_pvalue_adj', 1.0, 'feature', pngfile)

pngfile = outfile_prefix+'_bcr_nbr_graph_vs_gex_features_panels.pdf'
print('making:', pngfile)
conga.plotting.make_feature_panel_plots(adata, 'bcr', all_nbrs, results_df, pngfile, 'gex')

pval_threshold = 0.05
results = []
bcr_score_names = conga.tcr_scoring.all_tcr_scorenames # the BCR features to use
for nbr_frac in nbr_fracs:
    nbrs_gex, nbrs_bcr = all_nbrs[nbr_frac]
    results.append( conga.correlations.gex_nbrhood_rank_tcr_scores(adata, nbrs_gex, bcr_score_names, pval_threshold, verbose=False ))
    results[-1]['nbr_frac'] = nbr_frac
results_df = pd.concat(results, ignore_index=True)

tsvfile = outfile_prefix+'_gex_nbr_graph_vs_bcr_features.tsv'
print('making:', tsvfile)
results_df.to_csv(tsvfile, index=False, sep='\t')

pngfile = outfile_prefix+'_gex_nbr_graph_vs_bcr_features.pdf'
print('making:', pngfile)

conga.plotting.plot_ranked_strings_on_cells(
    adata, results_df, 'X_gex_2d', 'clone_index', 'mwu_pvalue_adj', 1.0, 'feature', pngfile,
    direction_column='ttest_stat')

pngfile = outfile_prefix+'_gex_nbr_graph_vs_bcr_features_panels.pdf'
print('making:', pngfile)
conga.plotting.make_feature_panel_plots(adata, 'gex', all_nbrs, results_df, pngfile, 'tcr')
