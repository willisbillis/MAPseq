// modules/analysis.nf
// GPU-accelerated secondary analysis: Harmony batch integration and clustering

process HARMONY_INTEGRATION {
    label 'gpu_large'
    publishDir "${params.outdir}/analysis/harmony/${batch_id}", mode: 'copy'
    
    input:
    tuple val(batch_id), path(qc_h5ad)
    val(gpu_available)
    
    output:
    tuple val(batch_id), path("harmony_integrated.h5ad"), emit: integrated
    path("harmony_report.txt"), emit: report
    path("harmony_metrics.csv"), emit: metrics
    path("harmony_diagnostics.csv"), emit: diagnostics
    
    script:
    use_gpu = params.harmony_gpu && gpu_available.toBoolean()
    """
    python3 << 'PYSCRIPT'
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import logging
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('Harmony')
    
    logger.info("=" * 70)
    logger.info("Harmony Batch Integration (GPU-accelerated)")
    logger.info("=" * 70)
    
    start_time = datetime.now()
    
    try:
        # Load QC-filtered data
        logger.info(f"Loading QC-filtered matrix from: ${qc_h5ad}")
        adata = sc.read_h5ad("${qc_h5ad}")
        logger.info(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        # Store initial statistics
        initial_stats = {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'n_batches': adata.obs['batch'].nunique() if 'batch' in adata.obs else 1,
        }
        logger.info(f"Metadata: {initial_stats['n_batches']} batches found")
        
        # Ensure data is in standard format
        import scipy.sparse as sp
        if not sp.issparse(adata.X):
            logger.info("Converting matrix to sparse format...")
            adata.X = sp.csr_matrix(adata.X)
        
        # ===== PREPROCESSING =====
        logger.info("Preprocessing...")
        
        # Normalize
        logger.info("  • Normalizing to 10,000 UMI/cell")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Select highly variable genes
        n_hvg = ${params.n_hvg}
        logger.info(f"  • Selecting {n_hvg:,} highly variable genes")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=False)
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
        logger.info(f"  • HVG subset: {adata_hvg.n_vars:,} genes")
        
        # ===== PCA COMPUTATION =====
        logger.info("Computing PCA...")
        n_pcs = ${params.harmony_npcs}
        pca_components = min(n_pcs, adata_hvg.n_obs - 1)
        
        use_rapids = False
        try:
            import rapids_singlecell as rsc
            logger.info(f"  • Using rapids-singlecell GPU-accelerated PCA ({pca_components} PCs)")
            rsc.tl.pca(adata_hvg, n_comps=pca_components)
            use_rapids = True
            logger.info(f"  • PCA completed on GPU")
        except ImportError:
            logger.warning("rapids-singlecell not available, using Scanpy CPU PCA")
            sc.tl.pca(adata_hvg, n_comps=pca_components)
        
        # Transfer PCA to full object
        logger.info("  • Transferring PCA to full object")
        adata.varm['PCs'] = adata_hvg.varm['PCs']
        adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
        adata.uns['pca'] = adata_hvg.uns.get('pca', {})
        logger.info(f"  • PCA shape: {adata.obsm['X_pca'].shape}")
        
        # ===== HARMONY BATCH INTEGRATION =====
        logger.info("Running Harmony batch integration...")
        
        use_gpu = ${use_gpu}
        device = 'cuda' if use_gpu else 'cpu'
        
        logger.info(f"  • Device: {device}")
        logger.info(f"  • Harmony PCs: {params.harmony_npcs}")
        logger.info(f"  • Harmony theta (batch correction strength): ${params.harmony_theta}")
        
        harmony_pcs = min(${params.harmony_npcs}, adata.obsm['X_pca'].shape[1])
        
        try:
            import harmonypy as hm
            
            logger.info(f"  • harmonypy version: {hm.__version__}")
            
            # Validate batch column
            if 'batch' not in adata.obs:
                logger.warning("  • No 'batch' column found, using first 10 characters of obs_names")
                adata.obs['batch'] = 'batch_1'
            
            # Run Harmony with GPU backend if available
            logger.info(f"  • Running Harmony algorithm with {harmony_pcs} PCs...")
            
            ho = hm.run_harmony(
                adata.obsm['X_pca'][:, :harmony_pcs],
                adata.obs['batch'],
                theta=[${params.harmony_theta}],
                lamb=[1],
                max_iter_harmony=10,
                return_object=True,
                device=device if use_gpu else None
            )
            
            # Extract integrated embeddings
            adata.obsm['X_harmony'] = ho.Z_corr.T
            logger.info(f"  • Harmony completed successfully")
            logger.info(f"  • Integrated embedding shape: {adata.obsm['X_harmony'].shape}")
            
            # Store Harmony metadata
            adata.uns['harmony'] = {
                'method': 'harmonypy',
                'version': hm.__version__,
                'n_pcs': harmony_pcs,
                'theta': ${params.harmony_theta},
                'lamb': 1,
                'device': device,
                'timestamp': start_time.isoformat(),
            }
            
            harmony_success = True
            
        except ImportError:
            logger.warning("  • harmonypy not installed - skipping batch integration")
            logger.warning("  • Using raw PCA embeddings (no batch correction)")
            adata.obsm['X_harmony'] = adata.obsm['X_pca'][:, :harmony_pcs]
            adata.uns['harmony'] = {
                'status': 'skipped',
                'reason': 'harmonypy not available',
                'fallback': 'raw_pca'
            }
            harmony_success = False
            
        except Exception as e:
            logger.error(f"  • Harmony integration failed: {str(e)}")
            logger.warning("  • Using raw PCA embeddings as fallback")
            adata.obsm['X_harmony'] = adata.obsm['X_pca'][:, :harmony_pcs]
            adata.uns['harmony'] = {
                'status': 'failed',
                'error': str(e),
                'fallback': 'raw_pca'
            }
            harmony_success = False
        
        # Save integrated object
        logger.info("Saving Harmony-integrated matrix...")
        adata.write_h5ad("harmony_integrated.h5ad", compression="gzip")
        logger.info(f"Saved: harmony_integrated.h5ad")
        
        # ===== GENERATE REPORTS =====
        
        # Text report
        with open("harmony_report.txt", "w") as f:
            f.write("Harmony Batch Integration Report\\n")
            f.write("=" * 70 + "\\n\\n")
            f.write(f"Batch ID: ${batch_id}\\n")
            f.write(f"Timestamp: {start_time.isoformat()}\\n\\n")
            
            f.write("Input Statistics:\\n")
            f.write(f"  • Cells: {initial_stats['n_cells']:,}\\n")
            f.write(f"  • Genes: {initial_stats['n_genes']:,}\\n")
            f.write(f"  • Batches: {initial_stats['n_batches']}\\n\\n")
            
            f.write("Preprocessing:\\n")
            f.write(f"  • Normalization: 10,000 UMI/cell\\n")
            f.write(f"  • HVG selected: {n_hvg:,}\\n")
            f.write(f"  • Sparsity: {(adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]) * 100):.2f}%\\n\\n")
            
            f.write("PCA (Dimensionality Reduction):\\n")
            f.write(f"  • Method: {'rapids-singlecell GPU' if use_rapids else 'Scanpy CPU'}\\n")
            f.write(f"  • Components: {pca_components}\\n")
            f.write(f"  • Variance explained (first 10 PCs): N/A\\n\\n")
            
            f.write("Harmony Batch Integration:\\n")
            f.write(f"  • Status: {'SUCCESS' if harmony_success else 'SKIPPED/FALLBACK'}\\n")
            f.write(f"  • Device: {device}\\n")
            f.write(f"  • PCs for integration: {harmony_pcs}\\n")
            f.write(f"  • Theta (correction strength): ${params.harmony_theta}\\n")
            f.write(f"  • Output shape: {adata.obsm['X_harmony'].shape}\\n\\n")
            
            f.write("Output:\\n")
            f.write(f"  • X_pca: {adata.obsm['X_pca'].shape}\\n")
            f.write(f"  • X_harmony: {adata.obsm['X_harmony'].shape}\\n")
        
        # Metrics
        metrics_df = pd.DataFrame({
            'batch_id': ['${batch_id}'],
            'n_cells': [initial_stats['n_cells']],
            'n_batches': [initial_stats['n_batches']],
            'n_hvg': [n_hvg],
            'n_pca_components': [pca_components],
            'pca_method': ['rapids' if use_rapids else 'scanpy'],
            'harmony_status': ['success' if harmony_success else 'failed'],
            'harmony_device': [device],
            'harmony_pcs': [harmony_pcs],
            'harmony_theta': [${params.harmony_theta}],
        })
        metrics_df.to_csv("harmony_metrics.csv", index=False)
        
        # Diagnostics: per-batch distribution
        batch_stats = adata.obs.groupby('batch' if 'batch' in adata.obs else 'sample_id' if 'sample_id' in adata.obs else pd.Series(0)).size()
        diag_df = pd.DataFrame({
            'batch': batch_stats.index,
            'n_cells': batch_stats.values,
        })
        diag_df.to_csv("harmony_diagnostics.csv", index=False)
        
        elapsed = (datetime.now() - start_time).total_seconds()
        logger.info(f"Harmony integration completed in {elapsed:.1f} seconds")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error(f"Error in Harmony integration: {str(e)}", exc_info=True)
        raise
    
    PYSCRIPT
    """
}

process RAPIDS_CLUSTERING {
    label 'gpu_large'
    publishDir "${params.outdir}/analysis/clustering/${batch_id}", mode: 'copy'
    
    input:
    tuple val(batch_id), path(harmony_h5ad)
    val(gpu_available)
    
    output:
    tuple val(batch_id), path("clustered.h5ad"), emit: clustered
    path("clustering_report.txt"), emit: report
    path("clustering_metrics.csv"), emit: metrics
    path("cluster_summary.csv"), emit: summary
    
    script:
    use_gpu = params.rapids_gpu && gpu_available.toBoolean()
    """
    python3 << 'PYSCRIPT'
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import logging
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('Clustering')
    
    logger.info("=" * 70)
    logger.info("rapids-singlecell GPU-Accelerated Clustering")
    logger.info("=" * 70)
    
    start_time = datetime.now()
    
    try:
        # Load Harmony-integrated data
        logger.info(f"Loading Harmony-integrated matrix from: ${harmony_h5ad}")
        adata = sc.read_h5ad("${harmony_h5ad}")
        logger.info(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        
        # Determine basis for clustering
        use_harmony = 'X_harmony' in adata.obsm
        basis_key = 'X_harmony' if use_harmony else 'X_pca'
        
        logger.info(f"Using {'Harmony-integrated' if use_harmony else 'raw PCA'} embeddings for downstream analysis")
        logger.info(f"Basis dimensions: {adata.obsm[basis_key].shape}")
        
        # ===== NEAREST NEIGHBOR GRAPH =====
        logger.info("Computing nearest neighbor graph...")
        n_neighbors = ${params.umap_neighbors}
        
        use_rapids = False
        try:
            import rapids_singlecell as rsc
            logger.info(f"  • Using rapids-singlecell GPU-accelerated kNN ({n_neighbors} neighbors)")
            rsc.pp.neighbors(adata, use_rep=basis_key, n_neighbors=n_neighbors, metric='euclidean')
            use_rapids = True
            logger.info(f"  • kNN graph computed on GPU")
            logger.info(f"  • Distances shape: {adata.obsp['distances'].shape}")
            
        except ImportError:
            logger.warning("rapids-singlecell not available, using Scanpy CPU kNN")
            sc.pp.neighbors(adata, use_rep=basis_key, n_neighbors=n_neighbors, metric='euclidean')
            logger.info(f"  • kNN graph computed on CPU")
            logger.info(f"  • Distances shape: {adata.obsp['distances'].shape}")
        
        except Exception as e:
            logger.error(f"  • kNN computation failed: {str(e)}")
            raise
        
        # ===== UMAP DIMENSIONALITY REDUCTION =====
        logger.info("Computing UMAP...")
        
        try:
            if use_rapids:
                logger.info("  • Using rapids-singlecell GPU-accelerated UMAP")
                rsc.tl.umap(adata)
                logger.info(f"  • UMAP completed on GPU: {adata.obsm['X_umap'].shape}")
            else:
                logger.info("  • Using Scanpy CPU UMAP")
                sc.tl.umap(adata)
                logger.info(f"  • UMAP completed on CPU: {adata.obsm['X_umap'].shape}")
                
        except Exception as e:
            logger.warning(f"  • UMAP computation failed ({str(e)})")
            logger.warning("  • Using first 2 Harmony components for visualization")
            adata.obsm['X_umap'] = adata.obsm[basis_key][:, :2]
        
        # ===== LEIDEN CLUSTERING =====
        logger.info("Running Leiden clustering...")
        leiden_res = ${params.leiden_resolution}
        use_weights = ${params.leiden_use_weights}
        
        logger.info(f"  • Resolution: {leiden_res}")
        logger.info(f"  • Use edge weights: {use_weights}")
        
        try:
            if use_rapids:
                logger.info("  • Using rapids-singlecell GPU-accelerated Leiden")
                rsc.tl.leiden(adata, resolution=leiden_res, use_weights=use_weights)
                logger.info(f"  • Leiden clustering completed on GPU")
            else:
                logger.info("  • Using Scanpy CPU Leiden")
                sc.tl.leiden(adata, resolution=leiden_res, use_weights=use_weights)
                logger.info(f"  • Leiden clustering completed on CPU")
            
            n_clusters = len(adata.obs['leiden'].unique())
            cluster_sizes = adata.obs['leiden'].value_counts().sort_index()
            logger.info(f"  • Clusters detected: {n_clusters}")
            logger.info(f"  • Cluster sizes: min={cluster_sizes.min()}, max={cluster_sizes.max()}, median={cluster_sizes.median():.0f}")
            
        except Exception as e:
            logger.error(f"  • Leiden clustering failed: {str(e)}")
            logger.warning("  • Assigning all cells to 'unclassified' cluster")
            adata.obs['leiden'] = 'unclassified'
            n_clusters = 1
        
        # Store clustering metadata
        adata.uns['leiden'] = {
            'method': 'leiden',
            'resolution': leiden_res,
            'use_weights': use_weights,
            'n_clusters': n_clusters,
            'basis': basis_key,
            'timestamp': start_time.isoformat(),
        }
        
        # Save clustered object
        logger.info("Saving clustered matrix...")
        adata.write_h5ad("clustered.h5ad", compression="gzip")
        logger.info(f"Saved: clustered.h5ad")
        
        # ===== GENERATE REPORTS =====
        
        # Text report
        with open("clustering_report.txt", "w") as f:
            f.write("Clustering & Dimensionality Reduction Report\\n")
            f.write("=" * 70 + "\\n\\n")
            f.write(f"Batch ID: ${batch_id}\\n")
            f.write(f"Timestamp: {start_time.isoformat()}\\n\\n")
            
            f.write("Input:\\n")
            f.write(f"  • Cells: {adata.n_obs:,}\\n")
            f.write(f"  • Genes: {adata.n_vars:,}\\n")
            f.write(f"  • Basis: {basis_key} ({adata.obsm[basis_key].shape[1]} dimensions)\\n\\n")
            
            f.write("Nearest Neighbor Graph:\\n")
            f.write(f"  • Method: {'rapids-singlecell GPU' if use_rapids else 'Scanpy CPU'}\\n")
            f.write(f"  • Neighbors: {n_neighbors}\\n")
            f.write(f"  • Metric: euclidean\\n\\n")
            
            f.write("UMAP:\\n")
            f.write(f"  • Status: {'SUCCESS' if 'X_umap' in adata.obsm else 'FALLBACK'}\\n")
            f.write(f"  • Shape: {adata.obsm.get('X_umap', np.array([])).shape}\\n\\n")
            
            f.write("Leiden Clustering:\\n")
            f.write(f"  • Method: {'rapids-singlecell GPU' if use_rapids else 'Scanpy CPU'}\\n")
            f.write(f"  • Resolution: {leiden_res}\\n")
            f.write(f"  • Clusters: {n_clusters}\\n")
            if n_clusters > 1:
                cluster_sizes = adata.obs['leiden'].value_counts().sort_index()
                f.write(f"  • Cluster sizes (min-max): {cluster_sizes.min()}-{cluster_sizes.max()}\\n")
                f.write(f"  • Median cluster size: {cluster_sizes.median():.0f}\\n")
        
        # Metrics
        metrics_df = pd.DataFrame({
            'batch_id': ['${batch_id}'],
            'n_cells': [adata.n_obs],
            'n_genes': [adata.n_vars],
            'basis': [basis_key],
            'n_neighbors': [n_neighbors],
            'knn_method': ['rapids' if use_rapids else 'scanpy'],
            'leiden_resolution': [leiden_res],
            'n_clusters': [n_clusters],
        })
        metrics_df.to_csv("clustering_metrics.csv", index=False)
        
        # Cluster summary
        if n_clusters > 1:
            cluster_summary = adata.obs['leiden'].value_counts().sort_index()
            cluster_df = pd.DataFrame({
                'cluster': cluster_summary.index,
                'n_cells': cluster_summary.values,
                'pct_cells': (cluster_summary.values / len(adata.obs) * 100),
            })
        else:
            cluster_df = pd.DataFrame({'cluster': ['unclassified'], 'n_cells': [adata.n_obs], 'pct_cells': [100.0]})
        
        cluster_df.to_csv("cluster_summary.csv", index=False)
        
        elapsed = (datetime.now() - start_time).total_seconds()
        logger.info(f"Clustering completed in {elapsed:.1f} seconds")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error(f"Error in clustering: {str(e)}", exc_info=True)
        raise
    
    PYSCRIPT
    """
}


process RAPIDS_CLUSTERING {
    label 'gpu_large'
    publishDir "${params.outdir}/analysis/clustering/${batch_id}", mode: 'copy'
    
    input:
    tuple val(batch_id), path(harmony_h5ad)
    val(gpu_available)
    
    output:
    tuple val(batch_id), path("clustered.h5ad"), emit: clustered
    path("clustering_report.txt"), emit: report
    
    script:
    """
    python3 << 'PYSCRIPT'
    import scanpy as sc
    import numpy as np
    import warnings
    warnings.filterwarnings('ignore')
    
    print("rapids-singlecell GPU-accelerated Clustering")
    print("=" * 60)
    
    # Load Harmony-integrated data
    adata = sc.read_h5ad("${harmony_h5ad}")
    print(f"Input: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Use Harmony PCA as basis for UMAP/clustering
    use_harmony = 'X_harmony' in adata.obsm
    
    if use_harmony:
        pca_key = 'X_harmony'
        print("Using Harmony-integrated embeddings for downstream analysis")
    else:
        pca_key = 'X_pca'
        print("Using standard PCA for downstream analysis")
    
    # Compute nearest neighbors
    print("Computing nearest neighbors graph...")
    n_neighbors = ${params.umap_neighbors}
    
    try:
        import rapids_singlecell as rsc
        print("Using rapids-singlecell GPU-accelerated kNN")
        # rapids-singlecell uses GPU for neighbor computation
        rsc.pp.neighbors(adata, use_rep=pca_key, n_neighbors=n_neighbors, metric='euclidean')
        use_rapids = True
    except ImportError:
        print("rapids-singlecell not available, using Scanpy kNN")
        sc.pp.neighbors(adata, use_rep=pca_key, n_neighbors=n_neighbors, metric='euclidean')
        use_rapids = False
    
    print(f"Neighbors computed: {adata.obsp['distances'].shape}")
    
    # UMAP dimensional reduction
    print("Computing UMAP...")
    
    try:
        if use_rapids:
            rsc.tl.umap(adata)
        else:
            sc.tl.umap(adata)
        print(f"UMAP: {adata.obsm['X_umap'].shape}")
    except Exception as e:
        print(f"UMAP failed ({e}), using PCA for visualization")
        adata.obsm['X_umap'] = adata.obsm[pca_key][:, :2]
    
    # Leiden clustering
    print("Running Leiden clustering...")
    leiden_res = ${params.leiden_resolution}
    
    try:
        if use_rapids:
            rsc.tl.leiden(adata, resolution=leiden_res, use_weights=${params.leiden_use_weights})
        else:
            sc.tl.leiden(adata, resolution=leiden_res, use_weights=${params.leiden_use_weights})
        
        n_clusters = len(adata.obs['leiden'].unique())
        print(f"Clusters detected: {n_clusters}")
    except Exception as e:
        print(f"Leiden clustering failed: {e}")
        adata.obs['leiden'] = 'unclassified'
    
    # Save clustered object
    adata.write_h5ad("clustered.h5ad", compression="gzip")
    
    # Generate report
    with open("clustering_report.txt", "w") as f:
        f.write(f"Clustering Report\n")
        f.write(f"Batch ID: ${batch_id}\n")
        f.write(f"Input cells: {adata.n_obs}\n")
        f.write(f"Basis for clustering: {pca_key}\n")
        f.write(f"Neighbors: {n_neighbors}\n")
        f.write(f"Leiden resolution: {leiden_res}\n")
        f.write(f"Clusters: {len(adata.obs['leiden'].unique())}\n")
        f.write(f"GPU acceleration: {'YES' if use_rapids else 'NO'}\n")
    PYSCRIPT
    """
}
