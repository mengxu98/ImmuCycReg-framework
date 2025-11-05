import os
import sys
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_deseq2_pydeseq2(tcga_cluster_df: pd.DataFrame, gtex_df: pd.DataFrame) -> pd.DataFrame:
    """Run DESeq2 analysis using PyDESeq2"""
    tcga_gtex = pd.concat([tcga_cluster_df, gtex_df], axis=1)
    group = ["C1"] * tcga_cluster_df.shape[1] + ["C2"] * gtex_df.shape[1]
    
    clinical_df = pd.DataFrame(
        {"condition": group},
        index=tcga_gtex.columns
    )
    
    counts_df = tcga_gtex.round().astype(int)
    counts_df = counts_df.fillna(0)
    
    try:
        dds = DeseqDataSet(
            counts=counts_df,
            clinical=clinical_df,
            design_factors="condition",
            refit_cooks=True,
            min_mu=0.5
        )
        dds.deseq2()
        
        stat_res = DeseqStats(dds, contrast=["condition", "C1", "C2"])
        stat_res.summary()
        res_df = stat_res.results_df
        
        res_df.index.name = 'gene'
        return res_df
    except Exception as e:
        print(f"Error in PyDESeq2: {e}")
        import traceback
        traceback.print_exc()
        raise


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, 'DESeq2'), exist_ok=True)
    os.makedirs(os.path.join(path_save, 'Figures'), exist_ok=True)
    
    genes_249_file = os.path.join(path_read, "Genes_249.txt")
    if os.path.exists(genes_249_file):
        genes_249 = pd.read_table(genes_249_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Warning: Genes_249.txt not found, using all genes")
        genes_249 = None
    
    samples_cluster = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"), header=None)
    
    tcga = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    gtex = pd.read_csv(os.path.join(path_save, "GTEx-LUAD.csv"), index_col=0)
    
    if genes_249 is not None:
        tcga = tcga.reindex(index=genes_249)
        gtex = gtex.reindex(index=genes_249)
    
    volcano_plot_list = []
    
    for j in range(1, 5):
        print(f"Cluster {j} DESeq2 start......")
        cls_samples = samples_cluster.loc[samples_cluster.iloc[:, 1] == j, 0].astype(str).tolist()
        cls_samples = [s for s in cls_samples if s in tcga.columns]
        if len(cls_samples) == 0:
            print(f"Warning: No samples found for cluster {j}")
            continue
        
        tcga_cluster = tcga.loc[:, cls_samples].copy().dropna()
        gtex_cluster = gtex.copy().dropna()
        
        common_genes = list(set(tcga_cluster.index) & set(gtex_cluster.index))
        if len(common_genes) == 0:
            print(f"Warning: No common genes for cluster {j}")
            continue
        
        tcga_cluster = tcga_cluster.loc[common_genes, :]
        gtex_cluster = gtex_cluster.loc[common_genes, :]
        
        try:
            res_df = run_deseq2_pydeseq2(tcga_cluster, gtex_cluster)
        except Exception as e:
            print(f"Error processing cluster {j}: {e}")
            continue
        
        res_df = res_df.sort_values(by='padj', na_last=True)
        
        padj_threshold = 0.05
        fold_change = 1.0
        res_df['reg'] = 'not DE'
        res_df.loc[(res_df['log2FoldChange'] >= 0.5), 'reg'] = 'min-up-regulated'
        res_df.loc[(res_df['log2FoldChange'] <= -0.5), 'reg'] = 'min-down-regulated'
        res_df.loc[(res_df['log2FoldChange'] >= fold_change) & (res_df['padj'] <= padj_threshold), 'reg'] = 'up-regulated'
        res_df.loc[(res_df['log2FoldChange'] <= -fold_change) & (res_df['padj'] <= padj_threshold), 'reg'] = 'down-regulated'
        
        if genes_249 is not None:
            res_df = res_df.reindex(index=genes_249)
        
        res_df = res_df.reset_index()
        res_df = res_df.rename(columns={'index': 'gene'})
        
        col_order = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'reg']
        res_df = res_df[[col for col in col_order if col in res_df.columns]]
        
        res_df.to_csv(os.path.join(path_save, f"DESeq2_results_cluster{j}.csv"), index=False)
        
        data_plot = res_df.copy()
        data_plot['change'] = 'no sig'
        data_plot.loc[(data_plot['pvalue'] < padj_threshold) & (data_plot['log2FoldChange'] >= fold_change), 'change'] = 'Up (log2FC > 1, P-value < 0.05)'
        data_plot.loc[(data_plot['pvalue'] < padj_threshold) & (data_plot['log2FoldChange'] <= -fold_change), 'change'] = 'Down (log2FC < -1, P-value < 0.05)'
        data_plot['label'] = data_plot.apply(lambda x: x['gene'] if (x['pvalue'] < padj_threshold and abs(x['log2FoldChange']) >= 1) else '', axis=1)
        
        fig, ax = plt.subplots(figsize=(6, 6))
        colors = {'Up (log2FC > 1, P-value < 0.05)': '#990033', 'no sig': '#d2dae2', 'Down (log2FC < -1, P-value < 0.05)': '#006699'}
        for change_type in ['Up (log2FC > 1, P-value < 0.05)', 'no sig', 'Down (log2FC < -1, P-value < 0.05)']:
            subset = data_plot[data_plot['change'] == change_type]
            if len(subset) > 0:
                ax.scatter(subset['log2FoldChange'], -np.log10(subset['pvalue'] + 1e-300), 
                          c=colors[change_type], alpha=0.5, s=30, label=change_type)
        
        labeled = data_plot[data_plot['label'] != '']
        if len(labeled) > 0:
            for _, row in labeled.head(20).iterrows():
                ax.annotate(row['label'], (row['log2FoldChange'], -np.log10(row['pvalue'] + 1e-300)),
                           fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        
        ax.axvline(x=-fold_change, linestyle='--', color='black', linewidth=0.8)
        ax.axvline(x=fold_change, linestyle='--', color='black', linewidth=0.8)
        ax.axhline(y=-np.log10(padj_threshold), linestyle='--', color='black', linewidth=0.8)
        ax.set_xlabel('Log2FoldChange')
        ax.set_ylabel('-Log10(P-value)')
        ax.set_title(f'Cluster {j}')
        ax.legend(loc='right')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        volcano_plot_list.append(fig)
        plt.close(fig)
        
        print(f"Cluster {j} DESeq2 done......")
    
    if len(volcano_plot_list) > 0:
        pdf_path = os.path.join(path_save, "Figures", "Figure 3.pdf")
        with PdfPages(pdf_path) as pdf:
            for fig in volcano_plot_list:
                pdf.savefig(fig)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
