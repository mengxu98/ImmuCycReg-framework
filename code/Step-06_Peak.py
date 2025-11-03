import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def check_dir(dirpath):
    os.makedirs(dirpath, exist_ok=True)
    return dirpath


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    
    tcga_luad = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    peak_luad = pd.read_csv(os.path.join(path_save, "Peak-LUAD.csv"), index_col=0)
    genes_adj_peak = pd.read_csv(os.path.join(path_save, "ATAC-LUAD.csv"), index_col=0)
    geneinfo_df = pd.read_csv(os.path.join(path_save, "Geneinfo_df.csv"), index_col=0)
    
    peak_name = genes_adj_peak.index.tolist()
    peak_gene_list = [str(genes_adj_peak.iloc[i, 0]).split(",") if pd.notna(genes_adj_peak.iloc[i, 0]) else [] 
                       for i in range(len(genes_adj_peak))]
    
    sample_atac = list(set(tcga_luad.columns) & set(peak_luad.columns))
    
    genes_17_file = os.path.join(path_read, "Genes_17.txt")
    if os.path.exists(genes_17_file):
        gene_list = pd.read_table(genes_17_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Warning: Genes_17.txt not found")
        return 1
    
    check_dir(os.path.join(path_save, "BED"))
    check_dir(os.path.join(path_save, "corrgram"))
    check_dir(os.path.join(path_save, "Figure"))
    
    peak_seq_pos_all = []
    volcano_plot_list = []
    
    for j in range(len(gene_list)):
        target_gene = gene_list[j]
        print(f"Running for {target_gene}......")
        
        peaks_around_gene = []
        for i in range(len(peak_name)):
            if target_gene in peak_gene_list[i]:
                peaks_around_gene.append(peak_name[i])
        
        if len(peaks_around_gene) == 0:
            continue
        
        check_dir(os.path.join(path_save, "corrgram", target_gene))
        
        mrna = tcga_luad.loc[target_gene, sample_atac].to_frame().T
        atac = peak_luad.loc[peaks_around_gene, sample_atac]
        
        peak_gene_corr = []
        peak_gene_pval = []
        for i in range(len(peaks_around_gene)):
            peak_values = atac.iloc[i, :].astype(float).to_numpy()
            mrna_values = mrna.iloc[0, :].astype(float).to_numpy()
            mask = np.isfinite(peak_values) & np.isfinite(mrna_values)
            if mask.sum() < 3:
                peak_gene_corr.append(np.nan)
                peak_gene_pval.append(1.0)
            else:
                r, p = spearmanr(peak_values[mask], mrna_values[mask])
                peak_gene_corr.append(r)
                peak_gene_pval.append(p)
        
        deg_data = pd.DataFrame({
            'Symbol': peaks_around_gene,
            'corr': peak_gene_corr,
            'pval': peak_gene_pval
        })
        
        deg_data['logP'] = -np.log10(deg_data['pval'] + 1e-300)
        deg_data['Group'] = 'Not-significant'
        deg_data.loc[(deg_data['pval'] < 0.05) & (deg_data['corr'] > 0.1), 'Group'] = 'Up-regulated'
        deg_data.loc[(deg_data['pval'] < 0.05) & (deg_data['corr'] < -0.1), 'Group'] = 'Down-regulated'
        
        deg_data['Label'] = ''
        deg_data = deg_data.sort_values('pval')
        
        up_peaks = deg_data.loc[deg_data['Group'] == 'Up-regulated', 'Symbol'].head(20).tolist()
        down_peaks = deg_data.loc[deg_data['Group'] == 'Down-regulated', 'Symbol'].head(20).tolist()
        not_sig_peaks = deg_data.loc[deg_data['Group'] == 'Not-significant', 'Symbol'].head(20).tolist()
        label_peaks = up_peaks + down_peaks + not_sig_peaks
        
        deg_data.loc[deg_data['Symbol'].isin(label_peaks), 'Label'] = deg_data.loc[deg_data['Symbol'].isin(label_peaks), 'Symbol']
        
        deg_data.to_csv(os.path.join(path_save, "corrgram", target_gene, f"{target_gene}_peakGeneCor.csv"), 
                       index=False)
        
        high_cor_peak = up_peaks[0] if len(up_peaks) > 0 else None
        down_cor_peak = down_peaks[0] if len(down_peaks) > 0 else None
        
        pos_col_names = ['chrom', 'chromStart', 'chromEnd', 'strand']
        temp_res = geneinfo_df[geneinfo_df['gene_name'] == target_gene].copy()
        temp_res = temp_res[temp_res['type'] == 'transcript'].copy()
        
        nearest_peak = None
        if len(temp_res) > 0:
            temp_res = temp_res.sort_values('start')
            most_start = temp_res.iloc[0]['start']
            most_end = temp_res.iloc[0]['end']
            most_mid = (most_start + most_end) / 2
            
            peaks_pos_temp = genes_adj_peak.loc[peaks_around_gene, pos_col_names].copy()
            peaks_pos_temp['mid'] = (peaks_pos_temp['chromStart'] + peaks_pos_temp['chromEnd']) / 2
            peaks_pos_temp['mid_distance'] = peaks_pos_temp['mid'] - most_mid
            peaks_pos_temp['start_distance'] = peaks_pos_temp['chromEnd'] - most_start
            distance_thres = -8000
            
            valid_peaks = peaks_pos_temp[
                (peaks_pos_temp['start_distance'] < 0) &
                (peaks_pos_temp['start_distance'] > distance_thres) &
                (peaks_pos_temp['chromEnd'] < most_start)
            ]
            
            if len(valid_peaks) > 0:
                nearest_peak = valid_peaks.index[-1]
        
        order_peaks = peaks_around_gene[np.argsort(peak_luad.loc[peaks_around_gene, sample_atac[0]].values)[::-1]]
        highest_score_peak = order_peaks[0] if len(order_peaks) > 0 else None
        
        candidate_peak = None
        if nearest_peak and nearest_peak in [high_cor_peak, down_cor_peak]:
            candidate_peak = nearest_peak
        elif highest_score_peak and highest_score_peak in [high_cor_peak, down_cor_peak]:
            candidate_peak = highest_score_peak
        elif high_cor_peak:
            candidate_peak = high_cor_peak
        elif down_cor_peak:
            candidate_peak = down_cor_peak
        
        if candidate_peak:
            peak_seq_pos = genes_adj_peak.loc[[candidate_peak], pos_col_names].copy()
            peak_seq_pos['chromStart'] = peak_seq_pos['chromStart'] - 1
            peak_seq_pos['gene'] = target_gene
            peak_seq_pos['peak'] = candidate_peak
            peak_seq_pos_all.append(peak_seq_pos)
            
            palette = {
                'Down-regulated': '#006699',
                'Not-significant': 'gray',
                'Up-regulated': '#EE3A8C'
            }
            
            fig, ax = plt.subplots(figsize=(6, 5))
            for group in ['Up-regulated', 'Down-regulated', 'Not-significant']:
                data = deg_data[deg_data['Group'] == group]
                ax.scatter(data['corr'], data['logP'], c=palette[group], label=group, s=30, alpha=0.6)
            
            labeled_data = deg_data[deg_data['Label'] != '']
            for _, row in labeled_data.iterrows():
                ax.annotate(row['Label'], (row['corr'], row['logP']), fontsize=8, 
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
            
            ax.axhline(y=1.5, linestyle='--', color='black', linewidth=0.8)
            ax.axvline(x=-0.3, linestyle='--', color='black', linewidth=0.8)
            ax.axvline(x=0.3, linestyle='--', color='black', linewidth=0.8)
            ax.set_xlabel(f'{target_gene} correlation')
            ax.set_ylabel('-Log10(P-value)')
            ax.set_title(target_gene)
            ax.legend(loc='bottom')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            
            volcano_plot_list.append(fig)
            plt.close(fig)
    
    if len(volcano_plot_list) > 0:
        pdf_path = os.path.join(path_save, "Figure", "Supplementary Figure 1-peaks.pdf")
        with PdfPages(pdf_path) as pdf:
            for fig in volcano_plot_list:
                pdf.savefig(fig)
    
    if len(peak_seq_pos_all) > 0:
        peak_seq_pos_all_df = pd.concat(peak_seq_pos_all, axis=0)
        
        bed_df = peak_seq_pos_all_df[['chrom', 'chromStart', 'chromEnd']].copy()
        bed_df.to_csv(os.path.join(path_save, "all_BED.bed"), sep="\t", header=False, index=False)
        
        peaks_csv = peak_seq_pos_all_df[['gene', 'peak', 'chrom', 'chromStart', 'chromEnd', 'strand']].copy()
        peaks_csv.to_csv(os.path.join(path_save, "all_peaks.csv"), index=False)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

