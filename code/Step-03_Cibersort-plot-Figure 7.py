import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.backends.backend_pdf import PdfPages


def check_dir(dirpath):
    os.makedirs(dirpath, exist_ok=True)
    return dirpath


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    cibersort_path = os.path.join(path_save, "Cibersort")
    os.makedirs(os.path.join(cibersort_path, "Figures"), exist_ok=True)
    
    cibersort_lm22 = pd.read_csv(os.path.join(cibersort_path, "CIBERSORT_Results_LM22.txt"), 
                                 sep="\t", index_col=0)
    cibersort_lm14 = pd.read_csv(os.path.join(cibersort_path, "CIBERSORT_Results_LM14.txt"), 
                                 sep="\t", index_col=0)
    
    cibersort_results = pd.concat([cibersort_lm22, cibersort_lm14], axis=1)
    if 'NK.cells.activated' in cibersort_results.columns and 'NK.cells.resting' in cibersort_results.columns:
        cibersort_results['NK'] = (cibersort_results['NK.cells.activated'] + 
                                   cibersort_results['NK.cells.resting'])
    
    sample_group = pd.read_csv(os.path.join(path_save, "Data_ALL_group.txt"), 
                               sep="\t", index_col=0)
    sample_group = sample_group.T
    sample_group = sample_group.loc[cibersort_results.index, :]
    cibersort_results['group'] = sample_group['group'].values
    
    comparisons = [
        ("Cluster1", "GTEx"),
        ("Cluster2", "GTEx"),
        ("Cluster3", "GTEx"),
        ("Cluster4", "GTEx")
    ]
    
    cluster_colors = ["#2874c5", "#008a00", "#c6524a", "#eabf00", "#696969"]
    
    cells = ["Th.cell", "Macrophages.M2", "NK.cells.activated"]
    
    height_list = [
        [0.7, 0.65, 0.6, 0.55],
        [0.7, 0.65, 0.6, 0.55],
        [0.23, 0.21, 0.19, 0.17]
    ]
    ylims = [0.75, 0.75, 0.25]
    
    fig_list = []
    
    for i, cell_type in enumerate(cells):
        if cell_type not in cibersort_results.columns:
            continue
        
        data = cibersort_results[[cell_type, 'group']].copy()
        data = data.rename(columns={cell_type: 'Composition'})
        
        fig, ax = plt.subplots(figsize=(8, 4))
        
        groups = ['Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'GTEx']
        group_data = [data[data['group'] == g]['Composition'].values for g in groups]
        
        positions = np.arange(len(groups))
        
        violin_parts = ax.violinplot(group_data, positions=positions, widths=0.6, 
                                    showmeans=False, showmedians=False, showextrema=False)
        
        for j, pc in enumerate(violin_parts['bodies']):
            pc.set_facecolor(cluster_colors[j])
            pc.set_alpha(0.6)
        
        bp = ax.boxplot(group_data, positions=positions, widths=0.12, 
                       patch_artist=True, showfliers=False)
        for patch in bp['boxes']:
            patch.set_facecolor('white')
            patch.set_edgecolor('black')
        
        for j, group in enumerate(groups):
            group_values = data[data['group'] == group]['Composition'].values
            if len(group_values) > 0:
                mean_val = np.mean(group_values)
                std_val = np.std(group_values)
                ax.errorbar(j, mean_val, yerr=std_val, fmt='o', color='black', 
                           capsize=3, capthick=1.5, markersize=5)
        
        for comp_idx, (group1, group2) in enumerate(comparisons):
            if group1 in groups and group2 in groups:
                data1 = data[data['group'] == group1]['Composition'].values
                data2 = data[data['group'] == group2]['Composition'].values
                
                if len(data1) > 0 and len(data2) > 0:
                    try:
                        stat, pval = mannwhitneyu(data1, data2, alternative='two-sided')
                        
                        if pval < 0.001:
                            sig = '***'
                        elif pval < 0.01:
                            sig = '**'
                        elif pval < 0.05:
                            sig = '*'
                        else:
                            sig = 'ns'
                        
                        idx1 = groups.index(group1)
                        idx2 = groups.index(group2)
                        y_pos = height_list[i][comp_idx] * ylims[i]
                        
                        ax.plot([idx1, idx2], [y_pos, y_pos], 'k-', linewidth=0.3)
                        ax.text((idx1 + idx2) / 2, y_pos, sig, ha='center', va='bottom',
                               fontsize=10, color='#6699cc', fontweight='bold')
                    except:
                        pass
        
        ax.set_xticks(positions)
        ax.set_xticklabels(groups, rotation=45, ha='right', fontsize=10)
        ax.set_ylabel('Percentage', fontsize=12)
        ax.set_xlabel('')
        ax.set_ylim(0, ylims[i])
        ax.set_title(cell_type)
        ax.grid(True, alpha=0.3, axis='y')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        fig_list.append(fig)
    
    if len(fig_list) > 0:
        fig_combined, axes = plt.subplots(1, len(fig_list), figsize=(9 * len(fig_list), 4))
        if len(fig_list) == 1:
            axes = [axes]
        
        for i, fig in enumerate(fig_list):
            fig.savefig(os.path.join(cibersort_path, "Figures", f"Figure_7_{cells[i]}.png"), 
                       dpi=300, bbox_inches='tight')
        
        pdf_path = os.path.join(cibersort_path, "Figures", "Figure 7.pdf")
        with PdfPages(pdf_path) as pdf:
            for fig in fig_list:
                pdf.savefig(fig, bbox_inches='tight')
        
        for fig in fig_list:
            plt.close(fig)
        plt.close(fig_combined)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
