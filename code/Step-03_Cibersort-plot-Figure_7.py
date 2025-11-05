import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.backends.backend_pdf import PdfPages


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    os.makedirs(os.path.join(path_save, "Cibersort", "Figures"), exist_ok=True)
    
    cibersort_lm22 = pd.read_table(os.path.join(path_save, "Cibersort", "CIBERSORT_Results_LM22.txt"), 
                                   sep="\t", header=0, index_col=0)
    cibersort_lm14 = pd.read_table(os.path.join(path_save, "Cibersort", "CIBERSORT_Results_LM14.txt"), 
                                   sep="\t", header=0, index_col=0)
    
    cibersort_results = pd.concat([cibersort_lm22, cibersort_lm14], axis=1)
    if 'NK.cells.activated' in cibersort_results.columns and 'NK.cells.resting' in cibersort_results.columns:
        cibersort_results['NK'] = cibersort_results['NK.cells.activated'] + cibersort_results['NK.cells.resting']
    
    sample_group = pd.read_table(os.path.join(path_save, "Data_ALL_group.txt"), 
                                 sep="\t", header=0, index_col=0).T
    sample_group = sample_group.loc[cibersort_results.index, :]
    cibersort_results['group'] = sample_group['group']
    
    comparisons = [
        ["Cluster1", "GTEx"],
        ["Cluster2", "GTEx"],
        ["Cluster3", "GTEx"],
        ["Cluster4", "GTEx"]
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
    
    for i, cell in enumerate(cells):
        if cell not in cibersort_results.columns:
            print(f"Warning: {cell} not found in CIBERSORT results")
            continue
        
        data_plot = cibersort_results[[cell, 'group']].copy()
        data_plot = data_plot.melt(id_vars=['group'], value_vars=[cell], 
                                   var_name='CellType', value_name='Composition')
        
        fig, ax = plt.subplots(figsize=(6, 4))
        
        sns.violinplot(data=data_plot, x='group', y='Composition', hue='group', 
                      palette=cluster_colors, ax=ax, inner=None)
        sns.boxplot(data=data_plot, x='group', y='Composition', width=0.12, ax=ax,
                   color='white', boxprops=dict(alpha=0.7))
        
        groups = data_plot['group'].unique()
        group_positions = {g: i for i, g in enumerate(sorted(groups))}
        
        for comp in comparisons:
            if comp[0] in groups and comp[1] in groups:
                g1_data = data_plot[data_plot['group'] == comp[0]]['Composition'].values
                g2_data = data_plot[data_plot['group'] == comp[1]]['Composition'].values
                
                if len(g1_data) > 0 and len(g2_data) > 0:
                    stat, pval = mannwhitneyu(g1_data, g2_data)
                    y_pos = height_list[i][comparisons.index(comp)]
                    
                    if pval < 0.001:
                        sig = '***'
                    elif pval < 0.01:
                        sig = '**'
                    elif pval < 0.05:
                        sig = '*'
                    else:
                        sig = 'ns'
                    
                    x1 = group_positions[comp[0]]
                    x2 = group_positions[comp[1]]
                    ax.plot([x1, x2], [y_pos, y_pos], 'k-', linewidth=0.5)
                    ax.text((x1 + x2) / 2, y_pos, sig, ha='center', va='bottom', 
                           fontsize=10, color='#6699cc')
        
        ax.set_ylim(0, ylims[i])
        ax.set_xlabel('')
        ax.set_ylabel('Percentage')
        ax.set_title(cell)
        ax.legend().remove()
        ax.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        fig_list.append(fig)
        plt.close(fig)
    
    if len(fig_list) > 0:
        fig_combined, axes = plt.subplots(1, 3, figsize=(15, 4))
        for i, fig in enumerate(fig_list):
            if i < len(axes):
                axes[i].remove()
                fig_combined.add_axes(axes[i].get_position())
                # Copy axes content would be complex, so we'll save individual plots
                fig.savefig(os.path.join(path_save, "Cibersort", "Figures", f"cibersort_{cells[i]}.png"))
        
        pdf_path = os.path.join(path_save, "Cibersort", "Figures", "Figure 7.pdf")
        with PdfPages(pdf_path) as pdf:
            for fig in fig_list:
                pdf.savefig(fig)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

