import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def format_regulation(input_data, target_gene, regulation):
    """Format regulation results from list of samples to DataFrame"""
    if regulation not in [1, "1", "positive", -1, "-1", "negative"]:
        raise ValueError("regulation must be 1/'positive' or -1/'negative'")
    
    tfs = []
    types = []
    target_genes = []
    samples = []
    
    if input_data and len(input_data) > 0:
        for sample_name, tf_list in input_data.items():
            if tf_list and len(tf_list) > 0:
                for tf in tf_list:
                    tfs.append(tf)
                    types.append(str(regulation))
                    target_genes.append(target_gene)
                    samples.append(sample_name)
        
        if len(tfs) > 0:
            return pd.DataFrame({
                'TF': tfs,
                'Type': types,
                'Gene': target_genes,
                'Sample': samples
            })
    
    return pd.DataFrame()


def frame_regulatory_table(regulatory_table, gene, databases_list):
    """Add database matching columns to regulatory table"""
    result_table = regulatory_table.copy()
    
    for database_name, database_df in databases_list.items():
        weights = []
        for tf in regulatory_table.iloc[:, 0]:
            matches = database_df[
                (database_df.iloc[:, 0] == gene) & 
                (database_df.iloc[:, 1] == tf)
            ]
            if len(matches) > 0:
                weight = matches.iloc[0, 2] if len(matches.columns) > 2 else 1
                weights.append(float(weight))
            else:
                weights.append(0)
        
        result_table[database_name] = weights
    
    return result_table


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    
    os.makedirs(os.path.join(path_save, "Regulatory_table"), exist_ok=True)
    os.makedirs(os.path.join(path_save, "TFs_list_selected"), exist_ok=True)
    os.makedirs(os.path.join(path_save, "T-test"), exist_ok=True)
    os.makedirs(os.path.join(path_save, "Figure"), exist_ok=True)
    
    string_data_file = os.path.join(path_read, "Databases", "STRING_9606.protein.v11.5.txt.gz")
    if os.path.exists(string_data_file):
        string_data = pd.read_csv(string_data_file, sep="\t", header=None)
    else:
        string_data = pd.DataFrame()
    
    databases = ["ASPAR", "ENCODE", "CHEA", "MotifMap", "TRANSFAC_Curated", "TRANSFAC_Predicted"]
    databases_list = {}
    for database in databases:
        db_file = os.path.join(path_read, "Databases", f"{database}.txt.gz")
        if os.path.exists(db_file):
            databases_list[database] = pd.read_csv(db_file, sep="\t", header=0)
        else:
            databases_list[database] = pd.DataFrame()
    
    tcga_luad = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    gtex_luad = pd.read_csv(os.path.join(path_save, "GTEx-LUAD.csv"), index_col=0)
    tcga_cnv = pd.read_csv(os.path.join(path_save, "CNV-LUAD.csv"), index_col=0)
    peak_luad = pd.read_csv(os.path.join(path_save, "Peak-LUAD.csv"), index_col=0)
    genes_adj_peak = pd.read_csv(os.path.join(path_save, "ATAC-LUAD.csv"), index_col=0)
    geneinfo_df = pd.read_csv(os.path.join(path_save, "Geneinfo_df.csv"), index_col=0)
    
    candidate_peaks = pd.read_csv(os.path.join(path_save, "all_peaks.csv"))
    
    negative_genes_file = os.path.join(path_read, "Genes_immune_cycle.txt")
    if os.path.exists(negative_genes_file):
        negative_genes_df = pd.read_table(negative_genes_file, header=0, sep="\t")
        negative_genes = negative_genes_df[negative_genes_df['Direction'] == 'negative'].iloc[:, 0].tolist()
    else:
        negative_genes = []
    
    genes_17_file = os.path.join(path_read, "Genes_17.txt")
    if os.path.exists(genes_17_file):
        target_genes_list = pd.read_table(genes_17_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Error: Genes_17.txt not found")
        return 1
    
    cnv_test = False
    
    format_result_all = []
    t_test_res_info_all = []
    num_reg_samples_all = []
    edge_num_all = []
    gene_tfs_node_num_all = []
    genes_samples_statistic = []
    
    for k in range(len(target_genes_list)):
        target_gene = target_genes_list[k]
        print(f"Running for gene: {target_gene} ......")
        
        candidate_peak_id = candidate_peaks[candidate_peaks['gene'] == target_gene]['peak'].tolist()
        if len(candidate_peak_id) == 0:
            continue
        
        candidate_peak_id = candidate_peak_id[0]
        
        tfs_file = os.path.join(path_read, "TFs", f"{target_gene}_TFs_list.txt")
        if not os.path.exists(tfs_file):
            print(f"Warning: {tfs_file} not found")
            continue
        
        tfs_list = pd.read_table(tfs_file, header=0).iloc[:, 0].astype(str).tolist()
        
        sample_list_file = os.path.join(path_read, "peak_luad_sample.csv")
        if os.path.exists(sample_list_file):
            sample_list = pd.read_csv(sample_list_file, header=0).iloc[:, 0].astype(str).tolist()
        else:
            sample_list = list(set(tcga_luad.columns) & set(peak_luad.columns))
        
        gtex_temp = gtex_luad.reindex(index=tfs_list).dropna()
        tcga_temp = tcga_luad.reindex(index=tfs_list, columns=sample_list).dropna()
        
        high_exp_samples = []
        down_exp_samples = []
        gene_sample_statistic = []
        
        for sample in sample_list:
            if sample not in tcga_luad.columns:
                continue
            
            target_sample_value = tcga_luad.loc[target_gene, sample]
            gtex_values = gtex_luad.loc[target_gene, :].astype(float).dropna()
            
            if len(gtex_values) < 2:
                continue
            
            t_stat, p_val = ttest_1samp(gtex_values, target_sample_value)
            
            if t_stat < 0 and p_val < 0.001:
                print(f"The expression of {target_gene} in sample: '{sample}' is higher expression than GTEx samples......")
                high_exp_samples.append(sample)
            
            if t_stat > 0 and p_val < 0.001:
                print(f"The expression of {target_gene} in sample: '{sample}' is lower expression than GTEx samples......")
                down_exp_samples.append(sample)
            
            gene_sample_statistic.append({
                'Gene': target_gene,
                'Sample': sample,
                'T': t_stat,
                'P-value': p_val
            })
        
        genes_samples_statistic.extend(gene_sample_statistic)
        
        if target_gene in negative_genes:
            if len(high_exp_samples) <= 0:
                continue
            
            t_test_res_info = []
            num_reg_samples = []
            
            results_tfs_up = {}
            results_rps_down = {}
            
            for sample in high_exp_samples:
                if sample not in peak_luad.columns:
                    continue
                
                up_thres = np.percentile(peak_luad.loc[:, sample].values, 75)
                peak_value = peak_luad.loc[candidate_peak_id, sample] if candidate_peak_id in peak_luad.index else 0
                
                print(f"{peak_value} <------> {up_thres}")
                
                up_reg_tfs = []
                down_reg_tfs = []
                t_test_res_samples = []
                
                if peak_value > up_thres:
                    for tf in tfs_list:
                        if tf not in gtex_temp.index or tf not in tcga_temp.index:
                            continue
                        
                        gtex_tf_values = gtex_temp.loc[tf, :].astype(float).dropna()
                        tcga_tf_value = tcga_temp.loc[tf, sample]
                        
                        if len(gtex_tf_values) < 2:
                            continue
                        
                        t_stat_temp, p_val_temp = ttest_1samp(gtex_tf_values, tcga_tf_value)
                        
                        if t_stat_temp < 0 and p_val_temp < 0.001:
                            up_reg_tfs.append(tf)
                        if t_stat_temp > 0 and p_val_temp < 0.001:
                            down_reg_tfs.append(tf)
                        
                        score = 0
                        if p_val_temp < 0.001 and len(string_data) > 0:
                            matches = string_data[
                                (string_data.iloc[:, 1] == target_gene) & 
                                (string_data.iloc[:, 0] == tf)
                            ]
                            if len(matches) > 0:
                                score = float(matches.iloc[0, 2]) if len(matches.columns) > 2 else 0
                        
                        t_test_res_samples.append({
                            'TF': tf,
                            'Gene': target_gene,
                            'Sample': sample,
                            'T-statistic': t_stat_temp,
                            'P-value': p_val_temp,
                            'Score': score
                        })
                
                if cnv_test:
                    up_reg_tfs = [tf for tf in up_reg_tfs if tf in tcga_cnv.index and tcga_cnv.loc[tf, sample] != 0]
                    down_reg_tfs = [tf for tf in down_reg_tfs if tf in tcga_cnv.index and tcga_cnv.loc[tf, sample] != 0]
                
                results_tfs_up[sample] = up_reg_tfs
                results_rps_down[sample] = down_reg_tfs
                
                if len(t_test_res_samples) > 0:
                    t_test_res_info.extend(t_test_res_samples)
                
                num_high = len(results_tfs_up[sample]) + len(results_rps_down[sample])
                num_reg_samples.append({
                    'Gene': target_gene,
                    'Sample': sample,
                    'Number': num_high
                })
            
            results_regulatory_frame = pd.concat([
                format_regulation(results_tfs_up, target_gene, 1),
                format_regulation(results_rps_down, target_gene, -1)
            ], ignore_index=True)
        else:
            t_test_res_info = []
            num_reg_samples = []
            
            results_tf_down = {}
            results_rp_up = {}
            
            if len(down_exp_samples) > 0:
                for sample in down_exp_samples:
                    if sample not in peak_luad.columns:
                        continue
                    
                    up_thres = np.percentile(peak_luad.loc[:, sample].values, 75)
                    peak_value = peak_luad.loc[candidate_peak_id, sample] if candidate_peak_id in peak_luad.index else 0
                    
                    print(f"{peak_value} <------> {up_thres}")
                    
                    up_reg_tfs = []
                    down_reg_tfs = []
                    t_test_res_samples = []
                    
                    if peak_value < up_thres:
                        for tf in tfs_list:
                            if tf not in gtex_temp.index or tf not in tcga_temp.index:
                                continue
                            
                            gtex_tf_values = gtex_temp.loc[tf, :].astype(float).dropna()
                            tcga_tf_value = tcga_temp.loc[tf, sample]
                            
                            if len(gtex_tf_values) < 2:
                                continue
                            
                            t_stat_temp, p_val_temp = ttest_1samp(gtex_tf_values, tcga_tf_value)
                            
                            if t_stat_temp < 0 and p_val_temp < 0.001:
                                up_reg_tfs.append(tf)
                            if t_stat_temp > 0 and p_val_temp < 0.001:
                                down_reg_tfs.append(tf)
                            
                            score = 0
                            if p_val_temp < 0.001 and len(string_data) > 0:
                                matches = string_data[
                                    (string_data.iloc[:, 1] == target_gene) & 
                                    (string_data.iloc[:, 0] == tf)
                                ]
                                if len(matches) > 0:
                                    score = float(matches.iloc[0, 2]) if len(matches.columns) > 2 else 0
                            
                            t_test_res_samples.append({
                                'TF': tf,
                                'Gene': target_gene,
                                'Sample': sample,
                                'T-statistic': t_stat_temp,
                                'P-value': p_val_temp,
                                'Score': score
                            })
                    
                    if cnv_test:
                        up_reg_tfs = [tf for tf in up_reg_tfs if tf in tcga_cnv.index and tcga_cnv.loc[tf, sample] != 0]
                        down_reg_tfs = [tf for tf in down_reg_tfs if tf in tcga_cnv.index and tcga_cnv.loc[tf, sample] != 0]
                    
                    results_rp_up[sample] = up_reg_tfs
                    results_tf_down[sample] = down_reg_tfs
                    
                    if len(t_test_res_samples) > 0:
                        t_test_res_info.extend(t_test_res_samples)
                    
                    num_down = len(results_tf_down[sample]) + len(results_rp_up[sample])
                    num_reg_samples.append({
                        'Gene': target_gene,
                        'Sample': sample,
                        'Number': num_down
                    })
            
            results_regulatory_frame = pd.concat([
                format_regulation(results_tf_down, target_gene, 1),
                format_regulation(results_rp_up, target_gene, -1)
            ], ignore_index=True)
        
        if len(results_regulatory_frame) > 0 and 'TF' in results_regulatory_frame.columns:
            freq_threshold = 0
            
            tf_counts = results_regulatory_frame['TF'].value_counts()
            sample_count = len(down_exp_samples) if target_gene not in negative_genes else len(high_exp_samples)
            if sample_count == 0:
                sample_count = 1
            
            results_regulatory_freq = pd.DataFrame({
                'TF': tf_counts.index,
                'Freq': tf_counts.values / sample_count
            })
            results_regulatory_freq = results_regulatory_freq[results_regulatory_freq['Freq'] >= freq_threshold]
            
            results_regulatory_table = pd.DataFrame({'TF': results_regulatory_freq['TF']})
            results_regulatory_table = frame_regulatory_table(results_regulatory_table, target_gene, databases_list)
            
            results_regulatory_table.to_csv(
                os.path.join(path_save, "Regulatory_table", f"{target_gene}_regulatory_table.csv"),
                index=False
            )
            
            db_cols = [col for col in results_regulatory_table.columns if col in databases_list]
            hit_mask = results_regulatory_table[db_cols].sum(axis=1) != 0
            hit_tfs = results_regulatory_table.loc[hit_mask, 'TF'].tolist()
            
            dataset_num = 0
            for tf in hit_tfs:
                dataset_num += len(results_regulatory_frame[results_regulatory_frame['TF'] == tf])
            
            type_counts = results_regulatory_frame['Type'].value_counts()
            positive_count = type_counts.get('1', 0) if '1' in type_counts.index else 0
            negative_count = type_counts.get('-1', 0) if '-1' in type_counts.index else 0
            
            edge_num_all.append({
                'Gene': target_gene,
                'Positive': positive_count,
                'Negative': negative_count,
                'ImmuCycReg_framework': len(results_regulatory_frame),
                'Database': dataset_num
            })
            
            with open(os.path.join(path_save, "TFs_list_selected", f"{target_gene}_TFs_list.txt"), 'w') as f:
                f.write("TF\n")
                for tf in hit_tfs:
                    f.write(f"{tf}\n")
            
            gene_tfs_node_num_all.append({
                'Gene': target_gene,
                'ImmuCyc_framework': len(results_regulatory_freq),
                'Database': len(hit_tfs)
            })
            
            if len(t_test_res_info) > 0:
                t_test_res_info_all.extend(t_test_res_info)
            
            if len(num_reg_samples) > 0:
                num_reg_samples_all.extend(num_reg_samples)
            
            format_result = results_regulatory_frame[results_regulatory_frame['TF'].isin(hit_tfs)]
            if len(format_result) > 0:
                format_result_all.append(format_result)
    
    if len(genes_samples_statistic) > 0:
        pd.DataFrame(genes_samples_statistic).to_csv(
            os.path.join(path_save, "ImmuCycRegFramework_Gene-Sample.csv"),
            index=False
        )
    
    if len(t_test_res_info_all) > 0:
        pd.DataFrame(t_test_res_info_all).to_csv(
            os.path.join(path_save, "T-test", "T-test_results_information.csv"),
            index=False
        )
    
    if len(edge_num_all) > 0:
        pd.DataFrame(edge_num_all).to_csv(
            os.path.join(path_save, "Gene_TFs_edge_num.csv"),
            index=False
        )
    
    if len(gene_tfs_node_num_all) > 0:
        pd.DataFrame(gene_tfs_node_num_all).to_csv(
            os.path.join(path_save, "Gene_TFs_node_num.csv"),
            index=False
        )
    
    if len(edge_num_all) > 0:
        edge_df = pd.DataFrame(edge_num_all)
        
        fig, axes = plt.subplots(1, 2, figsize=(15, 3))
        
        net_edge_info = edge_df[['Gene', 'ImmuCycReg_framework', 'Database']].melt(
            id_vars='Gene', var_name='Method', value_name='Count'
        )
        
        sns.barplot(data=net_edge_info, x='Gene', y='Count', hue='Method', ax=axes[0],
                   palette=['white', 'gray'])
        axes[0].set_title('Edge Count')
        axes[0].set_ylabel('Count')
        axes[0].tick_params(axis='x', rotation=45)
        
        net_reg_info = edge_df[['Gene', 'Negative', 'Positive']].melt(
            id_vars='Gene', var_name='Regulation', value_name='Count'
        )
        
        sns.barplot(data=net_reg_info, x='Gene', y='Count', hue='Regulation', ax=axes[1],
                   palette=['#6666ff', '#ff9933'])
        axes[1].set_title('Regulation Count')
        axes[1].set_ylabel('Count')
        axes[1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(path_save, "Figure", "Supplementary Figure 4_bc.pdf"))
        plt.close()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

