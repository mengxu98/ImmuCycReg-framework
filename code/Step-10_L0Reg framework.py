import os
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr, t as t_dist
from l0learn import L0Learn


def evaluate_model(y_true, y_pred):
    if len(y_true) != len(y_pred):
        return [0, 0, 0, 1]
    
    y_true = np.array(y_true).ravel()
    y_pred = np.array(y_pred).ravel()
    
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() < 2:
        return [0, 0, 0, 1]
    
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    rmse = np.sqrt(np.mean((y_true - y_pred) ** 2))
    
    corr, pval = pearsonr(y_true, y_pred) if len(y_true) > 1 else (0, 1)
    
    return [r2, rmse, corr, pval]


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/", "L0Reg-framework"))
    os.makedirs(path_save, exist_ok=True)
    
    penalty = "L0"
    
    tcga_luad = pd.read_csv(os.path.join(os.path.dirname(path_save), "TCGA-LUAD-tpm-tcga-t.csv"), index_col=0)
    scaler = StandardScaler()
    tcga_scale = pd.DataFrame(
        scaler.fit_transform(tcga_luad.T).T,
        index=tcga_luad.index,
        columns=tcga_luad.columns
    )
    
    samples_cluster = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"), header=None, index_col=0)
    
    genes_17_file = os.path.join(path_read, "Genes_17.txt")
    if os.path.exists(genes_17_file):
        gene_list = pd.read_table(genes_17_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Error: Genes_17.txt not found")
        return 1
    
    trn_cluster_filter_all = []
    evaluate_r2_all = []
    
    for j in range(1, 5):
        print(f"Running for cluster {j} ......")
        cluster_dir = os.path.join(path_save, f"cluster{j}")
        os.makedirs(cluster_dir, exist_ok=True)
        
        cluster_samples = samples_cluster[samples_cluster.iloc[:, 0] == j].index.tolist()
        cluster_samples = [s for s in cluster_samples if s in tcga_scale.columns]
        
        if len(cluster_samples) == 0:
            continue
        
        selected_tfs = []
        trn_cluster_filter = []
        evaluate_r2 = []
        
        for i, target_gene in enumerate(gene_list):
            tfs_file = os.path.join(path_read, "TFs", f"{target_gene}_TFs_list.txt")
            if not os.path.exists(tfs_file):
                print(f"No {target_gene} TF file......")
                continue
            
            print(f"Running for cluster: {j}, NO. {i+1} gene: {target_gene}......")
            
            tfs_list = pd.read_table(tfs_file, header=0).iloc[:, 0].astype(str).tolist()
            
            common_tfs = [tf for tf in tfs_list if tf in tcga_scale.index]
            common_samples = [s for s in cluster_samples if s in tcga_scale.columns]
            
            if len(common_tfs) == 0 or len(common_samples) == 0:
                continue
            
            test_data_x = tcga_scale.loc[common_tfs, common_samples].T.values
            test_data_y = tcga_scale.loc[target_gene, common_samples].values.reshape(-1, 1)
            
            try:
                model_l0 = L0Learn.fit(test_data_x, test_data_y, penalty=penalty, max_support_size=test_data_x.shape[1])
                
                model_info_list = []
                for lambda_idx in range(len(model_l0.lambda_)):
                    for gamma_idx in range(len(model_l0.gamma_)):
                        coefs = model_l0.coef_[lambda_idx][gamma_idx]
                        supp_size = np.sum(np.abs(coefs) > 1e-10)
                        model_info_list.append({
                            'lambda': model_l0.lambda_[lambda_idx],
                            'gamma': model_l0.gamma_[gamma_idx],
                            'suppSize': supp_size
                        })
                
                model_info_df = pd.DataFrame(model_info_list)
                model_info_df = model_info_df.sort_values('suppSize', ascending=False)
                
                if len(model_info_df) == 0:
                    continue
                
                lambda_l0 = model_info_df.iloc[0]['lambda']
                gamma_l0 = model_info_df.iloc[0]['gamma']
                
                lambda_idx = np.where(model_l0.lambda_ == lambda_l0)[0][0]
                gamma_idx = np.where(model_l0.gamma_ == gamma_l0)[0][0]
                
                coefs = model_l0.coef_[lambda_idx][gamma_idx].flatten()
                data_y_pre_l0 = model_l0.predict(test_data_x, lambda_l0, gamma_l0).flatten()
                
                lm_fit = LinearRegression()
                lm_fit.fit(test_data_x, test_data_y.ravel())
                
                lm_residuals = test_data_y.ravel() - lm_fit.predict(test_data_x)
                lm_mse = np.mean(lm_residuals ** 2)
                
                if lm_mse > 0:
                    lm_se = np.sqrt(lm_mse / (len(test_data_y) - len(common_tfs) - 1))
                    lm_t_stats = lm_fit.coef_ / (lm_se + 1e-10)
                    df = len(test_data_y) - len(common_tfs) - 1
                    lm_pvals = 2 * (1 - t_dist.cdf(np.abs(lm_t_stats), df))
                else:
                    lm_t_stats = np.zeros(len(common_tfs))
                    lm_pvals = np.ones(len(common_tfs))
                
                lm_fit_summ = pd.DataFrame({
                    'Estimate': coefs,
                    'T-statistic': lm_t_stats,
                    'P-value': lm_pvals,
                    'Pr(>|t|)': lm_pvals
                }, index=common_tfs)
                
                lm_fit_summ['Weight'] = np.abs(lm_fit_summ['Estimate'])
                weight_sum = lm_fit_summ['Weight'].sum()
                if weight_sum > 0:
                    lm_fit_summ['Weight'] = lm_fit_summ['Weight'] / weight_sum
                else:
                    lm_fit_summ['Weight'] = 0
                
                lm_fit_summ['Regulation'] = 1
                lm_fit_summ.loc[lm_fit_summ['Estimate'] < 0, 'Regulation'] = -1
                
                lm_fit_summ_filter = lm_fit_summ[lm_fit_summ['Pr(>|t|)'] < 0.05].copy()
                
                evaluate_result_l0 = evaluate_model(test_data_y.ravel(), data_y_pre_l0)
                
                if (evaluate_result_l0[3] < 0.05 and
                    evaluate_result_l0[2] > 0 and
                    len(lm_fit_summ_filter) > 0):
                    
                    single_gene_trn_filter = pd.DataFrame({
                        'TF': lm_fit_summ_filter.index,
                        'Gene': target_gene,
                        'Weight': lm_fit_summ_filter['Weight'].values,
                        'Regulation': lm_fit_summ_filter['Regulation'].values,
                        'Cluster': f"Cluster{j}"
                    })
                    
                    trn_cluster_filter.append(single_gene_trn_filter)
                    selected_tfs.extend(lm_fit_summ_filter.index.tolist())
                    
                    evaluate_r2.append({
                        'Gene': target_gene,
                        'Cluster': f"Cluster{j}",
                        'L0_Rsquare': evaluate_result_l0[0]
                    })
            
            except Exception as e:
                print(f"Error processing {target_gene}: {e}")
                continue
        
        if len(trn_cluster_filter) > 0:
            trn_cluster_filter_all.extend(trn_cluster_filter)
        if len(evaluate_r2) > 0:
            evaluate_r2_all.extend(evaluate_r2)
        
        if len(selected_tfs) > 0:
            unique_tfs = pd.DataFrame({'TF': list(set(selected_tfs))})
            unique_tfs.to_csv(os.path.join(cluster_dir, "TFs_list.csv"), index=False)
        
        print(f"Cluster {j} done......")
    
    if len(trn_cluster_filter_all) > 0:
        trn_all_df = pd.concat(trn_cluster_filter_all, ignore_index=True)
        trn_all_df.to_csv(os.path.join(path_save, "TRN_cluster_filter_all.csv"), index=False)
    
    if len(evaluate_r2_all) > 0:
        r2_all_df = pd.DataFrame(evaluate_r2_all)
        r2_all_df.to_csv(os.path.join(path_save, "TF_regulatory_value.csv"), index=False)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
