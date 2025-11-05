import os
import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from l0learn import L0Learn
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def genie3(matrix_data, n_cores=3, verbose=True):
    """GENIE3 implementation using Random Forest"""
    n_genes = matrix_data.shape[1]
    weight_matrix = np.zeros((n_genes, n_genes))
    
    for i in range(n_genes):
        if verbose and i % 10 == 0:
            print(f"Processing gene {i+1}/{n_genes}")
        
        target = matrix_data[:, i]
        predictors = np.delete(matrix_data, i, axis=1)
        
        rf = RandomForestRegressor(n_estimators=100, max_depth=10, random_state=42, n_jobs=n_cores)
        rf.fit(predictors, target)
        
        importances = rf.feature_importances_
        
        weight_matrix[:i, i] = importances[:i]
        weight_matrix[i+1:, i] = importances[i:]
    
    return weight_matrix


def get_link_list(weight_matrix, gene_names=None):
    """Convert weight matrix to link list"""
    n = weight_matrix.shape[0]
    links = []
    
    for i in range(n):
        for j in range(n):
            if i != j:
                source = gene_names[i] if gene_names is not None else f"Gene{i+1}"
                target = gene_names[j] if gene_names is not None else f"Gene{j+1}"
                weight = weight_matrix[i, j]
                links.append([source, target, weight])
    
    links_df = pd.DataFrame(links, columns=['Source', 'Target', 'Weight'])
    links_df = links_df.sort_values('Weight', ascending=False)
    return links_df


def l0reg_ensemble(expression_matrix, l0_pred_sample_min=20, l0_pred_sample_max=80,
                   l0_exp_sample_min=20, l0_exp_sample_max=80, l0_rank_threshold=5,
                   l0_ensemble_size=100):
    """L0REG ensemble implementation"""
    n_samples, n_genes = expression_matrix.shape
    
    l0_pred_sample_min = int(l0_pred_sample_min * n_genes / 100)
    l0_pred_sample_max = int(l0_pred_sample_max * n_genes / 100)
    l0_exp_sample_min = int(l0_exp_sample_min * n_samples / 100)
    l0_exp_sample_max = int(l0_exp_sample_max * n_samples / 100)
    l0_rank_threshold = int(l0_rank_threshold * n_genes / 100)
    
    weight_matrix = np.zeros((n_genes, n_genes))
    
    for target_idx in range(n_genes):
        target = expression_matrix[:, target_idx]
        predictor_indices = list(range(n_genes))
        predictor_indices.remove(target_idx)
        
        for _ in range(l0_ensemble_size):
            n_pred = np.random.randint(l0_pred_sample_min, l0_pred_sample_max + 1)
            n_exp = np.random.randint(l0_exp_sample_min, l0_exp_sample_max + 1)
            
            selected_predictors = np.random.choice(predictor_indices, size=n_pred, replace=False)
            selected_samples = np.random.choice(n_samples, size=n_exp, replace=False)
            
            X = expression_matrix[np.ix_(selected_samples, selected_predictors)]
            y = target[selected_samples]
            
            try:
                model_l0 = L0Learn.fit(X, y, penalty="L0", maxSuppSize=len(selected_predictors))
                
                model_info_list = []
                for i_idx in range(len(model_l0.lambda_)):
                    for j_idx in range(len(model_l0.gamma_)):
                        supp_size = model_l0.supp_size_[i_idx][j_idx]
                        model_info_list.append({
                            'lambda': model_l0.lambda_[i_idx],
                            'gamma': model_l0.gamma_[j_idx],
                            'suppSize': supp_size
                        })
                
                model_info_df = pd.DataFrame(model_info_list).sort_values('suppSize', ascending=False)
                
                if len(model_info_df) > 0:
                    lambda_l0 = model_info_df.iloc[0]['lambda']
                    gamma_l0 = model_info_df.iloc[0]['gamma']
                    
                    lambda_idx = model_l0.lambda_.index(lambda_l0) if lambda_l0 in model_l0.lambda_ else 0
                    gamma_idx = model_l0.gamma_.index(gamma_l0) if gamma_l0 in model_l0.gamma_ else 0
                    
                    coefs = model_l0.coef(lambda_idx, gamma_idx)
                    coefs = np.array(coefs).ravel()
                    
                    top_indices = np.argsort(np.abs(coefs))[-l0_rank_threshold:]
                    
                    for idx in top_indices:
                        pred_idx = selected_predictors[idx]
                        weight_matrix[pred_idx, target_idx] += 1
            except:
                continue
    
    return weight_matrix


def convert_sorted_rank_tsv_to_adj_matrix(tsv_file_or_df):
    """Convert sorted rank TSV to adjacency matrix"""
    if isinstance(tsv_file_or_df, str):
        df = pd.read_csv(tsv_file_or_df, sep="\t", header=None)
    else:
        df = tsv_file_or_df.copy()
    
    if len(df.columns) < 2:
        return pd.DataFrame()
    
    sources = df.iloc[:, 0].astype(str).unique()
    targets = df.iloc[:, 1].astype(str).unique()
    
    adj_matrix = pd.DataFrame(0, index=sources, columns=targets)
    
    for _, row in df.iterrows():
        source = str(row.iloc[0])
        target = str(row.iloc[1])
        if source in adj_matrix.index and target in adj_matrix.columns:
            adj_matrix.loc[source, target] = 1
    
    return adj_matrix


def prepare_eval(prediction_tsv, gold_tsv, total_predictions_accepted=100000):
    """Prepare evaluation object"""
    if isinstance(prediction_tsv, str):
        pred = pd.read_csv(prediction_tsv, sep="\t", header=None)
    else:
        pred = prediction_tsv.copy()
    
    gold = convert_sorted_rank_tsv_to_adj_matrix(gold_tsv)
    
    if len(pred) > total_predictions_accepted:
        pred = pred.iloc[:total_predictions_accepted, :]
    
    pred = pred[pred.iloc[:, 0] != pred.iloc[:, 1]].copy()
    
    p = gold.values.sum()
    n = gold.size - p - sum(gold.index.isin(gold.columns))
    t = p + n
    
    first_row = pred.iloc[:, 0].astype(str).values
    second_row = pred.iloc[:, 1].astype(str).values
    third_row = np.zeros(len(first_row), dtype=int)
    rank = np.zeros(len(first_row), dtype=int)
    
    correct = 0
    incorrect = 0
    
    for i in range(len(first_row)):
        source = first_row[i]
        target = second_row[i]
        if source in gold.index and target in gold.columns:
            if gold.loc[source, target] == 1:
                correct += 1
                third_row[i] = 1
                rank[i] = correct
            else:
                incorrect -= 1
                third_row[i] = 0
                rank[i] = incorrect
    
    if len(first_row) < t:
        odds = (p - correct) / (t - len(first_row))
    else:
        odds = 0
    
    random_positive = np.full(t - len(first_row), odds)
    random_negative = np.full(t - len(first_row), 1 - odds)
    
    positive = np.concatenate([third_row, random_positive])
    negative = np.concatenate([np.abs(third_row - 1), random_negative])
    
    tpk = np.cumsum(positive)
    fpk = np.cumsum(negative)
    
    k = np.arange(1, t + 1)
    
    tpr = tpk / p if p > 0 else np.zeros(len(tpk))
    fpr = fpk / n if n > 0 else np.zeros(len(fpk))
    rec = tpr
    prec = tpk / k
    
    return {
        'tpr': tpr,
        'fpr': fpr,
        'prec': prec,
        'rec': rec,
        'p': p,
        'n': n,
        'tpk': tpk,
        'fpk': fpk,
        'predictors': first_row,
        'targets': second_row,
        'rank': rank
    }


def calc_auroc(eval_obj):
    """Calculate AUROC"""
    return trapz(eval_obj['fpr'], eval_obj['tpr'])


def calc_aupr(eval_obj):
    """Calculate AUPR"""
    if eval_obj['p'] == 0:
        return 0
    return trapz(eval_obj['rec'], eval_obj['prec']) / (1 - 1 / eval_obj['p'])


def convert_adj_matrix_to_sorted_rank_tsv(adj_matrix, output_file=None, desc=False):
    """Convert adjacency matrix to sorted rank TSV"""
    sources = adj_matrix.index.tolist()
    targets = adj_matrix.columns.tolist()
    
    links = []
    for source in sources:
        for target in targets:
            if source != target:
                weight = adj_matrix.loc[source, target]
                links.append([source, target, weight])
    
    links_df = pd.DataFrame(links, columns=['Source', 'Target', 'Weight'])
    links_df = links_df.sort_values('Weight', ascending=not desc)
    
    if output_file:
        links_df.to_csv(output_file, sep="\t", header=False, index=False)
    
    return links_df


def bar_plot(data, title="", y_title="", bar_color=None):
    """Create bar plot"""
    fig, ax = plt.subplots(figsize=(8, 4))
    
    if bar_color is None:
        bar_color = ['white', 'gray']
    
    if isinstance(data, pd.DataFrame):
        if 'Dataset' in data.columns and 'Method' in data.columns:
            groups = data['Dataset'].unique()
            methods = data['Method'].unique()
            
            x = np.arange(len(groups))
            width = 0.8 / len(methods)
            
            for i, method in enumerate(methods):
                method_data = data[data['Method'] == method]
                values = []
                for group in groups:
                    group_data = method_data[method_data['Dataset'] == group]
                    if len(group_data) > 0:
                        values.append(group_data.iloc[0, 2])
                    else:
                        values.append(0)
                
                offset = (i - len(methods) / 2) * width + width / 2
                color = bar_color[i % len(bar_color)]
                ax.bar(x + offset, values, width, label=method, color=color, edgecolor='black')
            
            ax.set_xlabel('Dataset')
            ax.set_ylabel(y_title)
            ax.set_title(title)
            ax.set_xticks(x)
            ax.set_xticklabels(groups, rotation=45, ha='right')
            ax.legend()
            ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    return fig


def box_plot(data, box_color=None):
    """Create box plot"""
    fig, ax = plt.subplots(figsize=(6, 4))
    
    if box_color is None:
        box_color = ['white', 'gray']
    
    if isinstance(data, pd.DataFrame):
        if 'Method' in data.columns and 'AUROC' in data.columns:
            methods = data['Method'].unique()
            positions = np.arange(len(methods))
            
            for i, method in enumerate(methods):
                method_data = data[data['Method'] == method]['AUROC'].values
                bp = ax.boxplot([method_data], positions=[positions[i]], widths=0.6,
                               patch_artist=True, labels=[method])
                for patch in bp['boxes']:
                    patch.set_facecolor(box_color[i % len(box_color)])
    
    ax.set_ylabel('AUROC')
    ax.set_title('AUROC Comparison')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    return fig


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    grn_path = os.path.join(path_save, "GRN")
    os.makedirs(grn_path, exist_ok=True)
    os.makedirs(os.path.join(path_save, "Figure"), exist_ok=True)
    
    evaluation_gnws_list = []
    
    for n in range(1, 3):
        evaluation_gnws = []
        
        if n == 1:
            dataset_select = True
            dataset_num = 5
        else:
            dataset_select = False
            dataset_num = 10
        
        for i in range(1, dataset_num + 1):
            if dataset_select:
                path_way = os.path.join(path_read, "GRN_5_gold_networks", str(i))
                matrix_way = os.path.join(path_way, "measure", "measure.tsv")
                gold_way = os.path.join(path_way, "gold", "gold.tsv")
            else:
                path_way = os.path.join(path_read, "GRN_10_gold_networks", str(i))
                matrix_way = os.path.join(path_way, "Ecoli-20_nonoise_multifactorial.tsv")
                gold_way = os.path.join(path_way, "Ecoli-20_goldstandard.tsv")
            
            if not os.path.exists(matrix_way) or not os.path.exists(gold_way):
                print(f"Warning: Files not found for network {i}")
                continue
            
            matrix_data = pd.read_csv(matrix_way, sep="\t", header=0)
            matrix_data = matrix_data.T.values
            
            l0_output_file = os.path.join(grn_path, f"GNW{dataset_num}_L0_Net{i}.txt")
            
            print(f"Running L0REG for network {i}...")
            weight_matrix_l0 = l0reg_ensemble(
                matrix_data,
                l0_pred_sample_min=20,
                l0_pred_sample_max=80,
                l0_exp_sample_min=20,
                l0_exp_sample_max=80,
                l0_rank_threshold=5,
                l0_ensemble_size=100
            )
            
            gene_names = [f"Gene{j+1}" for j in range(matrix_data.shape[1])]
            adj_matrix_l0 = pd.DataFrame(weight_matrix_l0, index=gene_names, columns=gene_names)
            links_df_l0 = convert_adj_matrix_to_sorted_rank_tsv(adj_matrix_l0, desc=False)
            links_df_l0.to_csv(l0_output_file, sep="\t", header=False, index=False)
            
            print(f"Evaluating L0REG for network {i}...")
            eval_obj_l0 = prepare_eval(l0_output_file, gold_way)
            auroc_l0 = calc_auroc(eval_obj_l0)
            auprc_l0 = calc_aupr(eval_obj_l0)
            
            print(f"Running GENIE3 for network {i}...")
            weight_matrix_genie3 = genie3(matrix_data, n_cores=3, verbose=True)
            
            gene_names = [f"Gene{j+1}" for j in range(matrix_data.shape[1])]
            links_df_genie3 = get_link_list(weight_matrix_genie3, gene_names)
            
            genie3_output_file = os.path.join(grn_path, f"GNW{dataset_num}_GENIE3_Net{i}.txt")
            links_df_genie3.to_csv(genie3_output_file, sep="\t", header=False, index=False)
            
            print(f"Evaluating GENIE3 for network {i}...")
            eval_obj_genie3 = prepare_eval(genie3_output_file, gold_way)
            auroc_genie3 = calc_auroc(eval_obj_genie3)
            auprc_genie3 = calc_aupr(eval_obj_genie3)
            
            evaluation_gnws.append({
                'Dataset': f'Net{i}',
                'L0Reg_framework': auroc_l0,
                'GENIE3': auroc_genie3
            })
        
        if len(evaluation_gnws) > 0:
            eval_df = pd.DataFrame(evaluation_gnws)
            eval_df.to_csv(os.path.join(path_read, f"AUROC_Net{dataset_num}.csv"), index=False)
            evaluation_gnws_list.append(eval_df)
    
    if len(evaluation_gnws_list) == 2:
        net5_auroc = evaluation_gnws_list[0].melt(id_vars='Dataset', var_name='Method', value_name='AUROC')
        net10_auroc = evaluation_gnws_list[1].melt(id_vars='Dataset', var_name='Method', value_name='AUROC')
        all_net_data = pd.concat([net5_auroc, net10_auroc], ignore_index=True).drop(columns=['Dataset'])
        
        color = ['white', 'gray']
        
        p1 = bar_plot(net5_auroc, y_title="AUROC", bar_color=color)
        p2 = bar_plot(net10_auroc, y_title="AUROC", bar_color=color)
        p3 = box_plot(all_net_data, box_color=color)
        
        fig_combined, axes = plt.subplots(1, 3, figsize=(18, 4))
        axes[0].remove()
        axes[0] = p1.axes[0]
        axes[1].remove()
        axes[1] = p2.axes[0]
        axes[2].remove()
        axes[2] = p3.axes[0]
        
        plt.tight_layout()
        plt.savefig(os.path.join(path_save, "Figure", "Supplementary Figure 7.pdf"))
        plt.close()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

