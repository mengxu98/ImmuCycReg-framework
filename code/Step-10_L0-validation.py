import os
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
from l0learn import L0Learn
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def rmse(y_true, y_pred):
    y_true = np.array(y_true).ravel()
    y_pred = np.array(y_pred).ravel()
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    if mask.sum() < 2:
        return np.inf
    return np.sqrt(np.mean((y_true[mask] - y_pred[mask]) ** 2))


def randomize_matrix(X):
    X_rand = X.copy()
    for col in range(X_rand.shape[1]):
        np.random.shuffle(X_rand[:, col])
    return X_rand


def scatter_plot(data, title="", x_title="", y_title="", method="loess"):
    fig, ax = plt.subplots(figsize=(4, 4))
    
    ax.scatter(data['Raw'], data['Pre'], alpha=0.6, s=30)
    
    if method == "lm":
        z = np.polyfit(data['Raw'], data['Pre'], 1)
        p = np.poly1d(z)
        x_line = np.linspace(data['Raw'].min(), data['Raw'].max(), 100)
        ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2)
    
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


def main():
    path_read = os.path.abspath(os.path.join(os.path.dirname(__file__), "data/"))
    path_save = os.path.abspath(os.path.join(os.path.dirname(__file__), "results/"))
    validation_path = os.path.join(path_save, "Validation")
    os.makedirs(validation_path, exist_ok=True)
    os.makedirs(os.path.join(path_save, "Figures"), exist_ok=True)
    
    tcga_luad = pd.read_csv(os.path.join(path_save, "TCGA-LUAD.csv"), index_col=0)
    scaler = StandardScaler()
    tcga_scale = pd.DataFrame(
        scaler.fit_transform(tcga_luad.T).T,
        index=tcga_luad.index,
        columns=tcga_luad.columns
    )
    
    samples_cluster = pd.read_csv(os.path.join(path_read, "sample_cluster_4.csv"), header=None).iloc[:, 0].tolist()
    
    genes_17_file = os.path.join(path_read, "Genes_17.txt")
    if os.path.exists(genes_17_file):
        gene_list = pd.read_table(genes_17_file, header=0).iloc[:, 0].astype(str).tolist()
    else:
        print("Error: Genes_17.txt not found")
        return 1
    
    conditions = ["Random", "Randomize", "TCGA_LUAD"]
    
    evaluate_result_all = []
    result_train_cluster_list = []
    result_test_cluster_list = []
    
    for c, condition in enumerate(conditions):
        condition_path = os.path.join(validation_path, condition)
        os.makedirs(condition_path, exist_ok=True)
        
        result_test_cluster = []
        result_train_cluster = []
        correction_plot_list_test = []
        correction_plot_list_train = []
        
        for i, target_gene in enumerate(gene_list):
            tfs_file = os.path.join(path_read, "TFs", f"{target_gene}_TFs_list.txt")
            if not os.path.exists(tfs_file):
                print(f"The TFs file of {target_gene} not found......")
                continue
            
            tfs_list = pd.read_table(tfs_file, header=0).iloc[:, 0].astype(str).tolist()
            
            common_tfs = [tf for tf in tfs_list if tf in tcga_scale.index]
            common_samples = [s for s in samples_cluster if s in tcga_scale.columns]
            
            if len(common_tfs) == 0 or len(common_samples) == 0:
                continue
            
            X = tcga_scale.loc[common_tfs, common_samples].T.values
            Y = tcga_scale.loc[target_gene, common_samples].values.ravel()
            
            if condition == "TCGA_LUAD":
                print(f"Using '{condition}' data......")
            elif condition == "Randomize":
                print(f"Using '{condition}' data......")
                X = randomize_matrix(X)
                Y = randomize_matrix(Y.reshape(-1, 1)).ravel()
            elif condition == "Random":
                print(f"Using '{condition}' data......")
                Y = np.random.normal(0, 1, X.shape[0])
                for r in range(X.shape[1]):
                    X[:, r] = np.random.normal(0, 1, X.shape[0])
            
            np.random.seed(2022)
            
            train_idx, test_idx = train_test_split(range(len(X)), test_size=0.3, random_state=2022)
            train_data_x = X[train_idx, :]
            train_data_y = Y[train_idx]
            test_data_x = X[test_idx, :]
            test_data_y = Y[test_idx]
            
            try:
                model_l0 = L0Learn.fit(train_data_x, train_data_y.reshape(-1, 1), 
                                      penalty="L0", maxSuppSize=train_data_x.shape[1])
                
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
                
                train_data_y_pre = model_l0.predict(train_data_x, lambda_l0, gamma_l0).flatten()
                rmse_l0_train = rmse(train_data_y, train_data_y_pre)
                
                data_frame_train = pd.DataFrame({
                    'Raw': train_data_y,
                    'Pre': train_data_y_pre
                })
                
                cor_data_train, pval_train = pearsonr(data_frame_train['Raw'], data_frame_train['Pre'])
                
                correction_plot_list_train.append(
                    scatter_plot(data_frame_train,
                               title=f"Gene: {target_gene}",
                               x_title="True expression",
                               y_title="Prediction expression")
                )
                
                test_data_y_pre = model_l0.predict(test_data_x, lambda_l0, gamma_l0).flatten()
                rmse_l0_test = rmse(test_data_y, test_data_y_pre)
                
                data_frame_test = pd.DataFrame({
                    'Raw': test_data_y,
                    'Pre': test_data_y_pre
                })
                
                cor_data_test, pval_test = pearsonr(data_frame_test['Raw'], data_frame_test['Pre'])
                
                correction_plot_list_test.append(
                    scatter_plot(data_frame_test,
                               title=f"Gene: {target_gene}",
                               x_title="True expression",
                               y_title="Prediction expression")
                )
                
                result_train = {
                    'Gene': target_gene,
                    'Corr': cor_data_train,
                    'RMSD': rmse_l0_train,
                    'pval': pval_train
                }
                result_train_cluster.append(result_train)
                
                result_test = {
                    'Gene': target_gene,
                    'Corr': cor_data_test,
                    'RMSD': rmse_l0_test,
                    'pval': pval_test
                }
                result_test_cluster.append(result_test)
            
            except Exception as e:
                print(f"Error processing {target_gene}: {e}")
                continue
        
        if len(correction_plot_list_test) > 0:
            pdf_path_test = os.path.join(condition_path, f"Cor_{condition}_test.pdf")
            with PdfPages(pdf_path_test) as pdf:
                for fig in correction_plot_list_test:
                    pdf.savefig(fig, bbox_inches='tight')
                plt.close('all')
        
        if len(correction_plot_list_train) > 0:
            pdf_path_train = os.path.join(condition_path, f"Cor_{condition}_train.pdf")
            with PdfPages(pdf_path_train) as pdf:
                for fig in correction_plot_list_train:
                    pdf.savefig(fig, bbox_inches='tight')
                plt.close('all')
        
        result_train_df = pd.DataFrame(result_train_cluster)
        result_test_df = pd.DataFrame(result_test_cluster)
        
        result_train_filter = result_train_df[
            (result_train_df['Corr'] > 0.7) & (result_train_df['RMSD'] < 2)
        ].copy()
        
        result_test_filter = result_test_df[
            (result_test_df['Corr'] > 0.7) & (result_test_df['RMSD'] < 2)
        ].copy()
        
        result_test_df.to_csv(os.path.join(condition_path, f"{condition}_Cor-RMSD.csv"), index=False)
        result_test_filter.to_csv(os.path.join(condition_path, f"{condition}_Cor-RMSD_filter.csv"), index=False)
        
        evaluate_result = {
            'Condition': condition,
            'train_Corr': result_train_df['Corr'].mean() if len(result_train_df) > 0 else 0,
            'train_RMSD': result_train_df['RMSD'].mean() if len(result_train_df) > 0 else 0,
            'test_Corr': result_test_df['Corr'].mean() if len(result_test_df) > 0 else 0,
            'test_RMSD': result_test_df['RMSD'].mean() if len(result_test_df) > 0 else 0,
            'train_0.7_Corr': result_train_filter['Corr'].mean() if len(result_train_filter) > 0 else 0,
            'train_0.7_RMSD': result_train_filter['RMSD'].mean() if len(result_train_filter) > 0 else 0,
            'test_0.7_Corr': result_test_filter['Corr'].mean() if len(result_test_filter) > 0 else 0,
            'test_0.7_RMSD': result_test_filter['RMSD'].mean() if len(result_test_filter) > 0 else 0
        }
        evaluate_result_all.append(evaluate_result)
        
        result_train_cluster_list.append(result_train_df)
        result_test_cluster_list.append(result_test_df)
    
    if len(evaluate_result_all) > 0:
        pd.DataFrame(evaluate_result_all).to_csv(
            os.path.join(validation_path, "Cor-RMSD.csv"), index=False
        )
    
    if len(result_train_cluster_list) == 3 and len(result_test_cluster_list) == 3:
        cor_plot_data_train = pd.DataFrame({
            'Gene': result_train_cluster_list[0]['Gene'].values,
            'Random': result_train_cluster_list[0]['Corr'].values,
            'Randomize': result_train_cluster_list[1]['Corr'].values,
            'TCGA_LUAD': result_train_cluster_list[2]['Corr'].values
        }).melt(id_vars='Gene', var_name='Condition', value_name='Correlation')
        
        cor_plot_data_test = pd.DataFrame({
            'Gene': result_test_cluster_list[0]['Gene'].values,
            'Random': result_test_cluster_list[0]['Corr'].values,
            'Randomize': result_test_cluster_list[1]['Corr'].values,
            'TCGA_LUAD': result_test_cluster_list[2]['Corr'].values
        }).melt(id_vars='Gene', var_name='Condition', value_name='Correlation')
        
        rmsd_plot_data_train = pd.DataFrame({
            'Gene': result_train_cluster_list[0]['Gene'].values,
            'Random': result_train_cluster_list[0]['RMSD'].values,
            'Randomize': result_train_cluster_list[1]['RMSD'].values,
            'TCGA_LUAD': result_train_cluster_list[2]['RMSD'].values
        }).melt(id_vars='Gene', var_name='Condition', value_name='RMSD')
        
        rmsd_plot_data_test = pd.DataFrame({
            'Gene': result_test_cluster_list[0]['Gene'].values,
            'Random': result_test_cluster_list[0]['RMSD'].values,
            'Randomize': result_test_cluster_list[1]['RMSD'].values,
            'TCGA_LUAD': result_test_cluster_list[2]['RMSD'].values
        }).melt(id_vars='Gene', var_name='Condition', value_name='RMSD')
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 7))
        
        sns.barplot(data=cor_plot_data_train, x='Gene', y='Correlation', hue='Condition', ax=axes[0, 0])
        axes[0, 0].set_title('Train')
        axes[0, 0].set_ylabel('Correlation')
        axes[0, 0].tick_params(axis='x', rotation=45)
        
        sns.barplot(data=cor_plot_data_test, x='Gene', y='Correlation', hue='Condition', ax=axes[0, 1])
        axes[0, 1].set_title('Test')
        axes[0, 1].set_ylabel('Correlation')
        axes[0, 1].tick_params(axis='x', rotation=45)
        
        sns.barplot(data=rmsd_plot_data_train, x='Gene', y='RMSD', hue='Condition', ax=axes[1, 0])
        axes[1, 0].set_title('Train')
        axes[1, 0].set_ylabel('RMSD')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        sns.barplot(data=rmsd_plot_data_test, x='Gene', y='RMSD', hue='Condition', ax=axes[1, 1])
        axes[1, 1].set_title('Test')
        axes[1, 1].set_ylabel('RMSD')
        axes[1, 1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(validation_path, "Contrast_correlation&RMSD.pdf"), dpi=600)
        plt.close()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
