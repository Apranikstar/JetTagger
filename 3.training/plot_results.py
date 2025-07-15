import uproot
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_score_distribution(labels, scores, filename=None):
    plt.figure(figsize=(8,6))
    plt.hist(scores[labels==1], bins=50, alpha=0.6, label='Signal (top)', color='r', density=True)
    plt.hist(scores[labels==0], bins=50, alpha=0.6, label='Background', color='b', density=True)
    plt.xlabel('Score')
    plt.ylabel('Density')
    plt.legend()
    plt.title('Score Distribution')
    if filename:
        plt.savefig(filename)
    else:
        plt.show()

def plot_roc_curve(labels, scores, filename=None):
    from sklearn.metrics import roc_curve, auc
    fpr, tpr, _ = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)
    plt.figure(figsize=(8,6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
    plt.plot([0,1], [0,1], color='navy', lw=1, linestyle='--')
    plt.xlim([0.0,1.0])
    plt.ylim([0.0,1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    if filename:
        plt.savefig(filename)
    else:
        plt.show()

def plot_correlation_matrix(corr_data, variables, filename=None):
    data = np.vstack([corr_data[var] for var in variables])
    corr_matrix = np.corrcoef(data)
    plt.figure(figsize=(12,10))
    sns.heatmap(corr_matrix, xticklabels=variables, yticklabels=variables,
                cmap='coolwarm', annot=True, fmt=".2f", square=True)
    plt.title('Correlation Matrix')
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
    else:
        plt.show()

if __name__ == "__main__":
    file_path = "test_output.root"
    with uproot.open(file_path) as file:
        tree = file["Events;1"]
        
        # Read all needed arrays inside with block
        labels = tree["_label_"].array()
        score_top = tree["score_recojet_isT"].array()
        score_j = tree["score_recojet_isJ"].array()
        jet_p = tree["jet_p"].array()
        jet_mass = tree["jet_mass"].array()
        jet_phi = tree["jet_phi"].array()
        jet_theta = tree["jet_theta"].array()
        
        variables_to_corr = [
            'jet_p', 'jet_mass', 'jet_phi', 'jet_theta',
            'jet_nconst', 'jet_nmu', 'jet_nel', 'jet_nchad',
            'jet_ngamma', 'jet_nnhad', 'jet_npfcand',
            'score_recojet_isT', 'score_recojet_isJ'
        ]
        corr_data = {var: tree[var].array() for var in variables_to_corr}
    
    # Now plotting — file is closed but data is in memory
    
    plot_score_distribution(labels, score_top, filename="score_distribution.png")
    plot_roc_curve(labels, score_top, filename="roc_curve.png")
    plot_correlation_matrix(corr_data, variables_to_corr, filename="correlation_matrix.png")
    
    print("Plots saved: score_distribution.png, roc_curve.png, correlation_matrix.png")
