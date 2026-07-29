import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

adata = sc.read_h5ad("data/processed/GSE138852_annotated.h5ad")
donors = adata.obs["donor_id"].unique().tolist()
cell_types = sorted(adata.obs["cell_type"].unique().tolist())

comp_rows, expr_rows, labels = [], [], []
for donor in donors:
    sub = adata[adata.obs["donor_id"] == donor]
    props = sub.obs["cell_type"].value_counts(normalize=True).reindex(cell_types, fill_value=0.0)
    comp_rows.append(props.values)
    counts = np.asarray(sub.layers["counts"].sum(axis=0)).ravel()
    cpm = counts / max(counts.sum(), 1) * 1e6
    expr_rows.append(np.log1p(cpm))
    labels.append(sub.obs["diagnosis"].iloc[0])

comp_df = pd.DataFrame(comp_rows, index=donors, columns=[f"comp__{c}" for c in cell_types])
expr_df = pd.DataFrame(expr_rows, index=donors, columns=[f"expr__{g}" for g in adata.var_names])
expr_df = expr_df.loc[:, expr_df.std(axis=0) > 0]
X = pd.concat([comp_df, expr_df], axis=1)
y = pd.Series(labels, index=donors, name="diagnosis")
print("Feature matrix:", X.shape, "| donors:", donors)
print("Labels:", y.to_dict())

def make_pipeline(k_best):
    return Pipeline([
        ("scale", StandardScaler()),
        ("select", SelectKBest(score_func=f_classif, k=k_best)),
        ("clf", LogisticRegression(C=0.5, class_weight="balanced", max_iter=2000)),
    ])

def evaluate_loo(X, y, k_best=30):
    y_bin = (y == "AD").astype(int).values
    pipe = make_pipeline(min(k_best, X.shape[1]))
    proba = cross_val_predict(pipe, X.values, y_bin, cv=LeaveOneOut(), method="predict_proba")[:, 1]
    return {"accuracy": float(((proba >= 0.5).astype(int) == y_bin).mean()),
            "roc_auc": float(roc_auc_score(y_bin, proba)), "proba": proba}

result = evaluate_loo(X, y)
print("\n--- Leave-One-Library-Out (n=6) ---")
print("Accuracy:", result["accuracy"], "| ROC-AUC:", result["roc_auc"])
print(pd.DataFrame({"donor": X.index, "true_diagnosis": y.values, "predicted_proba_AD": result["proba"]}))

rng = np.random.default_rng(0)
null_aucs = []
for _ in range(500):
    y_shuf = pd.Series(rng.permutation(y.values), index=y.index)
    try:
        null_aucs.append(evaluate_loo(X, y_shuf)["roc_auc"])
    except ValueError:
        continue
null_aucs = np.array(null_aucs)
p_value = float((null_aucs >= result["roc_auc"]).mean())
print(f"\n--- Permutation null (500 label shuffles) ---")
print(f"Observed AUC: {result['roc_auc']:.3f} | Null mean: {null_aucs.mean():.3f} +/- {null_aucs.std():.3f}")
print(f"P-value: {p_value:.3f}")
