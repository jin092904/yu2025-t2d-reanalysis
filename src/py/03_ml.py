"""Stage 3 (Python): R 의 LASSO(glmnet) · SVM-RFE(caret rfe + kernlab) · RF(randomForest) 를 scikit-learn 으로 따로 구현한다.

R 과 같은 학습 시료·같은 겹(LASSO foldid, rfe 겹)·같은 λ 경로를 받아 와서, 차이가 구현에서만 나오게 한다.
- LASSO: glmnet 목적함수 (1/n)Σlogloss + λ‖β‖₁ (표준화, 모집단 SD) = sklearn C = 1/(nλ). CV 이탈도는 확률을 [1e-5, 1−1e-5] 로 자른다(glmnet 과 같음).
- SVM-RFE: caret caretFuncs 와 같은 순서 — 겹마다 단변량 ROC 면적 max(AUC, 1−AUC) 로 순위 → 크기별 RBF SVM 정확도 → 최고 정확도의 가장 작은 크기.
  kernlab 은 열을 표본 SD 로 표준화하고, 상수 열이 하나라도 있으면 표준화를 통째로 끈다 → 똑같이 한다. kernlab sigma = sklearn gamma.
- RF: 나무 500, mtry = floor(√210) = 14. 난수열이 R 과 달라 상위 10 겹침만 본다.
실행: python src/py/03_ml.py [tag]  (기본 tag = vst_table1_sample_no, 프로젝트 폴더에서)
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

warnings.filterwarnings("ignore", category=ConvergenceWarning)
ROOT = Path(__file__).resolve().parents[2]
tag = sys.argv[1] if len(sys.argv) > 1 else "vst_table1_sample_no"
input_name = tag.split("_table1")[0].split("_gsm")[0]

d = pd.read_csv(ROOT / f"data/processed/ml_input_{input_name}.tsv.gz", sep="\t").set_index("gsm")
split = pd.read_csv(ROOT / f"data/processed/03_split_{tag}.tsv", sep="\t").set_index("gsm")
d = d.loc[split.index]  # R 과 같은 학습 시료·같은 순서
y = (d.pop("condition") == "OT2D").astype(int).values
X, genes = d.values, d.columns.to_numpy()
n = len(y)

# 1) LASSO
lam = pd.read_csv(ROOT / f"results/tables/03_lasso_path_R_{tag}.tsv", sep="\t")["lambda"].values  # R 과 같은 λ 경로


def lasso_path(Xtr, ytr, lams):
    """λ 경로 전체를 푼다. liblinear 는 절편도 벌점을 주므로 intercept_scaling 을 크게 해 벌점을 사실상 없앤다
    (saga 는 같은 해를 내지만 거의 분리되는 자료에서 수십 배 느렸다)."""
    mu, sd = Xtr.mean(0), Xtr.std(0)  # glmnet: 모집단 SD
    sd[sd == 0] = 1
    m = LogisticRegression(penalty="l1", solver="liblinear", intercept_scaling=1e4, tol=1e-8, max_iter=100000)
    fits = []
    for lv in lams:
        m.set_params(C=1 / (len(ytr) * lv)).fit((Xtr - mu) / sd, ytr)
        fits.append((m.coef_[0].copy(), m.intercept_[0]))
    return fits, mu, sd


foldid = split.foldid_lasso.values
dev = np.zeros((len(lam), n))
for f in np.unique(foldid):
    te = foldid == f
    fits, mu, sd = lasso_path(X[~te], y[~te], lam)
    for k, (b, b0) in enumerate(fits):
        p = np.clip(1 / (1 + np.exp(-((X[te] - mu) / sd @ b + b0))), 1e-5, 1 - 1e-5)
        dev[k, te] = -2 * (y[te] * np.log(p) + (1 - y[te]) * np.log(1 - p))
cvm = dev.mean(1)
k_min = np.where(cvm <= cvm.min())[0].min()  # glmnet: 최소 CV 이탈도를 내는 가장 큰 λ
lam_min = lam[k_min]
fits, _, _ = lasso_path(X, y, lam[:k_min + 1])
lasso_genes = genes[fits[-1][0] != 0]

# 2) SVM-RFE (caret rfe)
SIGMA, C = 0.004400041, 0.5


def auc_imp(Xtr, ytr):
    """caret filterVarImp: 열마다 max(AUC, 1−AUC)."""
    r = np.apply_along_axis(rankdata, 0, Xtr)
    n1, n0 = ytr.sum(), (1 - ytr).sum()
    auc = (r[ytr == 1].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
    return np.maximum(auc, 1 - auc)


def svm_predict(Xtr, ytr, Xte):
    if (Xtr.std(0) > 0).all():  # kernlab: 상수 열이 하나라도 있으면 표준화하지 않는다
        mu, sd = Xtr.mean(0), Xtr.std(0, ddof=1)
        Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
    return SVC(kernel="rbf", gamma=SIGMA, C=C).fit(Xtr, ytr).predict(Xte)


sizes = np.arange(1, X.shape[1] + 1)
folds = np.unique(split.fold_rfe)
acc, imp = np.zeros((len(folds), len(sizes))), np.zeros((len(folds), X.shape[1]))
for i, f in enumerate(folds):
    te = split.fold_rfe.values == f
    imp[i] = auc_imp(X[~te], y[~te])
    order = np.argsort(-imp[i], kind="stable")  # 겹 안의 순위: 동점은 열 순서(caret 과 같음)
    for j, s in enumerate(sizes):
        keep = order[:s]
        acc[i, j] = (svm_predict(X[~te][:, keep], y[~te], X[te][:, keep]) == y[te]).mean()
best = sizes[np.argmax(acc.mean(0))]  # 최고 정확도 중 가장 작은 크기
# 최종 유전자: caret pickVars 처럼 겹별 중요도의 평균으로 순위, 동점은 유전자 이름순.
# 학습 20명에서는 AUC = 1 인 유전자가 14개라 이 동점 규칙이 어떤 유전자가 뽑히는지를 정한다.
svm_genes = genes[np.lexsort((genes, -imp.mean(0)))[:best]]

# 3) RF
rf = RandomForestClassifier(n_estimators=500, max_features=int(np.sqrt(X.shape[1])), random_state=100).fit(X, y)
rf_genes = genes[np.argsort(-rf.feature_importances_, kind="stable")[:10]]

hub = [g for g in lasso_genes if g in set(svm_genes) & set(rf_genes)]
out = pd.concat([pd.DataFrame({"method": meth, "rank": np.arange(1, len(g) + 1), "gene": g})
                 for meth, g in (("LASSO", lasso_genes), ("SVM-RFE", svm_genes), ("RF", rf_genes), ("hub", hub))])
out.to_csv(ROOT / f"results/tables/03_ml_py_{tag}.tsv", sep="\t", index=False)
pd.DataFrame({"Variables": sizes, "Accuracy": acc.mean(0)}).to_csv(
    ROOT / f"results/tables/03_rfe_sizes_py_{tag}.tsv", sep="\t", index=False)
pd.DataFrame({"lambda": lam, "cvm": cvm}).to_csv(ROOT / f"results/tables/03_lasso_path_py_{tag}.tsv", sep="\t", index=False)
print(f"[py {tag}] LASSO λmin {lam_min:.8f}, {len(lasso_genes)}개: {', '.join(lasso_genes)}")
print(f"SVM-RFE {best}개: {', '.join(svm_genes)}")
print(f"RF 상위 10: {', '.join(rf_genes)}")
print(f"세 방법 공통: {', '.join(hub)}")
