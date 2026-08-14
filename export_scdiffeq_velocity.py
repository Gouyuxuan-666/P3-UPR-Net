import scanpy as sc, scdiffeq as sdq, umap, numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")

print("Loading...")
adata = sc.read_h5ad("/root/autodl-tmp/tb_15samples.h5ad")
adata = adata[adata.obs["group"].isin(["HC","NY","YM"])].copy()
print("Cells:", adata.n_obs)
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger")
adata = adata[:, adata.var.highly_variable].copy()
sc.tl.pca(adata, n_comps=50)
adata.obs["t"] = adata.obs["group"].map({"HC":0,"YM":1,"NY":2}).astype(float)
UMAP = umap.UMAP(n_components=2, random_state=42)
adata.obsm["X_umap"] = UMAP.fit_transform(adata.obsm["X_pca"])

print("Fitting scDiffEq (500 epochs CPU)...")
model = sdq.scDiffEq(adata)
model.fit(train_epochs=500)
model.drift()
model.diffusion()
sdq.tl.velocity_graph(model.adata)

print("obsm keys:", list(model.adata.obsm.keys()))
print("uns keys:", list(model.adata.uns.keys()))
print("layers keys:", list(model.adata.layers.keys()))

umap_coords = np.array(model.adata.obsm["X_umap"])
group = np.array(model.adata.obs["group"].values)

vel = None; vel_key = None
for k in ["velocity", "drift", "X_drift", "X_velocity"]:
    if k in model.adata.obsm:
        vel = np.array(model.adata.obsm[k]); vel_key = k; break
if vel is None:
    for k in model.adata.uns.keys():
        if "velocity" in k.lower() or "drift" in k.lower():
            print("uns candidate:", k)

if vel is not None:
    print("velocity key:", vel_key, "shape:", vel.shape)
    df = pd.DataFrame({"umap1":umap_coords[:,0], "umap2":umap_coords[:,1],
                       "vel1":vel[:,0], "vel2":vel[:,1], "group":group})
    df.to_csv("/root/autodl-tmp/scdiffeq_velocity.csv", index=False)
    print("velocity 已导出, 行数:", len(df))
else:
    print("未找到 velocity, 检查 obsm 全键:", {k: model.adata.obsm[k].shape for k in model.adata.obsm.keys()})
