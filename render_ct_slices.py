import os, glob
import numpy as np
import pydicom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/root/autodl-tmp/radiomics_ct/DICOM"
OUT = "/root/autodl-tmp/R_figures/CT_slices"
os.makedirs(OUT, exist_ok=True)

# 全部 7 个患者：3 耐药(NY) + 4 药敏(YM)
patients = {
    "0094802": "NY3", "0119476": "NY4", "0112284": "NY5",
    "0119359": "YM2", "0119397": "YM3", "0119181": "YM4", "0119202": "YM5",
}

def read_axial_series(pid):
    pd = [d for d in os.listdir(BASE) if d.startswith(pid)][0]
    pd_path = os.path.join(BASE, pd)
    ct_dir = None
    for root, dirs, files in os.walk(pd_path):
        for d in dirs:
            if d.endswith("_CT"):
                ct_dir = os.path.join(root, d)
                break
        if ct_dir:
            break
    fs = glob.glob(os.path.join(ct_dir, "*.dcm"))
    slices = []
    for f in fs:
        ds = pydicom.dcmread(f)
        it = getattr(ds, "ImageType", [])
        if len(it) > 2 and it[2] == "AXIAL":
            loc = getattr(ds, "SliceLocation", None)
            if loc is None:
                ip = getattr(ds, "ImagePositionPatient", [0,0,0])
                loc = ip[2]
            slices.append((float(loc), ds.pixel_array.astype(np.float32)))
    slices.sort(key=lambda x: x[0])
    return np.stack([s[1] for s in slices])

def lung_window(arr):
    wl, ww = -600.0, 1500.0
    low, high = wl - ww/2, wl + ww/2
    arr = np.clip(arr, low, high)
    arr = (arr - low) / (high - low) * 255.0
    return arr.astype(np.uint8)

for pid, label in patients.items():
    print(f"渲染 {pid} ({label})...")
    arr = read_axial_series(pid)
    print(f"  轴位切片: {arr.shape[0]}")

    lung_mask = (arr > -1000) & (arr < -400)
    lung_area = lung_mask.sum(axis=(1, 2))
    # 选肺面积 top 3 的切片（病灶最明显层面）
    top3 = np.argsort(lung_area)[-3:][::-1]
    print(f"  top3 切片: {top3.tolist()}")

    for rank, sl in enumerate(top3, 1):
        slice_arr = lung_window(arr[sl])
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(slice_arr, cmap="gray", vmin=0, vmax=255)
        ax.axis("off")
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.savefig(os.path.join(OUT, f"CT_{label}_{pid}_{rank}.png"), dpi=300, bbox_inches="tight", pad_inches=0)
        plt.close(fig)

print("=== 完成, 总数:", len(os.listdir(OUT)), "===")
