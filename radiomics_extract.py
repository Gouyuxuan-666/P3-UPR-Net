import os, glob, csv
import numpy as np
import SimpleITK as sitk
from radiomics import featureextractor

BASE = "/root/autodl-tmp/radiomics_ct/DICOM"
OUT = "/root/autodl-tmp/radiomics_features.csv"

# 样本映射
mapping = {
    "0094802": ("NY3", "NY"),
    "0112284": ("NY5", "NY"),
    "0119476": ("NY4", "NY"),
    "0119359": ("YM2", "YM"),
    "0119397": ("YM3", "YM"),
    "0119181": ("YM4", "YM"),
    "0119202": ("YM5", "YM"),
}

def read_series(patient_dir):
    """读 DICOM 序列为 SimpleITK Image"""
    reader = sitk.ImageSeriesReader()
    dcm_dir = None
    for root, dirs, files in os.walk(patient_dir):
        for d in dirs:
            if d.endswith("_CT"):
                dcm_dir = os.path.join(root, d)
                break
        if dcm_dir: break
    if not dcm_dir: return None
    dicom_names = reader.GetGDCMSeriesFileNames(dcm_dir)
    reader.SetFileNames(dicom_names)
    return reader.Execute()

def segment_lung(image):
    """肺分割：HU 阈值 + 保留最大 2 连通域（左右肺）"""
    arr = sitk.GetArrayFromImage(image)
    # 肺实质 HU 范围
    lung = ((arr > -1000) & (arr < -400)).astype(np.uint8)
    # 连通域，保留最大 2 个
    cc = sitk.GetImageFromArray(lung)
    cc = sitk.ConnectedComponent(cc)
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(cc)
    labels = [l for l in stats.GetLabels() if stats.GetPhysicalSize(l) > 0]
    labels.sort(key=lambda l: stats.GetPhysicalSize(l), reverse=True)
    keep = labels[:2]  # 最大 2 个 = 左右肺
    mask = np.zeros_like(lung)
    lab_arr = sitk.GetArrayFromImage(cc)
    for l in keep:
        mask[lab_arr == l] = 1
    return mask

def main():
    extractor = featureextractor.RadiomicsFeatureExtractor(binWidth=25)
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("firstorder")

    rows = []
    all_feature_names = None
    for pid, (sample, group) in mapping.items():
        patient_dir = os.path.join(BASE, [d for d in os.listdir(BASE) if d.startswith(pid)][0])
        image = read_series(patient_dir)
        if image is None:
            print(f"{pid} ({sample}): 无 CT, 跳过"); continue
        mask = segment_lung(image)
        mask_img = sitk.GetImageFromArray(mask)
        mask_img.CopyInformation(image)

        # PyRadiomics 提取
        try:
            result = extractor.execute(image, mask_img, label=1)
        except Exception as e:
            print(f"{pid} ({sample}): 提取失败 {e}"); continue

        row = {"sample": sample, "group": group, "patient_id": pid}
        row.update(result)
        rows.append(row)
        if all_feature_names is None:
            all_feature_names = list(result.keys())
        n = sum(1 for k in result if k.startswith("original_") or k.startswith("log-") or k.startswith("wavelet-"))
        print(f"{pid} ({sample}, {group}): {n} 特征, 肺体积 {mask.sum()} 体素")

    if rows:
        keys = ["sample", "group", "patient_id"] + [k for k in all_feature_names if k.startswith("original_")]
        with open(OUT, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in keys})
        print(f"\n特征已保存: {OUT}, {len(rows)} 例")
    else:
        print("无数据")

if __name__ == "__main__":
    main()
