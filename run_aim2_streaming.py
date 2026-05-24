#!/usr/bin/env python3
"""Stream-extract each Aim 2 model from the zip, generate its lookup CSV via
Rscript build_lookup_aim2.R, then delete the .RData to keep peak disk low."""

import os, shutil, subprocess, sys, time, zipfile

ZIP_PATH   = "../all aim2.zip"
MODEL_DIR  = "aim2_models"
SKIP_EXISTING = True  # skip if CSV already exists

CANCERS = ["breast","prostate","colon","rectum","urine","esophagu",
           "melanoma","liver","kidney","ovary","retroper","testis","lnsc","lsc"]
SITES   = ["bone","brain","liver","lung"]

# Locate each model inside the zip by trying known path patterns.
def find_entry(zf, cancer, site):
    fname = f"prostate_decision_{site}.RData" if cancer == "prostate" else f"{cancer}_{site}.RData"
    candidates = [n for n in zf.namelist() if n.endswith("/" + fname) or n == fname]
    if not candidates:
        return None
    return candidates[0]

def extract_one(zf, entry, dst):
    with zf.open(entry) as src_f, open(dst, "wb") as dst_f:
        shutil.copyfileobj(src_f, dst_f, length=64*1024*1024)

def main():
    os.makedirs(MODEL_DIR, exist_ok=True)
    os.makedirs("lookup_tables_aim2", exist_ok=True)

    pending = [(c, s) for c in CANCERS for s in SITES]
    total_t0 = time.time()
    done = 0
    failed = []

    with zipfile.ZipFile(ZIP_PATH) as zf:
        for cancer, site in pending:
            tag = f"{cancer}_{site}"
            csv_path = f"lookup_tables_aim2/{tag}_lookup.csv.gz"
            if SKIP_EXISTING and os.path.exists(csv_path):
                print(f"[{tag}] SKIP — CSV already exists ({os.path.getsize(csv_path)/1e6:.2f} MB)")
                done += 1
                continue

            entry = find_entry(zf, cancer, site)
            if not entry:
                print(f"[{tag}] MISSING in zip", flush=True)
                failed.append(tag)
                continue

            rdata_path = os.path.join(MODEL_DIR, os.path.basename(entry))
            print(f"[{tag}] extracting {entry} ({zf.getinfo(entry).file_size/1e6:.0f} MB)", flush=True)
            t0 = time.time()
            extract_one(zf, entry, rdata_path)
            t_ext = time.time() - t0

            t0 = time.time()
            res = subprocess.run(
                ["Rscript", "build_lookup_aim2.R", tag],
                capture_output=True, text=True
            )
            t_pred = time.time() - t0

            try:
                os.remove(rdata_path)
            except OSError:
                pass

            if res.returncode != 0:
                print(f"[{tag}] FAILED (extract {t_ext:.1f}s, predict {t_pred:.1f}s)")
                print(res.stdout[-2000:])
                print(res.stderr[-2000:], file=sys.stderr)
                failed.append(tag)
                continue

            sz = os.path.getsize(csv_path)/1e6 if os.path.exists(csv_path) else 0
            print(f"[{tag}] done — extract {t_ext:.1f}s, predict+write {t_pred:.1f}s, csv {sz:.2f} MB",
                  flush=True)
            done += 1

    total_t = time.time() - total_t0
    print(f"\nFinished: {done}/{len(pending)} models, total {total_t/60:.1f} min")
    if failed:
        print("FAILED:", failed)
        sys.exit(1)

if __name__ == "__main__":
    main()
