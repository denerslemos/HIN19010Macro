import glob
import os
import shutil
import subprocess

# 1. Count the number of .root files in the 'rootfiles' folder
n_files = len(glob.glob("rootfiles/*.root"))

os.makedirs("results", exist_ok=True)

# 2. Outer loop for the first argument (0 to n_files - 1)
for file_idx in range(n_files):
    # 3. Inner loop for the second argument (0 to 2)
    for i in range(3):
        subprocess.run(
            ["root", "-l", "-q", f"Standard_N.C({file_idx},{i},0)"], check=True
        )

os.makedirs("results/pdf", exist_ok=True)
os.makedirs("results/root", exist_ok=True)

for f in glob.glob("results/*.pdf"):
    shutil.move(f, "results/pdf/")
for f in glob.glob("results/*.root"):
    shutil.move(f, "results/root/")

shutil.rmtree("results_LAL", ignore_errors=True)
shutil.move("results", "results_LAL")