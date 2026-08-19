import glob
import os
import shutil
import subprocess

# Count the number of .root files
n_files = len(glob.glob("rootfiles/*.root"))

os.makedirs("results", exist_ok=True)

# Loop over all ROOT files
for a in range(n_files):

    # Multiplicity bins b = 3 to 17
    # kT bin c = 0
    for b in range(3, 17):
        subprocess.run(
            ["root", "-l", "-q", "-b", f"Standard_N.C({a},{b},0)"],
            check=True
        )

    # Multiplicity b = 1
    # kT bins c = 1 to 5
    for c in range(1, 6):
        subprocess.run(
            ["root", "-l", "-q", "-b", f"Standard_N.C({a},1,{c})"],
            check=True
        )

    # Multiplicity b = 2
    # kT bins c = 1 to 5
    for c in range(1, 6):
        subprocess.run(
            ["root", "-l", "-q", "-b", f"Standard_N.C({a},2,{c})"],
            check=True
        )

# Create output directories
os.makedirs("results/pdf", exist_ok=True)
os.makedirs("results/root", exist_ok=True)

# Move PDFs
for f in glob.glob("results/*.pdf"):
    shutil.move(f, "results/pdf/")

# Move ROOT files
for f in glob.glob("results/*.root"):
    shutil.move(f, "results/root/")

# Remove old results_KK and rename
shutil.rmtree("results_KK", ignore_errors=True)
shutil.move("results", "results_KK")