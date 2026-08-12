import os, subprocess, shutil, glob

os.makedirs("results", exist_ok=True)

for i in range(3):
    subprocess.run(["root", "-l", "-q", f"Standard_N.C(0,{i},0)"], check=True)

os.makedirs("results/pdf", exist_ok=True)
os.makedirs("results/root", exist_ok=True)

for f in glob.glob("results/*.pdf"):
    shutil.move(f, "results/pdf/")
for f in glob.glob("results/*.root"):
    shutil.move(f, "results/root/")

shutil.rmtree("results_LAL", ignore_errors=True)
shutil.move("results", "results_LAL")
