"""00_env_setup.py — verify and install dependencies."""
import subprocess,sys
for p in ["httpx","anndata","scanpy","pyarrow","pandas","numpy","scipy","matplotlib","seaborn"]:
    subprocess.check_call([sys.executable,"-m","pip","install",p,"-q"])
print("Environment ready.")
