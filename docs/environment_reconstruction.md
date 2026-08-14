# Environment reconstruction record

Do not copy the research-lab Conda environment, virtual environment, site-packages,
CUDA installation or caches. Reconstruct an isolated environment and compare it
with `artifacts/phase7_environment_snapshot.json`.

## Recorded software

| Component | Recorded version |
|---|---:|
| Python | 3.9.18 |
| PyYAML | 6.0.1 |
| NumPy | 1.26.4 |
| SciPy | 1.12.0 |
| h5py | 3.1.0 |
| matplotlib | 3.7.5 |
| PyTorch | 2.8.0 |
| Transformers | 4.57.6 |
| Triton | 3.4.0 |
| R | 4.5.2 |
| SCPA | 1.6.2 |
| multicross | 2.1.0 |
| Matrix | 1.7.6 |
| jsonlite | 2.0.0 |
| rhdf5 | 2.54.1 |
| Seurat | 5.5.1 |
| SeuratObject | 5.4.0 |

`accelerate`, `kernels`, and `openai-harmony` were not installed in the recorded
lab environment. Phase 7 declares `accelerate` and `kernels` as future runtime
requirements; their versions must be frozen only after an approved runtime setup.
The current gate must therefore remain unsupported. No model download is part of
environment reconstruction.

## Verification

```bash
python --version
python -m pip install -e .
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
Rscript -e 'pkgs <- c("SCPA","multicross","Matrix","jsonlite","rhdf5"); for (p in pkgs) cat(p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT_INSTALLED", "\n")'
```

Install the CUDA-compatible PyTorch stack using the destination GPU/driver-specific
reviewed procedure. Do not infer CUDA binaries from a copied environment.
