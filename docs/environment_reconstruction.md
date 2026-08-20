# Environment reconstruction record

Do not copy the research-lab Conda environment, virtual environment, site-packages,
CUDA installation or caches. Reconstruct an isolated environment from the project
metadata and verify it with the repository tests.

## Recorded software

| Component | Recorded version |
|---|---:|
| Python | 3.9.18 |
| PyYAML | 6.0.1 |
| NumPy | 1.26.4 |
| SciPy | 1.12.0 |
| h5py | 3.1.0 |
| matplotlib | 3.7.5 |
| R | 4.5.2 |
| SCPA | 1.6.2 |
| multicross | 2.1.0 |
| Matrix | 1.7.6 |
| jsonlite | 2.0.0 |
| rhdf5 | 2.54.1 |
| Seurat | 5.5.1 |
| SeuratObject | 5.4.0 |

Phase 7 is CPU-only and does not require a model runtime, PyTorch, CUDA, or GPU.
CUDA-compatible packages are therefore outside the Phase 7 reconstruction contract.

## Verification

```bash
python --version
python -m pip install -e .
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
Rscript -e 'pkgs <- c("SCPA","multicross","Matrix","jsonlite","rhdf5"); for (p in pkgs) cat(p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT_INSTALLED", "\n")'
```

Do not add a GPU stack merely to run Phase 7; SCPA/multicross runs on CPU.
