# GPU runtime environments

MetaG scores reactions with UMA through the official FairChem SDK, pinned to commit
`645af7d03e91138ef07e51de387316c9e2e152dc`, with Torch 2.13. Install Torch from the wheel index that
matches your hardware first, then install MetaG with the matching constraints file.

| Hardware | Torch wheel index | Constraints | Runtime profile |
|---|---|---|---|
| V100 (sm_70) | `https://download.pytorch.org/whl/cu126` | `qm-official-cu126.txt` | default; Nitro is unavailable before sm_80, so FP32 is used throughout |
| H100 / recent GPUs, CUDA 13 driver | `https://download.pytorch.org/whl/cu130` | `qm-official-cu130.txt` | `nitro` (default) |

```bash
python -m pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu126
python -m pip install -c environments/qm-official-cu126.txt -e '.[qm]'
export XTB_BIN=/absolute/path/to/xtb          # xtb >= 6.7, single-threaded per process
```

The CUDA 13 Torch build has no `sm_70` kernels, so V100 needs the CUDA 12.6 environment.

The profiles differ only in where UMA Nitro may be used. `nitro` uses it for supported single-system
optimisation steps. Every energy that enters a Boltzmann sum, every convergence check and every Hessian
uses official FP32 in both profiles. The profile, checkpoint hash and SDK identity are recorded in each
result's `config.producer`. Species caches from different producers are never shared.

The UMA weights are downloaded by FairChem from Hugging Face (`facebook/UMA`). This requires accepting
the FAIR Chemistry License and authenticating, for example with `huggingface-cli login`.
