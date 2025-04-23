# DeepVRegulome

DeepVRegulome is an end‑to‑end framework for predicting the functional impact of small somatic variants in non‑coding regulatory regions (splice sites and transcription‑factor‑binding sites) using fine‑tuned DNABERT models.

## Features
* Task‑specific DNABERT classifiers (acceptor, donor, and ~700 TFBS models)
* Variant‑effect scoring (Δp and log₂ odds ratio)
* Streamlit dashboard for interactive exploration
* Modular codebase for retraining on new genomes or diseases

## Repository structure
```
src/                 Core Python package
  deepvregulome/
    __init__.py
    models.py        # model loading utilities
    scoring.py       # Δp & log‑odds functions
    data.py          # dataset loaders
models/              Pre‑trained checkpoints or download script
dashboard/           Streamlit app
notebooks/           Example Jupyter notebooks
data/                Sample input files (VCF, BED)
docs/                Documentation stubs
requirements.txt     Python dependencies
setup.py             Install script (optional)
```
## Installation
```bash
git clone https://github.com/DavuluriLab//DeepVRegulome.git
cd DeepVRegulome
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

## Model checkpoints
Full DNABERT fine-tuned weights (acceptor, donor, and 700 TFBS models) will be deposited in Zenodo and made publicly available immediately upon journal acceptance.
In the meantime, researchers may request access by emailing pratik.dutta@stonybrook.edu and ramana.davuluri@stonybrookmedicine.edu  with a brief statement of intended use.




## Quick start (command‑line)
```bash
python -m deepvregulome.scoring \
    --vcf example/variants.vcf \
    --region-bed example/splice_sites.bed \
    --checkpoint models/splice_acceptor.pt
```
## Live Demo

An interactive instance of the DeepVRegulome dashboard is hosted here:

➡️ **https://davuluri-lab-brainved.streamlit.app/**

The deployed app lets you browse model performance metrics and variant-effect predictions without installing any software locally.



## Citation
If you use DeepVRegulome in your research, please cite:



## Licence
MIT. See [LICENSE](LICENSE) for details.
