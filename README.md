# DeepVRegulome

DeepVRegulome is an end‑to‑end framework for predicting the functional impact of small somatic variants in non‑coding regulatory regions (splice sites and transcription‑factor‑binding sites) using fine‑tuned DNABERT models.

## Features
* Task‑specific DNABERT classifiers (acceptor, donor, and ~700 TFBS models)
* Variant‑effect scoring (Δp and log₂ odds ratio)
* Streamlit dashboard for interactive exploration
* Modular codebase for retraining on new genomes or diseases

## Repository structure
```
DeepVRegulome/
├── .devcontainer/
├── .streamlit/
├── data/
│   └── Brain/
├── figures/
│   └── attention/
│       ├── CTCFL/
│       └── ZNF384/
├── notebooks/
│   ├── 01_parse_and_merge_vcfs.ipynb
│   ├── 02_tfbs_intersection.ipynb
│   |── 03_dnabert_input_generation.ipynb
│   ├── 04_scoring_candidate_variants.ipynb
│   └── 05_tfbs_attention_motif_visualization.ipynb
├── scripts/
│   ├── run_prediction_tfbs.sh
│   └── run_prediction_splice_acceptor.sh
├── src/
│   └── deepvregulome/
│       ├── __init__.py
│       ├── data_loader.py
│       ├── model.py
│       ├── utils.py
│       └── config.yaml
├── streamlit_app/
│   └── app_variant_clinical_dashboard.py
├── LICENSE
├── README.md
├── requirements.txt
└── .gitignore

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
