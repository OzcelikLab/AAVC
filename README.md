## AAVC: Automated ACMG-based Variant Classifier
AAVC is a computational framework for the clinical interpretation of genomic variants based on the standards and guidelines of the American College of Medical Genetics and Genomics (ACMG) and ClinGen specifications. By integrating large-scale public databases and _in silico_ prediction tools, AAVC provides a robust, automated solution for variant classification.

### Key Features
- Implements ACMG/AMP guidelines and ClinGen refinements
- Achieves 99.3% concordance with FDA-approved variant classifications
- Resolves ~50% of variants of uncertain significance (≈710,000 of 1.38M in ClinVar)
- Supports population-scale curation (e.g., >300 novel variants identified in the Turkish Variome)
- Highlights clinically actionable findings (1 in 15 individuals carried an actionable genotype)

### Installation

**1. Get prerequisites**

```bash
sudo apt-get update
sudo apt-get install -y git git-lfs
```

**2. Clone the repository**

```bash
git clone https://github.com/OzcelikLab/AAVC.git
cd AAVC
git lfs pull
```

**3. Run the installer**

```bash
sudo ./install.sh
```

### Usage

This script (`aavc.py`) runs the AAVC on a list of variants or a single variant ID. It can also process pre-annotated VCF files.

Remember to activate the virtual environment before running AAVC:

```bash
source aavc_env/bin/activate
```

**1. Process a text file with variant IDs**

`input.txt` should contain one variant per line in the format `chr-pos-ref-alt` (e.g., `12-106992962-T-G`).

```bash
python aavc.py input.txt
```

**2. Process a single variant**

Provide variant in the format `chr-pos-ref-alt` (e.g., `12-106992962-T-G`).

```bash
python aavc.py 12-106992962-T-G
```

**3. Process a pre-calculated VCF file**

The `input.vcf` file should be annotated with the following: variant frequency, number of homozygous individuals, phyloP conservation scores, REVEL or BayesDel pathogenicity scores, and SpliceAI maximum delta scores. Update `vcf_config.txt` to set the annotation options.

```bash
python aavc.py input.vcf --vcf_mode
```

Optionally keep the `INFO` column in the output:

```bash
python aavc.py input.vcf --vcf_mode --keep_info
```

