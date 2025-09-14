## AAVC: Automated ACMG-based Variant Classifier
AAVC is a computational framework for the clinical interpretation of genomic variants based on the standards and guidelines of the American College of Medical Genetics and Genomics (ACMG) and ClinGen specifications. By integrating large-scale public databases and in silico prediction tools, AAVC provides a robust, automated solution for variant classification.

#### Key Features
- Implements ACMG/AMP guidelines and ClinGen refinements
- Achieves 99.3% concordance with FDA-approved variant classifications
- Resolves ~50% of variants of uncertain significance (≈710,000 of 1.38M in ClinVar)
- Supports population-scale curation (e.g., >300 novel variants identified in the Turkish Variome)
- Highlights clinically actionable findings (1 in 15 individuals carried an actionable genotype)

### Installation

**1 Prerequisites

- Operating System: Ubuntu/Debian
- Account with sudo privileges.
- Internet connection.

**2 Clone the repository

```bash
git clone https://github.com/ardainn/AAVC.git
cd AAVC
ls -la
```

**3 Run the installer

```bash
./install.sh
```
