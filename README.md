# 🧬 GC Content Analyzer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/BrunaGil25/gc-content-analyzer/blob/main/LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen.svg)](https://pre-commit.com/)
[![Tests](https://img.shields.io/badge/tests-passing-success.svg)](https://github.com/BrunaGil25/gc-content-analyzer/actions)
[![Streamlit](https://img.shields.io/badge/Streamlit-app-orange.svg)](http://localhost:8501/)
[![Docs](https://img.shields.io/badge/docs-MkDocs-lightgrey.svg)](https://github.com/BrunaGil25/gc-content-analyzer/tree/main/docs)


A Python-based toolkit for analyzing GC content and sequence metadata from NCBI files. Supports `.fna`, `.gbff`, and `.faa` formats, with outputs in CSV and visual plots. Includes a Streamlit app for interactive exploration and MkDocs-powered documentation.


---

## Table of Contents

- [Project Structure](#-project-structure)
- [Getting Started](#-getting-started)
  - [Installation](#installation)
  - [Running the App](#running-the-app)
- [Testing](#-testing)
- [Code Quality](#-code-quality)
- [Documentation](#-documentation)
- [Features](#-features)
- [Usage](#usage)
  - [Jupyter Notebook](#1-jupyter-notebook)
  - [Streamlit App](#2-streamlit-app)
- [Ideas for Extension](#ideas-for-extension)
- [License](#license)
- [Author](#author)

## Project Structure

gc-content-analyzer/ 
├── src/                  # Core package code 
│   └── gc_analysis/ 
├── tests/                # Unit tests 
├── data/                 # Input files (.fna, .gbff, .faa) 
├── results/              # Output CSVs and plots 
├── ncbi_temp/            # Temporary files 
├── notebooks/            # Jupyter exploration 
├── docs/                 # MkDocs documentation 
├── app.py                # Streamlit app entry point 
├── pyproject.toml        # Packaging config 
├── .pre-commit-config.yaml 
├── README.m

---

## Getting Started

### Installation

Clone the repo and install dependencies:

#```bash
git clone https://github.com/your-username/gc-content-analyzer.git
cd gc-content-analyzer
pip install -e .

Running the App
Launch the Streamlit interface:
streamlit run app.py

## Testing
Run unit tests with:
pytest

## Code Quality
This project uses pre-commit hooks for formatting and linting:
pip install pre-commit
pre-commit install

Hooks include:
- black for formatting
- flake8 for linting
- isort for import sorting
- check-yaml for YAML validation

##  Documentation
Browse the full user guide and API reference locally:
mkdocs serve
To deploy online: 
mkdocs gh-deploy


## Features

- GC content analysis from .fna, .gbff, and .faa files
- CSV export of results
- Plot generation (GC distribution, sequence length)
- Streamlit app for interactive use
- Jupyter notebooks for exploration
- MkDocs documentation site


---

## Installation

1. Clone the repository  
  
   git clone https://github.com/your-username/gc-content-analyzer.git
   cd gc-content-analyzer

2. (Optional) Create and activate a virtual environment
 
   python -m venv venv
   source venv/bin/activate   # macOS/Linux
   venv\Scripts\activate      # Windows

3. Install dependencies
 
   pip install -r requirements.txt

## Usage

1. Jupyter Notebook
1.1. Place your data files (.fna, .gbff, .faa) in data/.
1.2. Launch the notebook: jupyter notebook notebooks/GC_Content_Analyzer.ipynb
1.3. Run all cells to generate plots and CSV outputs in results/.

2. Streamlit App
2.1. (Optional) Set your NCBI email as a secret or environment variable:
export STREAMLIT_SECRETS='{"NCBI_EMAIL": "your.email@example.com"}'

2.2. Run the app:

    streamlit run app.py

2.3. In the browser:
- Switch between “Upload files” and “Download from NCBI” tabs
- Provide files or an accession (e.g. NC_000913.3, GCF_001592505.1) + your email
- Click Run Analysis and view results instantly

## Ideas for Extension- Compare multiple genomes side-by-side
- Add codon usage bias analysis
- Generate a downloadable PDF report
- Integrate Plotly/Altair for interactive charts

## License MIT 
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Author
Made by Bruna Gil. Data-driven, clean, and powerful.
🔗https://github.com/BrunaGil25 | 🔗 https://www.linkedin.com/in/bruna-gil-garcia-80656069/
