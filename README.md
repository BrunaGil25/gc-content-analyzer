# 🧬 GC Content Analyzer

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/BrunaGil25/gc-content-analyzer/blob/main/LICENSE)  
[![Streamlit App](https://img.shields.io/badge/Streamlit-Live_App-orange.svg)](https://gc-content-analyzer-brunagil.streamlit.app/)  
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/downloads/)  
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen.svg)](https://pre-commit.com/)  
[![Tests](https://img.shields.io/badge/tests-passing-success.svg)](https://github.com/BrunaGil25/gc-content-analyzer/actions)  
[![Docs](https://img.shields.io/badge/docs-MkDocs-lightgrey.svg)](https://github.com/BrunaGil25/gc-content-analyzer/tree/main/docs)

A Python toolkit for exploring GC content and amino acid composition from genomic data.  
Supports `.fna`, `.gbff`, and `.faa` formats. Offers CSV export, publication-quality plots, and an interactive Streamlit interface.

🔗 **Live demo:** https://gc-content-analyzer-brunagil.streamlit.app/

---

## Table of Contents

- [Features](#features)  
- [Project Structure](#project-structure)  
- [Installation](#installation)  
- [Usage](#usage)  
  - [Streamlit App](#streamlit-app)  
  - [Jupyter Notebook](#jupyter-notebook)  
- [Contributing](#contributing)
- [Ideas for Extension](#Ideas-for-extension)
- [License](#license)  
- [Author](#author)  

---

## Features

- Compute GC content per sequence and per gene  
- Analyze amino acid composition from protein sequences  
- Export results as CSV and generate plots  
- Fetch genomes directly from NCBI by accession  
- Interactive Streamlit app for instant visualization  
- MkDocs documentation and Jupyter notebooks for deeper exploration  

---

## Project Structure
 ```text
gc-content-analyzer/ 
├── src/
│   └── gc_analysis/          # Core analysis functions
├── tests/                    # Unit tests
├── data/                     # Example inputs (.fna, .gbff, .faa)
├── results/                  # Generated CSVs & plots
├── ncbi_temp/                # Temporary NCBI downloads
├── notebooks/                # Jupyter exploration
├── docs/                     # MkDocs site
├── app.py                    # Streamlit entry point
├── pyproject.toml            # Packaging & metadata
├── requirements.txt          # Python dependencies
├── .pre-commit-config.yaml   # Git hooks & linters
└── README.md                 # This file

---

## Getting Started

### Installation

Clone the repo and install dependencies:

1. Clone the repository
 ```bash  
git clone https://github.com/BrunaGil25/gc-content-analyzer.git
cd gc-content-analyzer
 ```
2. (Optional) Create and activate a virtual environment
 ```bash
   python -m venv venv
   source venv/bin/activate   # macOS/Linux
   venv\Scripts\activate      # Windows
 ```
4. Install dependencies
    ```bash  
   pip install -r requirements.txt
 ```
5. Install pre-commit hooks
 ```bash
pip install pre-commit
  pre-commit install
 ```
## Usage
### Streamlit App
Launch locally:
   streamlit run app.py

Or visit the hosted version:
https://gc-content-analyzer-brunagil.streamlit.app/
In the browser:
- Switch between “Upload Files” and “Fetch from NCBI” tabs
- Provide your .fna/.gbff/.faa or an NCBI accession + email
- Click Run Analysis
- View tables & plots, then download CSVs with the download buttons

### Jupyter Notebook
- Place your data files in data/.
- Launch the notebook:
   jupyter notebook notebooks/GC_Content_Analyzer.ipynb
- Run all cells to generate results in results/.

## Contributing
We welcome improvements!
- Fork the repo
- Create a branch:
   git checkout -b feature/your-feature
- Make your changes, then run tests:
   pytest
- Commit with a clear message:
   git commit -m "Add feature: description"
- Push and open a Pull Request.

Please follow the existing code style and ensure all linting/tests pass via pre-commit hooks.

## Ideas for Extension
- Compare multiple genomes side-by-side
- Add codon usage bias analysis
- Generate a downloadable PDF report
- Integrate Plotly/Altair for interactive charts

## License MIT 
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Author
Made by Bruna Gil. Data-driven, clean, and powerful.
🔗https://github.com/BrunaGil25 | 🔗 https://www.linkedin.com/in/bruna-gil-garcia-80656069/
