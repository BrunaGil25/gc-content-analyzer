# tests/test_gc_analysis.py
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from gc_analysis import gc_content, run_analysis


import os
import pandas as pd
from gc_analysis import gc_content, run_analysis

def test_gc_content_basic():
    assert gc_content("GCGCGC") == 100.0
    assert gc_content("ATATAT") == 0.0
    assert gc_content("ATGC") == 50.0
    assert gc_content("") == 0.0

def test_run_analysis_outputs():
    # Use small test files from your data/ folder
    fna = "data/ecoli_genomic.fna"
    gbff = "data/ecoli_genomic.gbff"
    faa = "data/ecoli_proteins.faa"

    results = run_analysis(fna, gbff, faa, results_dir="results_test")

    # Check keys exist
    assert "sequence_gc" in results
    assert "sequence_plot" in results
    assert "gene_gc" in results
    assert "gene_plot" in results
    assert "aa_comp" in results
    assert "aa_plot" in results

    # Check DataFrame structure
    assert isinstance(results["sequence_gc"], pd.DataFrame)
    assert "gc_percent" in results["sequence_gc"].columns

    # Clean up
    if os.path.exists("results_test"):
        for f in os.listdir("results_test"):
            os.remove(os.path.join("results_test", f))
        os.rmdir("results_test")


