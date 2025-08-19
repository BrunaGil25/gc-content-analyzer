import os
import sys

import pandas as pd
import streamlit as st
from Bio import Entrez, SeqIO

# Add gc_analysis.py to import path
sys.path.insert(0, os.path.abspath("."))
from gc_analysis import run_analysis

# ───────────────────────────────────────────────────────────────
# Page Setup
st.set_page_config(page_title="GC Content Analyzer", layout="wide")

# ───────────────────────────────────────────────────────────────
# Sidebar Guide
st.sidebar.markdown(
    """
## 🧭 How to Use This App

1. Upload your own genome files **or** enter an NCBI accession
2. Click **Run Analysis** to process your data
3. View GC content plots and tables
4. Download results for further exploration

Supported formats: `.fna`, `.gbff`, `.faa`
"""
)

# ───────────────────────────────────────────────────────────────
# Welcome Message
st.markdown(
    """
# 🧬 GC Content Analyzer

GC content — the proportion of guanine (G) and cytosine (C) bases in DNA — is a key metric in genomics.
It affects gene expression, DNA stability, and evolutionary patterns.
This tool helps you explore GC content and amino acid composition across sequences from your own files or directly from NCBI.
"""
)


# ───────────────────────────────────────────────────────────────
# NCBI Download Function
def download_ncbi(accession: str, email: str):
    Entrez.email = email
    acc = accession.strip()
    prefix = acc.split("_", 1)[0].upper()

    if prefix in ("NC", "NM", "MT"):
        fna_records = Entrez.efetch(
            db="nuccore", id=acc, rettype="fasta", retmode="text"
        ).read()
        gb_records = Entrez.efetch(
            db="nuccore", id=acc, rettype="gbwithparts", retmode="text"
        ).read()
        nuccore_ids = [acc]

    elif prefix in ("GCF", "GCA"):
        search_results = Entrez.read(
            Entrez.esearch(db="assembly", term=acc, retmode="xml")
        )
        uids = search_results.get("IdList", [])
        if not uids:
            raise ValueError(f"Assembly accession '{acc}' not found.")
        asm_uid = uids[0]

        link_results = Entrez.read(
            Entrez.elink(
                dbfrom="assembly", db="nuccore", id=asm_uid, linkname="assembly_nuccore"
            )
        )
        linksets = link_results[0].get("LinkSetDb", [])
        if not linksets or not linksets[0].get("Link"):
            raise ValueError(f"No nuccore records linked to assembly '{acc}'.")
        nuccore_ids = [lnk["Id"] for lnk in linksets[0]["Link"]]

        ids_csv = ",".join(nuccore_ids)
        fna_records = Entrez.efetch(
            db="nuccore", id=ids_csv, rettype="fasta", retmode="text"
        ).read()
        gb_records = Entrez.efetch(
            db="nuccore", id=ids_csv, rettype="gbwithparts", retmode="text"
        ).read()

    else:
        raise ValueError("Invalid accession prefix. Use NC_/NM_/MT_ or GCF_/GCA_.")

    with open("ncbi_temp.fna", "w") as f:
        f.write(fna_records)
    with open("ncbi_temp.gbff", "w") as f:
        f.write(gb_records)

    proteins = []
    for rec in SeqIO.parse("ncbi_temp.gbff", "genbank"):
        for feat in rec.features:
            if feat.type == "CDS" and "translation" in feat.qualifiers:
                locus = feat.qualifiers.get("locus_tag", ["unknown"])[0]
                aa_seq = feat.qualifiers["translation"][0]
                proteins.append(f">{locus}\n{aa_seq}\n")

    with open("ncbi_temp.faa", "w") as f:
        f.writelines(proteins)

    return "ncbi_temp.fna", "ncbi_temp.gbff", "ncbi_temp.faa"


# ───────────────────────────────────────────────────────────────
# Tabs: Upload vs NCBI
tab1, tab2 = st.tabs(["🔼 Upload Files", "🌐 Fetch from NCBI"])

# ──────────────── Tab 1: Upload ────────────────
with tab1:
    st.header("Upload Your Genome Files")
    fna_file = st.file_uploader("Genome FASTA (.fna)", type="fna")
    gbff_file = st.file_uploader("GenBank (.gbff)", type="gbff")
    faa_file = st.file_uploader("Protein FASTA (.faa)", type="faa")

    if st.button("Run Analysis (Upload)"):
        if not fna_file:
            st.warning("Please upload at least a .fna file.")
        else:
            with open("temp.fna", "wb") as f:
                f.write(fna_file.getvalue())
            gb_path = faa_path = None
            if gbff_file:
                with open("temp.gbff", "wb") as f:
                    f.write(gbff_file.getvalue())
                gb_path = "temp.gbff"
            if faa_file:
                with open("temp.faa", "wb") as f:
                    f.write(faa_file.getvalue())
                faa_path = "temp.faa"

            results = run_analysis("temp.fna", gb_path, faa_path)

            st.markdown("### 📈 GC Content per Sequence")
            st.dataframe(results["sequence_gc"], use_container_width=True)
            st.download_button(
                "📥 Download GC Content per Sequence",
                results["sequence_gc"].to_csv(index=False).encode("utf-8"),
                "gc_content_per_sequence.csv",
                "text/csv",
            )
            st.pyplot(results["sequence_plot"])

            if "gene_gc" in results:
                st.markdown("### 🧬 GC Content per Gene")
                st.dataframe(results["gene_gc"], use_container_width=True)
                st.download_button(
                    "📥 Download GC Content per Gene",
                    results["gene_gc"].to_csv(index=False).encode("utf-8"),
                    "gc_content_per_gene.csv",
                    "text/csv",
                )
                st.pyplot(results["gene_plot"])

            if "aa_comp" in results:
                st.markdown("### 🔠 Amino Acid Composition")
                st.dataframe(results["aa_comp"], use_container_width=True)
                st.download_button(
                    "📥 Download Amino Acid Composition",
                    results["aa_comp"].to_csv(index=False).encode("utf-8"),
                    "amino_acid_composition.csv",
                    "text/csv",
                )
                st.pyplot(results["aa_plot"])

# ──────────────── Tab 2: NCBI ────────────────
with tab2:
    st.header("Download Genome from NCBI")
    acc = st.text_input(
        "NCBI Accession", placeholder="e.g. NC_000913.3 or GCF_001592505.1"
    )
    email = st.text_input("Your Email", help="Required by NCBI for data access.")

    if st.button("Download & Run Analysis"):
        if not acc or not email:
            st.warning("Both accession and email are required.")
        else:
            with st.spinner("Fetching data from NCBI..."):
                try:
                    fna_path, gbff_path, faa_path = download_ncbi(acc, email)
                    results = run_analysis(fna_path, gbff_path, faa_path)
                except Exception as e:
                    st.error(f"Error: {e}")
                    st.stop()

            st.success(f"Results for accession '{acc}':")

            st.markdown("### 📈 GC Content per Sequence")
            st.dataframe(results["sequence_gc"], use_container_width=True)
            st.download_button(
                "📥 Download GC Content per Sequence",
                results["sequence_gc"].to_csv(index=False).encode("utf-8"),
                "gc_content_per_sequence.csv",
                "text/csv",
            )
            st.pyplot(results["sequence_plot"])

            if "gene_gc" in results:
                st.markdown("### 🧬 GC Content per Gene")
                st.dataframe(results["gene_gc"], use_container_width=True)
                st.download_button(
                    "📥 Download GC Content per Gene",
                    results["gene_gc"].to_csv(index=False).encode("utf-8"),
                    "gc_content_per_gene.csv",
                    "text/csv",
                )
                st.pyplot(results["gene_plot"])

            if "aa_comp" in results:
                st.markdown("### 🔠 Amino Acid Composition")
                st.dataframe(results["aa_comp"], use_container_width=True)
                st.download_button(
                    "📥 Download Amino Acid Composition",
                    results["aa_comp"].to_csv(index=False).encode("utf-8"),
                    "amino_acid_composition.csv",
                    "text/csv",
                )
                st.pyplot(results["aa_plot"])

# ───────────────────────────────────────────────────────────────
# Footer
st.markdown("---")
st.markdown("Made with ❤️ by Bruna Gil — Data-driven, clean, and powerful.")
st.markdown(
    "[🔗 GitHub](https://github.com/BrunaGil25) | [🔗 LinkedIn](https://www.linkedin.com/in/bruna-gil-garcia-80656069/)"
)
