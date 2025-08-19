# gc_analysis.py
# 🧬 GC Content Analyzer – modular version
# Author: Bruna Gil

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO

sns.set(style="whitegrid")


def gc_content(seq):
    """Compute GC% for a SeqRecord or raw sequence."""
    seq = str(seq).upper()
    return (
        round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)
        if len(seq) > 0
        else 0
    )


def run_analysis(fna_path, gbff_path=None, faa_path=None, results_dir="results"):
    """
    Perform three analyses:
      1) GC% per FASTA sequence
      2) GC% distribution across CDS from GenBank
      3) Amino acid composition from protein FASTA

    Saves CSVs to results_dir and returns a dict:
      {
        "sequence_gc": DataFrame,
        "sequence_plot": matplotlib.Figure,
        "gene_gc": DataFrame,
        "gene_plot": Figure,
        "aa_comp": DataFrame,
        "aa_plot": Figure
      }
    """
    os.makedirs(results_dir, exist_ok=True)
    results = {}

    # --- 1) GC% per sequence (FASTA) ---
    seq_records = list(SeqIO.parse(fna_path, "fasta"))
    seq_data = [
        {"id": rec.id, "gc_percent": gc_content(rec.seq)} for rec in seq_records
    ]
    df_seq = pd.DataFrame(seq_data)
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    sns.barplot(
        data=df_seq, x="id", y="gc_percent", hue="id", palette="viridis", ax=ax1
    )
    ax1.tick_params(axis="x", labelrotation=90)
    ax1.set_ylabel("GC Content (%)")
    ax1.set_xlabel("Sequence ID")
    ax1.set_title("GC Content per Sequence")
    fig1.tight_layout()
    df_seq.to_csv(os.path.join(results_dir, "gc_content_per_sequence.csv"), index=False)
    results["sequence_gc"] = df_seq
    results["sequence_plot"] = fig1

    # --- 2) GC% by gene (GenBank) ---
    if gbff_path:
        gene_entries = []
        for rec in SeqIO.parse(gbff_path, "genbank"):
            for feat in rec.features:
                if feat.type == "CDS":
                    seq = feat.extract(rec.seq)
                    gene_entries.append(
                        {
                            "locus_tag": feat.qualifiers.get("locus_tag", ["N/A"])[0],
                            "gene": feat.qualifiers.get("gene", ["N/A"])[0],
                            "gc_percent": gc_content(seq),
                            "length_bp": len(seq),
                        }
                    )
        df_genes = pd.DataFrame(gene_entries)
        fig2, ax2 = plt.subplots(figsize=(12, 6))
        sns.histplot(df_genes["gc_percent"], bins=30, color="teal", ax=ax2)
        ax2.set_xlabel("GC Content (%)")
        ax2.set_ylabel("Number of Genes")
        ax2.set_title("Distribution of GC Content in Genes")
        fig2.tight_layout()
        df_genes.to_csv(
            os.path.join(results_dir, "gc_content_per_gene.csv"), index=False
        )
        results["gene_gc"] = df_genes
        results["gene_plot"] = fig2

    # --- 3) Amino acid composition (protein FASTA) ---
    if faa_path:
        aa_counts = {}
        for rec in SeqIO.parse(faa_path, "fasta"):
            for aa in str(rec.seq):
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
        df_aa = pd.DataFrame.from_dict(aa_counts, orient="index", columns=["Count"])
        df_aa.index.name = "AminoAcid"
        df_aa = df_aa.reset_index()
        df_aa["Frequency"] = df_aa["Count"] / df_aa["Count"].sum()
        df_aa.sort_values("Frequency", ascending=False, inplace=True)

        fig3, ax3 = plt.subplots(figsize=(10, 5))
        sns.barplot(
            x="AminoAcid",
            y="Frequency",
            data=df_aa,
            hue="AminoAcid",
            palette="mako",
            ax=ax3,
        )
        ax3.set_ylabel("Relative Frequency")
        ax3.set_xlabel("Amino Acid")
        ax3.set_title("Amino Acid Composition")
        fig3.tight_layout()
        df_aa.to_csv(
            os.path.join(results_dir, "amino_acid_composition.csv"), index=False
        )
        results["aa_comp"] = df_aa
        results["aa_plot"] = fig3

    return results
