# app.py

import sys
import os
from Bio import Entrez, SeqIO
import streamlit as st

# ensure gc_analysis.py is in the import path
sys.path.insert(0, os.path.abspath("."))
from gc_analysis import run_analysis

st.set_page_config(layout="wide")
st.title("Genome GC Content Analyzer 🧬")

def download_ncbi(accession: str, email: str):
    """
    Given an NCBI accession (nucleotide or assembly) and user email,
    fetch genome FASTA, GenBank (with translations), extract proteins,
    and return the three temp file paths.
    """
    Entrez.email = email

    # decide which fetch flow to use
    acc = accession.strip()
    prefix = acc.split("_", 1)[0].upper()

    if prefix in ("NC", "NM", "MT"):
        # --- direct nuccore fetch ---
        fna_records = Entrez.efetch(
            db="nuccore", id=acc, rettype="fasta", retmode="text"
        ).read()
        gb_records = Entrez.efetch(
            db="nuccore", id=acc, rettype="gbwithparts", retmode="text"
        ).read()
        nuccore_ids = [acc]

    elif prefix in ("GCF", "GCA"):
        # --- assembly → nuccore workflow ---
        # 1. find assembly UID
        search_handle = Entrez.esearch(
            db="assembly", term=acc, retmode="xml"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        uids = search_results.get("IdList", [])
        if not uids:
            raise ValueError(f"Assembly accession '{acc}' not found in NCBI.")
        asm_uid = uids[0]

        # 2. link assembly → nuccore IDs
        link_handle = Entrez.elink(
            dbfrom="assembly",
            db="nuccore",
            id=asm_uid,
            linkname="assembly_nuccore"
        )
        link_results = Entrez.read(link_handle)
        link_handle.close()

        linksets = link_results[0].get("LinkSetDb", [])
        if not linksets or not linksets[0].get("Link"):
            raise ValueError(f"No nuccore records linked to assembly '{acc}'.")
        nuccore_ids = [lnk["Id"] for lnk in linksets[0]["Link"]]

        # 3. fetch combined FASTA & GenBank
        ids_csv = ",".join(nuccore_ids)
        fna_records = Entrez.efetch(
            db="nuccore", id=ids_csv, rettype="fasta", retmode="text"
        ).read()
        gb_records = Entrez.efetch(
            db="nuccore", id=ids_csv, rettype="gbwithparts", retmode="text"
        ).read()

    else:
        raise ValueError(
            "Invalid accession prefix. "
            "Must begin with NC_/NM_/MT_ (nucleotide) or GCF_/GCA_ (assembly)."
        )

    # write out fasta (.fna)
    fna_path = "ncbi_temp.fna"
    with open(fna_path, "w") as f:
        f.write(fna_records)

    # write out GenBank (.gbff)
    gbff_path = "ncbi_temp.gbff"
    with open(gbff_path, "w") as f:
        f.write(gb_records)

    # parse proteins from GenBank translations
    proteins = []
    for rec in SeqIO.parse(gbff_path, "genbank"):
        for feat in rec.features:
            if feat.type == "CDS" and "translation" in feat.qualifiers:
                locus = feat.qualifiers.get("locus_tag", ["unknown"])[0]
                aa_seq = feat.qualifiers["translation"][0]
                proteins.append(f">{locus}\n{aa_seq}\n")

    faa_path = "ncbi_temp.faa"
    with open(faa_path, "w") as f:
        f.writelines(proteins)

    return fna_path, gbff_path, faa_path


# UI tabs
tab1, tab2 = st.tabs(["🔼 Upload files", "🌐 Download from NCBI"])

with tab1:
    st.header("Upload your own files")
    fna_file = st.file_uploader("Genome FASTA (.fna)", type="fna")
    gbff_file = st.file_uploader("GenBank (.gbff)", type="gbff")
    faa_file = st.file_uploader("Protein FASTA (.faa)", type="faa")

    if st.button("Run Analysis (Upload)"):
        if not fna_file:
            st.warning("Please upload at least a .fna file.")
        else:
            # save uploads to temp files
            with open("temp.fna", "wb") as f:
                f.write(fna_file.getvalue())
            gb_path = None
            faa_path = None
            if gbff_file:
                with open("temp.gbff", "wb") as f:
                    f.write(gbff_file.getvalue())
                gb_path = "temp.gbff"
            if faa_file:
                with open("temp.faa", "wb") as f:
                    f.write(faa_file.getvalue())
                faa_path = "temp.faa"

            results = run_analysis("temp.fna", gb_path, faa_path)
            # display results
            st.subheader("GC Content per Sequence")
            st.dataframe(results["sequence_gc"], use_container_width=True)
            st.pyplot(results["sequence_plot"])

            if "gene_gc" in results:
                st.subheader("GC Content per Gene")
                st.dataframe(results["gene_gc"], use_container_width=True)
                st.pyplot(results["gene_plot"])

            if "aa_comp" in results:
                st.subheader("Amino Acid Composition")
                st.dataframe(results["aa_comp"], use_container_width=True)
                st.pyplot(results["aa_plot"])

with tab2:
    st.header("Fetch directly from NCBI")
    acc = st.text_input(
        "NCBI accession",
        placeholder="e.g. NC_000913.3 or GCF_001592505.1"
    )
    email = st.text_input(
        "Your email (required by NCBI)",
        help="NCBI policy requires an email address with each request."
    )

    if st.button("Download & Run Analysis"):
        if not acc or not email:
            st.warning("Both accession and email are required.")
        else:
            with st.spinner("Downloading and parsing from NCBI..."):
                try:
                    fna_path, gbff_path, faa_path = download_ncbi(acc, email)
                    results = run_analysis(fna_path, gbff_path, faa_path)
                except Exception as e:
                    st.error(f"Error fetching data: {e}")
                    st.stop()

            st.success(f"Results for accession '{acc}':")
            st.subheader("GC Content per Sequence")
            st.dataframe(results["sequence_gc"], use_container_width=True)
            st.pyplot(results["sequence_plot"])

            if "gene_gc" in results:
                st.subheader("GC Content per Gene")
                st.dataframe(results["gene_gc"], use_container_width=True)
                st.pyplot(results["gene_plot"])

            if "aa_comp" in results:
                st.subheader("Amino Acid Composition")
                st.dataframe(results["aa_comp"], use_container_width=True)
                st.pyplot(results["aa_plot"])

# ───────────────────────────────────────────────────────────────
# Footer: author credit + links
st.markdown("---")
st.markdown("Made by Bruna Gil. Data-driven, clean, and powerful.")
st.markdown(
    "[🔗 GitHub](https://github.com/BrunaGil25) | "
    "[🔗 LinkedIn](https://www.linkedin.com/in/bruna-gil-garcia-80656069/)"
)
