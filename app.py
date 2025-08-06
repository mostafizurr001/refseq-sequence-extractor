import os
import re
import io
import zipfile
import smtplib
import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
from email.message import EmailMessage

# ---------------------------
# Setup
# ---------------------------
Entrez.email = os.getenv("NCBI_EMAIL", "your_email@example.com")
# Entrez.api_key = os.getenv("NCBI_API_KEY")

st.set_page_config(page_title="RefSeq Sequence Extractor with Output Format", layout="centered")

st.title("🧬 RefSeq Sequence Extractor (output formatted like example)")
st.markdown("Fetch sequences, compute nucleotide statistics including Z-DNA propensity, and export in the target CSV/Excel schema.")

# Inputs
gene_input = st.text_area("Enter Gene Symbols or RefSeq Accessions (comma/newline separated)", height=150)
organism = st.text_input("Organism (scientific name, e.g., Homo sapiens)", "Homo sapiens")
sequence_type = st.selectbox("Sequence Type", ["Nucleotide", "Protein"])
output_format = st.selectbox("Output Format", ["FASTA", "GenBank"])
send_email = st.checkbox("Email results (requires SMTP env vars)")
user_email = st.text_input("Recipient email") if send_email else ""

# Helpers
def is_refseq_accession(token: str) -> bool:
    return bool(re.match(r"^[A-Z]{2,3}_\\d+(\\.\\d+)?$", token))

def compute_base_stats(seq: str):
    seq = seq.upper()
    A = seq.count("A")
    T = seq.count("T")
    G = seq.count("G")
    C = seq.count("C")
    total = len(seq)
    at = A + T
    gc = G + C
    def pct(x): return (x / total * 100) if total > 0 else 0.0
    return {
        "Length": total,
        "A Count": A,
        "T Count": T,
        "G Count": G,
        "C Count": C,
        "A %": round(pct(A), 3),
        "T %": round(pct(T), 3),
        "G %": round(pct(G), 3),
        "C %": round(pct(C), 3),
        "GC %": round(pct(gc), 3),
        "AT %": round(pct(at), 3),
    }

def zdna_propensity(seq: str):
    s = seq.upper()
    n = len(s)
    covered = [False] * n
    patterns = ["CG", "GC"]
    # check alternating patterns of length 6,8,10,12
    for length in [6, 8, 10, 12]:
        for i in range(0, n - length + 1):
            window = s[i:i+length]
            # verify alternating pattern: e.g., CGCGCG, GCGCGC...
            ok = True
            for j in range(0, length, 2):
                pair = window[j:j+2]
                if pair not in patterns:
                    ok = False
                    break
            if ok:
                for j in range(i, i+length):
                    covered[j] = True
    covered_bases = sum(1 for x in covered if x)
    propensity = (covered_bases / n * 100) if n > 0 else 0.0
    return round(propensity, 3)

def fetch_record(gene, organism, db, rettype):
    term = f"{gene}[Gene] AND {organism}[Organism] AND RefSeq[Filter]"
    if db == "nuccore":
        term += " AND biomol_mrna[PROP]"
    with Entrez.esearch(db=db, term=term, retmax=5) as handle:
        rec = Entrez.read(handle)
    ids = rec.get("IdList", [])
    if not ids:
        return None, None
    chosen = ids[0]
    data = None
    with Entrez.efetch(db=db, id=chosen, rettype=rettype.lower(), retmode="text") as fh:
        data = fh.read()
    return chosen, data

def send_email(zip_bytes, csv_bytes, recipient):
    sender = os.getenv("SMTP_SENDER", "")
    pwd = os.getenv("SMTP_PASSWORD", "")
    host = os.getenv("SMTP_HOST", "smtp.gmail.com")
    port = int(os.getenv("SMTP_PORT", "465"))
    if not sender or not pwd:
        st.error("SMTP credentials not set in environment")
        return
    msg = EmailMessage()
    msg["Subject"] = "RefSeq Output"
    msg["From"] = sender
    msg["To"] = recipient
    msg.set_content("Attached sequence outputs and metadata.")
    msg.add_attachment(zip_bytes, maintype="application", subtype="zip", filename="sequences.zip")
    msg.add_attachment(csv_bytes, maintype="text", subtype="csv", filename="output.csv")
    with smtplib.SMTP_SSL(host, port) as smtp:
        smtp.login(sender, pwd)
        smtp.send_message(msg)

# Main logic
results = []
zip_buffer = io.BytesIO()

if st.button("Fetch and Generate Output"):
    genes = [g.strip() for g in gene_input.replace(",", "\\n").splitlines() if g.strip()]
    if not genes:
        st.error("No genes provided")
        st.stop()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zipf:
        for gene in genes:
            st.markdown(f"Processing {gene}...")
            db = "nuccore" if sequence_type == "Nucleotide" else "protein"
            rettype = "GenBank" if output_format == "GenBank" else "fasta"
            rec_id, raw = fetch_record(gene, organism, db, rettype)
            if raw is None:
                st.warning(f"No RefSeq found for {gene} in {organism}")
                continue
            # Save per-gene file structure
            safe_gene = re.sub(r\"[^A-Za-z0-9._-]\", \"_\", gene)
            ext = "gb" if output_format == "GenBank" else "fasta"
            zipf.writestr(f\"{safe_gene}/{safe_gene}_full.{ext}\", raw)
            # Initialize sequences
            seq_str = ""
            utr3 = ""
            if output_format == "FASTA":
                lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
                seq_str = \"\".join(ln for ln in lines if not ln.startswith(">"))
            else:
                try:
                    gbrec = SeqIO.read(io.StringIO(raw), "genbank")
                    seq_str = str(gbrec.seq)
                    if db == "nuccore":
                        cds_end = None
                        for feat in gbrec.features:
                            if feat.type == "CDS":
                                cds_end = int(feat.location.end)
                                break
                        if cds_end is not None and cds_end < len(gbrec.seq):
                            utr3 = str(gbrec.seq[cds_end:])
                            # write utr fasta
                            hdr = f\">{gene}|{organism}|3'UTR_from_{cds_end+1}\\n\"
                            wrapped = \"\\n\".join([utr3[i:i+70] for i in range(0, len(utr3), 70)])
                            zipf.writestr(f\"{safe_gene}/{safe_gene}_utr.fasta\", hdr + wrapped + \"\\n\")
                except Exception:
                    seq_str = \"\"
                    utr3 = \"\"
            # Compute stats for full
            full_stats = compute_base_stats(seq_str)
            zdna_full = zdna_propensity(seq_str)
            full_entry = {
                "Gene Name": gene,
                "Sequence": seq_str,
                **full_stats,
                "Z-DNA Propensity %": zdna_full
            }
            results.append(full_entry)
            # UTR entry
            if utr3:
                utr_stats = compute_base_stats(utr3)
                zdna_utr = zdna_propensity(utr3)
                utr_entry = {
                    "Gene Name": gene + " (3'UTR)",
                    "Sequence": utr3,
                    **utr_stats,
                    "Z-DNA Propensity %": zdna_utr
                }
                results.append(utr_entry)
    # Output
    if results:
        out_df = pd.DataFrame(results)
        st.markdown("### Output Table")
        st.dataframe(out_df, use_container_width=True)
        csv_bytes = out_df.to_csv(index=False).encode("utf-8")
        st.download_button("Download output.csv", csv_bytes, "output.csv", "text/csv")
        xlsx_buf = io.BytesIO()
        with pd.ExcelWriter(xlsx_buf, engine="xlsxwriter") as writer:
            out_df.to_excel(writer, index=False, sheet_name="Output")
        xlsx_bytes = xlsx_buf.getvalue()
        st.download_button("Download output.xlsx", xlsx_bytes,
                           "output.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        # add manifest to zip
        with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zipf:
            zipf.writestr("output.csv", csv_bytes)
        zip_buffer.seek(0)
        st.download_button("Download all (ZIP)", zip_buffer, "full_package.zip", "application/zip")
        if send_email and user_email:
            send_email(zip_buffer.getvalue(), csv_bytes, user_email)
