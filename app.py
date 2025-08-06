
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
    return bool(re.match(r"^[A-Z]{2,3}_d+(.d+)?$", token))

def compute_base_stats(seq: str):
    seq = seq.upper()
    A = seq.count("A")
    T = seq.count("T")
    G = seq.count("G")
    C = seq.count("C")
    total = len(seq)
    at = A + T
    gc = G + C
    def pct(x): return (x / total * 100) if total>0 else 0.0
    return {
        "Length": total,
        "A Count": A,
        "T Count": T,
        "G Count": G,
        "C Count": C,
        "A %": round(pct(A),3),
        "T %": round(pct(T),3),
        "G %": round(pct(G),3),
        "C %": round(pct(C),3),
        "GC %": round(pct(gc),3),
        "AT %": round(pct(at),3),
    }

def zdna_propensity(seq: str):
    # Simple heuristic: identify alternating CG or GC runs of length>=6 (i.e., 3 repeats)
    s = seq.upper()
    n = len(s)
    covered = [False]*n
    patterns = ["CG","GC"]
    for i in range(n-5):  # window size at least 6
        window = s[i:i+6]
        # check if window is alternating CGCGCG or GCGCGC
        if window in ("CGCGCG","GCGCGC"):
            for j in range(i,i+6):
                covered[j]=True
    # extend to longer alternating beyond 6
    # check longer even lengths up to 12
    for length in [8,10,12]:
        for i in range(n-length+1):
            window=s[i:i+length]
            # build expected patterns
            if all(window[k:k+2] in patterns for k in range(0,length,2)) and all(window[k]!=window[k+1] for k in range(length-1)):
                for j in range(i,i+length):
                    covered[j]=True
    covered_bases = sum(1 for x in covered if x)
    propensity = (covered_bases / n * 100) if n>0 else 0.0
    return round(propensity,3)

def fetch_record(gene, organism, db, rettype):
    # search by gene and organism
    term = f"{gene}[Gene] AND {organism}[Organism] AND RefSeq[Filter]"
    if db=="nuccore":
        term += " AND biomol_mrna[PROP]"
    handle = Entrez.esearch(db=db, term=term, retmax=5)
    rec = Entrez.read(handle)
    ids = rec.get("IdList",[])
    if not ids:
        return None, None
    chosen = ids[0]
    data = None
    with Entrez.efetch(db=db, id=chosen, rettype=rettype.lower(), retmode="text") as fh:
        data=fh.read()
    return chosen, data

def send_email(zip_bytes, csv_bytes, recipient):
    sender=os.getenv("SMTP_SENDER","")
    pwd=os.getenv("SMTP_PASSWORD","")
    host=os.getenv("SMTP_HOST","smtp.gmail.com")
    port=int(os.getenv("SMTP_PORT","465"))
    if not sender or not pwd:
        st.error("SMTP credentials not set in environment")
        return
    msg=EmailMessage()
    msg["Subject"]="RefSeq Output"
    msg["From"]=sender
    msg["To"]=recipient
    msg.set_content("Attached sequence outputs and metadata.")
    msg.add_attachment(zip_bytes, maintype="application", subtype="zip", filename="sequences.zip")
    msg.add_attachment(csv_bytes, maintype="text", subtype="csv", filename="output.csv")
    with smtplib.SMTP_SSL(host, port) as smtp:
        smtp.login(sender, pwd)
        smtp.send_message(msg)

# Main logic
results = []
zip_buffer=io.BytesIO()
with zipfile.ZipFile(zip_buffer,"a",zipfile.ZIP_DEFLATED) as zipf:
    if st.button("Fetch and Generate Output"):
        genes=[g.strip() for g in gene_input.replace(",","\\n").splitlines() if g.strip()]
        if not genes:
            st.error("No genes provided")
            st.stop()
        for gene in genes:
            st.markdown(f"Processing {gene}...")
            db="nuccore" if sequence_type=="Nucleotide" else "protein"
            rettype="GenBank" if output_format=="GenBank" else "fasta"
            rec_id, raw = fetch_record(gene, organism, db, rettype)
            if raw is None:
                st.warning(f"No RefSeq found for {gene} in {organism}")
                continue
            # Save per-gene file
            safe_gene=re.sub(r"[^A-Za-z0-9._-]","_", gene)
            ext="gb" if output_format=="GenBank" else "fasta"
            zipf.writestr(f"{safe_gene}{safe_gene}_full.{ext}", raw)
            # Extract sequence string
            seq_str=""
            if output_format=="FASTA":
                lines=[ln.strip() for ln in raw.splitlines() if ln.strip()]
                seq_str="".join(ln for ln in lines if not ln.startswith(">"))
            else:
                try:
                    gbrec=SeqIO.read(io.StringIO(raw),"genbank")
                    seq_str=str(gbrec.seq)
                    # also extract 3'UTR if nucleotide
                    if db=="nuccore":
                        # find CDS end
                        cds_end=None
                        for feat in gbrec.features:
                            if feat.type=="CDS":
                                cds_end=int(feat.location.end)
                                break
                        if cds_end and cds_end < len(gbrec.seq):
                            utr3=str(gbrec.seq[cds_end:])
                            # write utr file
                            hdr=f">{gene}|{organism}|3'UTR_from_{cds_end+1}n"
                            wrapped="n".join([utr3[i:i+70] for i in range(0,len(utr3),70)])
                            zipf.writestr(f"{safe_gene}/{safe_gene}_utr.fasta", hdr+wrapped+"n")
                        else:
                            utr3=""
                    else:
                        utr3=""
                except Exception:
                    seq_str=""
                    utr3=""
            # Compute stats for full and utr
            full_stats=compute_base_stats(seq_str)
            zdna=zdna_propensity(seq_str)
            # add z-dna to full
            full_entry={
                "Gene Name": gene,
                "Sequence": seq_str,
                **full_stats,
                "Z-DNA Propensity %": zdna
            }
            results.append(full_entry)
            # utr stats if exists
            if output_format=="GenBank" and db=="nuccore" and 'utr3' in locals() and utr3:
                utr_stats=compute_base_stats(utr3)
                utr_zdna=zdna_propensity(utr3)
                utr_entry={
                    "Gene Name": gene + " (3'UTR)",
                    "Sequence": utr3,
                    **utr_stats,
                    "Z-DNA Propensity %": utr_zdna
                }
                results.append(utr_entry)
# After button press, present
if results:
    out_df=pd.DataFrame(results)
    # Calculate percentages columns if missing
    st.markdown("### Output Table")
    st.dataframe(out_df, use_container_width=True)
    csv_bytes=out_df.to_csv(index=False).encode("utf-8")
    st.download_button("Download output.csv", csv_bytes, "output.csv", "text/csv")
    # Package progress-like workbook with output as sheet
    xlsx_buf=io.BytesIO()
    with pd.ExcelWriter(xlsx_buf, engine="xlsxwriter") as writer:
        out_df.to_excel(writer, index=False, sheet_name="Output")
    xlsx_bytes=xlsx_buf.getvalue()
    st.download_button("Download output.xlsx", xlsx_bytes,
                       "output.xlsx","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    # ZIP already has sequence files; also add a manifest CSV inside ZIP
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zipf:
        zipf.writestr("output.csv", csv_bytes)
    # Final downloadable ZIP
    zip_buffer.seek(0)
    st.download_button("Download all (ZIP)", zip_buffer, "full_package.zip", "application/zip")
    if send_email and user_email:
        send_email(zip_buffer.getvalue(), csv_bytes, user_email)
