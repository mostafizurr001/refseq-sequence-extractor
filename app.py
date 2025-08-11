mport os
import re
import io
import zipfile
import smtplib
import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
from email.message import EmailMessage

# ---------------------------
# NCBI Entrez Setup
# ---------------------------
# IMPORTANT: set your email (required by NCBI)
Entrez.email = os.getenv("NCBI_EMAIL", "your_email@example.com")
# Optional: use an API key to increase rate limits
# Entrez.api_key = os.getenv("NCBI_API_KEY")

st.set_page_config(page_title="NCBI RefSeq Multi-Gene Extractor", layout="centered")

st.title("🧬 NCBI RefSeq Multi-Gene Sequence Extractor")
st.markdown("""
Search and download sequences from NCBI RefSeq for multiple genes with metadata and annotations.  
This tool prefers **RefSeq transcript variant 1** and **avoids 'transcript variant X'** when possible.  
It can automatically extract the **3′UTR** from mRNA (GenBank format) and compute nucleotide composition & GC%.
""")

# ---------------------------
# Inputs
# ---------------------------
gene_input = st.text_area(
    "🔍 Enter Gene Symbols or RefSeq Accessions (comma or newline separated)",
    height=150,
    help="Examples: BRCA1, TP53, NM_000546. One per line or comma-separated."
)

common_organisms = [
    "Homo sapiens (Human)",
    "Mus musculus (Mouse)",
    "Rattus norvegicus (Rat)",
    "Danio rerio (Zebrafish)",
    "Drosophila melanogaster (Fruit fly)",
    "Caenorhabditis elegans (Nematode)",
    "Saccharomyces cerevisiae (Yeast)",
    "Escherichia coli (E. coli)",
    "Arabidopsis thaliana (Thale cress)",
    "Bos taurus (Cow)",
    "Oryza sativa (Rice)",
    "Zea mays (Maize)",
    "Schizosaccharomyces pombe (Fission yeast)",
    "Plasmodium falciparum (Malaria parasite)",
    "Chlamydomonas reinhardtii (Green algae)",
    "Other (type manually)"
]
organism_choice = st.selectbox(
    "🌱 Organism Name",
    common_organisms,
    help="Select a common organism or choose 'Other (type manually)'."
)
if organism_choice == "Other (type manually)":
    organism = st.text_input("Type Organism Name (e.g., Homo sapiens)", "")
else:
    organism = organism_choice.split(" (")[0]

sequence_type = st.selectbox("📦 Sequence Type", ["Nucleotide", "Protein"])
output_format = st.selectbox("📄 Output Format", ["GenBank", "FASTA"], help="GenBank is required for 3′UTR extraction.")
selected_types = st.multiselect(
    "🔖 Filter Feature Types (GenBank only)",
    ["CDS", "gene", "exon", "mRNA", "3'UTR", "5'UTR"],
    default=["CDS", "gene"]
)

prefer_variant1 = st.checkbox("Prefer transcript variant 1 (recommended)", value=True)
exclude_variant_x = st.checkbox("Exclude 'transcript variant X*' records", value=True)

send_email = st.checkbox("✉️ Email me the ZIP and metadata (uses SMTP)")
if send_email:
    user_email = st.text_input("📧 Your Email")
else:
    user_email = ""

# ---------------------------
# Helpers
# ---------------------------
def is_refseq_accession(token: str) -> bool:
    # Accepts RefSeq-like accessions: NM_, NR_, NP_, NG_, NC_, XM_, XR_, XP_, etc.
    return bool(re.match(r"^[A-Z]{2,3}_\d+(\.\d+)?$", token))

def nucleotide_counts(seq: str):
    seq = seq.upper()
    return {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C"),
        "Length": len(seq),
        "GC%": (seq.count("G") + seq.count("C")) / len(seq) * 100 if len(seq) > 0 else 0.0
    }

def choose_variant1_id(id_list, db="nuccore"):
    """
    Given a list of nuccore IDs, fetch their summaries and pick one that
    contains 'transcript variant 1' (and does not contain 'transcript variant X').
    If none match, return the first ID.
    """
    if not id_list:
        return None
    try:
        handle = Entrez.esummary(db=db, id=",".join(id_list))
        summaries = Entrez.read(handle)
        handle.close()
        best_id = None
        for s in summaries:
            title = s.get("Title", "") or s.get("Caption", "")
            uid = s.get("Id") or s.get("Uid")
            title_l = title.lower()
            if "transcript variant 1" in title_l:
                if exclude_variant_x and "transcript variant x" in title_l:
                    continue
                best_id = uid
                break
        if best_id:
            return best_id
        if exclude_variant_x:
            for s in summaries:
                title = s.get("Title", "") or s.get("Caption", "")
                uid = s.get("Id") or s.get("Uid")
                if "transcript variant x" not in (title or "").lower():
                    return uid
        return summaries[0].get("Id") or summaries[0].get("Uid")
    except Exception:
        return id_list[0]

def search_ids_for_gene(gene: str, organism: str, db: str, outfmt: str):
    """
    Search nuccore or protein for records matching gene + organism + RefSeq.
    For nucleotide mRNA GenBank, try to target transcript variant 1.
    """
    base = f"{gene}[Gene] AND {organism}[Organism] AND RefSeq[Filter]"
    if db == "nuccore":
        base += " AND biomol_mrna[PROP]"
        if prefer_variant1:
            base += ' AND "transcript variant 1"[Title]'
        if exclude_variant_x:
            base += ' NOT "transcript variant X"[Title]'
    with Entrez.esearch(db=db, term=base, retmax=20) as handle:
        rec = Entrez.read(handle)
    ids = rec.get("IdList", [])
    if not ids and prefer_variant1:
        with Entrez.esearch(db=db, term=f"{gene}[Gene] AND {organism}[Organism] AND RefSeq[Filter] AND biomol_mrna[PROP]", retmax=20) as handle:
            rec = Entrez.read(handle)
        ids = rec.get("IdList", [])
    if not ids:
        with Entrez.esearch(db=db, term=f"{gene}[All Fields] AND {organism}[Organism] AND RefSeq[Filter]", retmax=20) as handle:
            rec = Entrez.read(handle)
        ids = rec.get("IdList", [])
    if ids and db == "nuccore":
        chosen = choose_variant1_id(ids, db="nuccore")
        return [chosen] if chosen else ids[:1]
    return ids[:1]

def fetch_text_record(db: str, rec_id: str, rettype: str):
    with Entrez.efetch(db=db, id=rec_id, rettype=rettype.lower(), retmode="text") as h:
        return h.read()

def extract_utr3_from_gb(gb_text: str):
    """
    Return (utr3_seq, cds_start, cds_end, full_seq) from a GenBank text.
    Uses CDS feature end to define 3'UTR start (end+1 through end of sequence).
    """
    record_io = io.StringIO(gb_text)
    gb = SeqIO.read(record_io, "genbank")
    full_seq = str(gb.seq)
    cds_end = None
    cds_start = None
    for feat in gb.features:
        if feat.type == "CDS":
            cds_start = int(feat.location.start)
            cds_end = int(feat.location.end)
            break
    utr3 = ""
    if cds_end is not None and cds_end < len(gb.seq):
        utr3 = str(gb.seq[cds_end:])
    return utr3, cds_start, cds_end, full_seq, gb

def send_email_with_attachment(to_email, zip_bytes, csv_bytes=None, xlsx_bytes=None, filename="RefSeq_sequences.zip"):
    sender = os.getenv("SMTP_SENDER", "")
    passwd = os.getenv("SMTP_PASSWORD", "")
    smtp_host = os.getenv("SMTP_HOST", "smtp.gmail.com")
    smtp_port = int(os.getenv("SMTP_PORT", "465"))
    if not sender or not passwd:
        raise RuntimeError("SMTP credentials not set. Please set SMTP_SENDER and SMTP_PASSWORD environment variables.")
    msg = EmailMessage()
    msg["Subject"] = "Your RefSeq Sequence Results"
    msg["From"] = sender
    msg["To"] = to_email
    msg.set_content("Please find attached your RefSeq sequences and metadata.")
    msg.add_attachment(zip_bytes, maintype="application", subtype="zip", filename=filename)
    if csv_bytes:
        msg.add_attachment(csv_bytes, maintype="text", subtype="csv", filename="metadata.csv")
    if xlsx_bytes:
        msg.add_attachment(xlsx_bytes, maintype="application", subtype="vnd.openxmlformats-officedocument.spreadsheetml.sheet", filename="metadata.xlsx")
    with smtplib.SMTP_SSL(smtp_host, smtp_port) as smtp:
        smtp.login(sender, passwd)
        smtp.send_message(msg)

# ---------------------------
# Main
# ---------------------------
metadata_rows = []
zip_buffer = io.BytesIO()

if st.button("🚀 Fetch Sequences"):
    genes = [g.strip() for g in gene_input.replace(",", "\n").splitlines() if g.strip()]
    if not genes:
        st.error("Please enter at least one gene symbol or RefSeq accession.")
        st.stop()

    db = "nuccore" if sequence_type == "Nucleotide" else "protein"
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED) as zipf:
        progress = st.progress(0)
        for idx, gene in enumerate(genes):
            st.markdown(f"🔎 Searching for `{gene}`...")
            try:
                ids = []
                if is_refseq_accession(gene):
                    search_terms = [
                        f"{gene}[Accession]",
                        f"{gene}[All Fields] AND RefSeq[Filter]",
                        gene
                    ]
                    for q in search_terms:
                        with Entrez.esearch(db=db, term=q, retmax=1) as handle:
                            rec = Entrez.read(handle)
                            if rec["IdList"]:
                                ids = rec["IdList"]
                                break
                else:
                    ids = search_ids_for_gene(gene, organism, db, output_format)

                if not ids:
                    st.warning(f"⚠️ No sequence found for `{gene}`.")
                    progress.progress((idx + 1) / len(genes))
                    continue

                rec_id = ids[0]
                data_text = fetch_text_record(db, rec_id, rettype=output_format)

                filename = f"{gene}_{organism.replace(' ', '_')}.{output_format.lower()}"
                zipf.writestr(filename, data_text)

                if db == "nuccore" and output_format == "GenBank":
                    utr3, cds_start, cds_end, full_seq, gb = extract_utr3_from_gb(data_text)
                    full_stats = nucleotide_counts(full_seq)
                    utr_stats = nucleotide_counts(utr3) if utr3 else {"A":0,"T":0,"G":0,"C":0,"Length":0,"GC%":0.0}

                    if utr3:
                        utr_fname = f"{gene}_{organism.replace(' ', '_')}_3UTR.fasta"
                        fasta_hdr = f">{gene}|{organism}|3'UTR_from_{cds_end+1 if cds_end is not None else 'NA'}\n"
                        wrapped = "\n".join([utr3[i:i+70] for i in range(0, len(utr3), 70)])
                        zipf.writestr(utr_fname, fasta_hdr + wrapped + "\n")

                    with st.expander(f"📘 Annotations for `{gene}`"):
                        for feat in gb.features:
                            if not selected_types or feat.type in selected_types:
                                st.code(f"{feat.type} | {feat.location} | Qualifiers: {feat.qualifiers}")

                    metadata_rows.append({
                        "Gene": gene,
                        "Accession": gb.id,
                        "Region": "FULL",
                        "Sequence": full_seq,
                        **full_stats
                    })
                    metadata_rows.append({
                        "Gene": gene,
                        "Accession": gb.id,
                        "Region": "3'UTR",
                        "Sequence": utr3,
                        **utr_stats
                    })
                else:
                    if output_format == "FASTA":
                        lines = [ln.strip() for ln in data_text.splitlines() if ln.strip()]
                        seq = "".join(ln for ln in lines if not ln.startswith(">"))
                        stats = nucleotide_counts(seq)
                        sequence_out = seq
                    else:
                        stats = {"A":0,"T":0,"G":0,"C":0,"Length":0,"GC%":0.0}
                        sequence_out = ""

                    metadata_rows.append({
                        "Gene": gene,
                        "Accession": rec_id,
                        "Region": "FULL",
                        "Sequence": sequence_out,
                        **stats
                    })

            except Exception as e:
                st.error(f"❌ Error for `{gene}`: {e}")
            finally:
                progress.progress((idx + 1) / len(genes))

    st.success("✅ All available sequences retrieved.")
    zip_buffer.seek(0)
    st.download_button("📦 Download All Sequences as ZIP", zip_buffer, "RefSeq_sequences.zip", "application/zip")

    if metadata_rows:
        st.markdown("📊 **Metadata & Nucleotide Composition**")
        df = pd.DataFrame(metadata_rows)
        st.dataframe(df, use_container_width=True)

        csv_bytes = df.to_csv(index=False).encode("utf-8")
        xlsx_buf = io.BytesIO()
        with pd.ExcelWriter(xlsx_buf, engine="xlsxwriter") as writer:
            df.to_excel(writer, index=False, sheet_name="Metadata")
        xlsx_bytes = xlsx_buf.getvalue()

        st.download_button("⬇️ Download Metadata as CSV", csv_bytes, "metadata.csv", "text/csv")
        st.download_button("⬇️ Download Metadata as Excel", xlsx_bytes, "metadata.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        if send_email and user_email:
            try:
                send_email_with_attachment(
                    to_email=user_email,
                    zip_bytes=zip_buffer.getvalue(),
                    csv_bytes=csv_bytes,
                    xlsx_bytes=xlsx_bytes
                )
                st.success(f"📧 Results sent to {user_email}")
            except Exception as e:
                st.error(f"❌ Failed to send email: {e}")
