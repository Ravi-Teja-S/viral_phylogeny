import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from Bio import AlignIO, SeqIO, Phylo
from io import StringIO
from pathlib import Path

BASE = Path(__file__).resolve().parent

st.set_page_config(
    page_title="SARS-CoV-2 Phylogenomic Analysis",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

# ---- Global CSS ----
st.markdown("""
<style>
    [data-testid="stSidebar"] { background-color: #0f1117; }
    [data-testid="stSidebar"] * { color: #ffffff !important; }
    .metric-card {
        background: linear-gradient(135deg, #1e3a5f, #0d1b2a);
        border-radius: 12px;
        padding: 20px;
        text-align: center;
        border: 1px solid #2a4a7f;
    }
    .metric-card h2 { color: #4fc3f7; font-size: 2rem; margin: 0; }
    .metric-card p  { color: #90caf9; margin: 4px 0 0 0; font-size: 0.9rem; }
    .section-header {
        background: linear-gradient(90deg, #1a237e, #0d47a1);
        padding: 12px 20px;
        border-radius: 8px;
        color: white;
        font-size: 1.1rem;
        font-weight: 600;
        margin-bottom: 16px;
    }
    .variant-badge {
        display: inline-block;
        background: #1565c0;
        color: white;
        padding: 3px 10px;
        border-radius: 20px;
        font-size: 0.8rem;
        margin: 2px;
    }
    .stTabs [data-baseweb="tab"] { font-size: 1rem; font-weight: 500; }
    div[data-testid="stMetricValue"] { font-size: 2rem; color: #4fc3f7; }
</style>
""", unsafe_allow_html=True)

# ---- Variant metadata ----
VARIANTS_INFO = {
    "Wuhan":         {"lineage": "Original",   "accession": "NC_045512.2", "wave": "2019", "color": "#4caf50"},
    "Alpha_B117":    {"lineage": "B.1.1.7",    "accession": "MW422255",    "wave": "2020", "color": "#2196f3"},
    "Beta_B1351":    {"lineage": "B.1.351",    "accession": "MW598419",    "wave": "2020", "color": "#9c27b0"},
    "Gamma_P1":      {"lineage": "P.1",        "accession": "MT126808",    "wave": "2021", "color": "#ff9800"},
    "Delta_B16172":  {"lineage": "B.1.617.2",  "accession": "MZ208926",    "wave": "2021", "color": "#f44336"},
    "Omicron_BA1":   {"lineage": "BA.1",       "accession": "OM570283",    "wave": "2021", "color": "#00bcd4"},
    "Omicron_BA2":   {"lineage": "BA.2",       "accession": "OM617939",    "wave": "2022", "color": "#00acc1"},
    "Omicron_BA5":   {"lineage": "BA.5",       "accession": "OP093341",    "wave": "2022", "color": "#0097a7"},
    "Lambda_C37":    {"lineage": "C.37",       "accession": "MZ068159",    "wave": "2021", "color": "#ff5722"},
    "Mu_B1621":      {"lineage": "B.1.621",    "accession": "MZ344997",    "wave": "2021", "color": "#795548"},
    "Bat_outgroup":  {"lineage": "Outgroup",   "accession": "MN996532",    "wave": "—",    "color": "#607d8b"},
}

# ---- Data loaders ----
@st.cache_data
def load_mutations():
    alignment = AlignIO.read(BASE / "data/aligned/spike_proteins.aln", "fasta")
    seqs = {rec.id: rec.seq for rec in alignment}
    ref_name = next((n for n in seqs if "Wuhan" in n or "NC_045512" in n), None)
    ref = seqs[ref_name]
    rbd_start, rbd_end = 318, 541
    data = {}
    for name, seq in seqs.items():
        if name == ref_name:
            continue
        mutations = []
        for i in range(rbd_start, rbd_end):
            w, v = ref[i], seq[i]
            if w != v and w not in "-X" and v not in "-X":
                mutations.append({
                    "Position": i + 1,
                    "Wuhan AA": str(w),
                    "Variant AA": str(v),
                    "Mutation": f"{w}{i+1}{v}"
                })
        data[name] = mutations
    return ref_name, data

@st.cache_data
def load_sequences(path, fmt):
    return list(SeqIO.parse(path, fmt))

@st.cache_data
def load_rbd_matrix():
    alignment = AlignIO.read(BASE / "data/aligned/spike_sequences.aln", "fasta")
    names = [rec.id for rec in alignment]
    RBD_START, RBD_END = 957, 1623
    n = len(alignment)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            s1 = alignment[i].seq[RBD_START:RBD_END]
            s2 = alignment[j].seq[RBD_START:RBD_END]
            m = v = 0
            for a, b in zip(s1, s2):
                if a == "-" or b == "-":
                    continue
                v += 1
                if a == b:
                    m += 1
            matrix[i, j] = round(100 * m / v, 2) if v > 0 else 0
    return names, matrix

ref_name, mutation_data = load_mutations()

# ---- Sidebar navigation ----
with st.sidebar:
    st.markdown("## 🧬 SARS-CoV-2")
    st.markdown("### Phylogenomic Analysis")
    st.markdown("---")
    page = st.radio(
        "Navigate",
        ["🏠 Home", "🔍 Sequences", "🌳 Phylogenetic Tree", "🔥 RBD Heatmap", "🧪 RBD Mutations"],
        label_visibility="collapsed"
    )
    st.markdown("---")
    st.markdown("**Reference:** Wuhan (NC_045512.2)")
    st.markdown("**Variants analysed:** 10 + 1 outgroup")
    st.markdown("**RBD region:** AA 319–541")
    st.markdown("---")
    st.caption("Built with Biopython + Streamlit")


# ══════════════════════════════════════════════
# PAGE 1 — HOME
# ══════════════════════════════════════════════
if page == "🏠 Home":
    st.title("🧬 SARS-CoV-2 Phylogenomic Analysis")
    st.markdown("A comprehensive genomic analysis of SARS-CoV-2 variants — tracking evolutionary divergence through the **Spike protein's Receptor Binding Domain (RBD)**.")
    st.markdown("---")

    # ---- Top metrics ----
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.markdown('<div class="metric-card"><h2>11</h2><p>Sequences Analysed</p></div>', unsafe_allow_html=True)
    with c2:
        total_mut = sum(len(v) for v in mutation_data.values())
        st.markdown(f'<div class="metric-card"><h2>{total_mut}</h2><p>Total RBD Mutations</p></div>', unsafe_allow_html=True)
    with c3:
        max_var = max(mutation_data, key=lambda x: len(mutation_data[x]))
        st.markdown(f'<div class="metric-card"><h2>{len(mutation_data[max_var])}</h2><p>Most Mutated ({max_var.split("_")[0]})</p></div>', unsafe_allow_html=True)
    with c4:
        st.markdown('<div class="metric-card"><h2>223</h2><p>RBD Positions Scanned</p></div>', unsafe_allow_html=True)

    st.markdown("---")

    # ---- Pipeline overview ----
    st.markdown('<div class="section-header">📋 Analysis Pipeline</div>', unsafe_allow_html=True)
    steps = [
        ("1️⃣", "Sequence Retrieval",      "Fetch Spike CDS from NCBI GenBank for all 11 variants"),
        ("2️⃣", "Translation",             "Convert nucleotide sequences to amino acid (protein) sequences"),
        ("3️⃣", "Multiple Sequence Alignment", "Align all sequences using MAFFT to enable position-by-position comparison"),
        ("4️⃣", "Phylogenetic Tree",        "Build Neighbor-Joining tree using BLOSUM62 distance matrix"),
        ("5️⃣", "RBD Divergence Heatmap",  "Compute pairwise % identity across the RBD nucleotide region"),
        ("6️⃣", "RBD Mutation Landscape",  "Visualise per-position mutations vs Wuhan across all variants"),
        ("7️⃣", "Mutation Table",           "List exact amino acid changes per variant in the RBD"),
    ]
    for icon, title, desc in steps:
        st.markdown(f"{icon} &nbsp; **{title}** — {desc}")

    st.markdown("---")

    # ---- Variants table ----
    st.markdown('<div class="section-header">🦠 Variants in This Study</div>', unsafe_allow_html=True)
    rows = []
    for name, info in VARIANTS_INFO.items():
        mut_count = len(mutation_data.get(name, []))
        rows.append({
            "Variant": name,
            "Lineage": info["lineage"],
            "Accession": info["accession"],
            "Wave": info["wave"],
            "RBD Mutations vs Wuhan": mut_count if name != "Wuhan" else "—"
        })
    df_variants = pd.DataFrame(rows)
    st.dataframe(df_variants, use_container_width=True, hide_index=True)


# ══════════════════════════════════════════════
# PAGE 2 — SEQUENCES
# ══════════════════════════════════════════════
elif page == "🔍 Sequences":
    st.title("🔍 Spike Protein Sequences")
    st.markdown("Raw nucleotide CDS and translated amino acid sequences fetched from NCBI GenBank.")
    st.markdown("---")

    tab1, tab2 = st.tabs(["🧬 Nucleotide (DNA)", "🔤 Protein (Amino Acids)"])

    with tab1:
        st.markdown('<div class="section-header">Spike CDS — Nucleotide Sequences</div>', unsafe_allow_html=True)
        records = load_sequences(BASE / "data/raw/spike_sequences.fasta", "fasta")
        for rec in records:
            info = VARIANTS_INFO.get(rec.id, {})
            with st.expander(f"**{rec.id}**  |  {info.get('lineage','—')}  |  {info.get('accession','—')}  |  Length: {len(rec.seq)} nt"):
                col1, col2 = st.columns([1, 3])
                with col1:
                    st.metric("Length (nt)", len(rec.seq))
                    st.metric("Accession", info.get("accession", "—"))
                    st.metric("Wave", info.get("wave", "—"))
                with col2:
                    st.code(str(rec.seq), language=None)

    with tab2:
        st.markdown('<div class="section-header">Spike Protein — Amino Acid Sequences</div>', unsafe_allow_html=True)
        records = load_sequences(BASE / "data/proteins/spike_proteins.fasta", "fasta")
        for rec in records:
            info = VARIANTS_INFO.get(rec.id, {})
            with st.expander(f"**{rec.id}**  |  {info.get('lineage','—')}  |  Length: {len(rec.seq)} aa"):
                col1, col2 = st.columns([1, 3])
                with col1:
                    st.metric("Length (aa)", len(rec.seq))
                    st.metric("Accession", info.get("accession", "—"))
                with col2:
                    st.code(str(rec.seq), language=None)


# ══════════════════════════════════════════════
# PAGE 3 — PHYLOGENETIC TREE
# ══════════════════════════════════════════════
elif page == "🌳 Phylogenetic Tree":
    st.title("🌳 Phylogenetic Tree")
    st.markdown("Evolutionary relationships between SARS-CoV-2 variants rooted at the **Bat coronavirus outgroup**.")
    st.markdown("---")

    VARIANT_COLORS = {
        "Wuhan":        "#4caf50",
        "Alpha_B117":   "#2196f3",
        "Beta_B1351":   "#ce93d8",
        "Gamma_P1":     "#ffb300",
        "Delta_B16172": "#ef5350",
        "Omicron_BA1":  "#00e5ff",
        "Omicron_BA2":  "#40c4ff",
        "Omicron_BA5":  "#80deea",
        "Lambda_C37":   "#ff7043",
        "Mu_B1621":     "#bcaaa4",
        "Bat_outgroup": "#78909c",
    }

    tree = Phylo.read(str(BASE / "results/phylogeny_tree.nwk"), "newick")
    n_leaves = len(tree.get_terminals())
    leaves   = tree.get_terminals()
    depths   = tree.depths(unit_branch_lengths=False)

    fig, ax = plt.subplots(figsize=(14, max(8, n_leaves * 0.8)))

    def label_colors(name):
        return VARIANT_COLORS.get(name, "black")

    Phylo.draw(tree, axes=ax, do_show=False, label_colors=label_colors)

    # overlay colored dots on leaf nodes
    for idx, leaf in enumerate(leaves):
        ax.scatter(
            depths[leaf], idx + 1,
            s=60, color=VARIANT_COLORS.get(leaf.name, "#333333"),
            zorder=5, linewidths=0.8, edgecolors="white"
        )

    plt.tight_layout()

    col_tree, col_info = st.columns([3, 1])

    with col_tree:
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)

        st.markdown("**Node Info** — select a variant:")
        selected_node = st.selectbox(
            "",
            ["— select —"] + [lf.name for lf in leaves],
            label_visibility="collapsed",
            key="tree_node_sel"
        )
        if selected_node and selected_node != "— select —":
            info = VARIANTS_INFO.get(selected_node, {})
            muts = mutation_data.get(selected_node, [])
            vcolor = VARIANT_COLORS.get(selected_node, "#333")
            st.markdown(
                f'<div style="border-left: 4px solid {vcolor}; padding: 10px 16px; '
                f'background:#f9f9f9; border-radius:6px; margin-top:8px">'
                f'<b style="font-size:1rem">{selected_node}</b><br>'
                f'Lineage: <b>{info.get("lineage", "—")}</b><br>'
                f'Accession: {info.get("accession", "—")}<br>'
                f'Wave: {info.get("wave", "—")}<br>'
                f'RBD Mutations vs Wuhan: <b>{len(muts)}</b>'
                f'</div>',
                unsafe_allow_html=True
            )

    with col_info:
        st.markdown('<div class="section-header">🎨 Variant Legend</div>',
                    unsafe_allow_html=True)
        for vname, vcolor in VARIANT_COLORS.items():
            vinfo = VARIANTS_INFO.get(vname, {})
            st.markdown(
                f'<div style="display:flex;align-items:center;margin-bottom:6px">'
                f'<span style="color:{vcolor};font-size:1.3rem;margin-right:8px">●</span>'
                f'<span style="font-size:0.82rem"><b>{vname}</b><br>'
                f'<span style="color:#78909c">{vinfo.get("lineage","—")}</span>'
                f'</span></div>',
                unsafe_allow_html=True
            )
        st.markdown("---")
        st.markdown('<div class="section-header">📐 Method</div>',
                    unsafe_allow_html=True)
        st.markdown("**Algorithm:** Neighbor-Joining  ")
        st.markdown("**Distance:** BLOSUM62  ")
        st.markdown("**Root:** Bat coronavirus outgroup  ")
        st.markdown("**Labels:** Leaf = variant name, Inner = internal node")

    st.markdown("---")
    st.markdown('<div class="section-header">🔍 Evolutionary Groupings</div>',
                unsafe_allow_html=True)
    g1, g2, g3 = st.columns(3)
    with g1:
        st.markdown("**🟢 Closest to Wuhan**")
        st.markdown("- Gamma (P.1)\n- Wuhan (reference)")
    with g2:
        st.markdown("**🟡 Intermediate Divergence**")
        st.markdown("- Alpha (B.1.1.7)\n- Mu (B.1.621)\n- Delta (B.1.617.2)\n- Beta (B.1.351)\n- Lambda (C.37)")
    with g3:
        st.markdown("**🔴 Most Diverged**")
        st.markdown("- Omicron BA.1\n- Omicron BA.2\n- Omicron BA.5")


# ══════════════════════════════════════════════
# PAGE 4 — RBD HEATMAP
# ══════════════════════════════════════════════
elif page == "🔥 RBD Heatmap":
    st.title("🔥 RBD Divergence Heatmap")
    st.markdown("Visual analysis of the **Receptor Binding Domain** — the region of the Spike protein that binds to human ACE2 receptors.")
    st.markdown("---")

    tab1, tab2 = st.tabs(["🧬 Mutation Landscape (vs Wuhan)", "📊 Pairwise % Identity"])

    with tab1:
        st.markdown('<div class="section-header">Per-position RBD Mutations Relative to Wuhan (AA 319–541)</div>', unsafe_allow_html=True)
        st.markdown("Each **red cell** = a mutation at that amino acid position in that variant compared to Wuhan.")
        heatmap1 = BASE / "figures/rbd_heatmap.png"
        if heatmap1.exists():
            st.image(str(heatmap1), use_container_width=True)
        else:
            st.warning("rbd_heatmap.png not found. Run script 6 to generate it.")

        st.markdown("---")
        st.markdown("""
**Key observations:**
- **Omicron BA.5** has the most spread-out mutations across the RBD
- **Alpha** has a single concentrated mutation (N501Y)
- **Omicron variants** share many mutations in the 417–505 region
        """)

    with tab2:
        st.markdown('<div class="section-header">Pairwise % Identity Across RBD (Nucleotide)</div>', unsafe_allow_html=True)
        heatmap2 = BASE / "figures/rbd_heatmap1.png"
        if heatmap2.exists():
            st.image(str(heatmap2), use_container_width=True)
        else:
            with st.spinner("Computing pairwise identity matrix..."):
                names, matrix = load_rbd_matrix()
            fig, ax = plt.subplots(figsize=(9, 7))
            im = ax.imshow(matrix, cmap="YlOrRd_r", vmin=85, vmax=100)
            ax.set_xticks(range(len(names)))
            ax.set_yticks(range(len(names)))
            ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(names, fontsize=9)
            for i in range(len(names)):
                for j in range(len(names)):
                    ax.text(j, i, f"{matrix[i,j]:.1f}", ha="center", va="center", fontsize=7,
                            color="black" if matrix[i,j] > 90 else "white")
            plt.colorbar(im, ax=ax, label="% Identity")
            ax.set_title("RBD Pairwise % Identity", fontsize=13, fontweight="bold")
            plt.tight_layout()
            st.pyplot(fig)

        st.markdown("""
**How to read:**
- **Darker = more similar** to each other
- Diagonal is always 100% (self-comparison)
- Omicron variants cluster together with high mutual identity
        """)


# ══════════════════════════════════════════════
# PAGE 5 — RBD MUTATIONS
# ══════════════════════════════════════════════
elif page == "🧪 RBD Mutations":
    st.title("🧪 RBD Mutation Tracker")
    st.markdown("Exact amino acid changes in the **Receptor Binding Domain (AA 319–541)** of each variant compared to Wuhan.")
    st.markdown("---")

    # ---- Summary bar chart ----
    st.markdown('<div class="section-header">📊 Mutation Count — All Variants vs Wuhan</div>', unsafe_allow_html=True)
    counts = {v: len(m) for v, m in mutation_data.items()}
    chart_df = pd.DataFrame.from_dict(counts, orient="index", columns=["RBD Mutations"]).sort_values("RBD Mutations", ascending=False)
    st.bar_chart(chart_df, color="#4fc3f7")

    st.markdown("---")

    # ---- Variant selector ----
    col_sel, col_info = st.columns([1, 2])
    with col_sel:
        st.markdown('<div class="section-header">Select Variant</div>', unsafe_allow_html=True)
        selected = st.selectbox("", list(mutation_data.keys()), label_visibility="collapsed")

    mutations = mutation_data[selected]
    info = VARIANTS_INFO.get(selected, {})

    with col_info:
        st.markdown('<div class="section-header">Variant Info</div>', unsafe_allow_html=True)
        i1, i2, i3 = st.columns(3)
        i1.metric("Lineage", info.get("lineage", "—"))
        i2.metric("Accession", info.get("accession", "—"))
        i3.metric("RBD Mutations", len(mutations))

    st.markdown("---")

    # ---- Mutation detail ----
    st.markdown(f'<div class="section-header">🔬 {selected} — RBD Mutations vs Wuhan</div>', unsafe_allow_html=True)

    if mutations:
        df = pd.DataFrame(mutations)

        # colour Variant AA column red
        def highlight(row):
            return ["", "", "background-color:#fff0f0; color:#c0392b; font-weight:bold", "background-color:#fdecea; color:#b71c1c; font-weight:bold"]

        st.dataframe(df.style.apply(highlight, axis=1), use_container_width=True, hide_index=True)

        st.markdown("---")
        st.markdown("**All mutations:** " + "  ".join([f"`{m['Mutation']}`" for m in mutations]))

        # ---- Position distribution chart ----
        st.markdown("---")
        st.markdown('<div class="section-header">📍 Mutation Position Distribution</div>', unsafe_allow_html=True)
        pos_df = pd.DataFrame({"Position": [m["Position"] for m in mutations], "Count": [1] * len(mutations)})
        pos_df = pos_df.set_index("Position")
        st.bar_chart(pos_df, color="#ef5350")

    else:
        st.success(f"✅ No RBD mutations detected in **{selected}** compared to Wuhan.")

    st.markdown("---")

    # ---- All variants comparison ----
    st.markdown('<div class="section-header">📋 All Variants — Full Mutation Summary</div>', unsafe_allow_html=True)
    for variant, muts in mutation_data.items():
        vinfo = VARIANTS_INFO.get(variant, {})
        with st.expander(f"**{variant}**  |  {vinfo.get('lineage','—')}  |  {len(muts)} mutations"):
            if muts:
                st.dataframe(pd.DataFrame(muts), use_container_width=True, hide_index=True)
                st.markdown("**Mutations:** " + "  ".join([f"`{m['Mutation']}`" for m in muts]))
            else:
                st.success("No RBD mutations vs Wuhan.")
