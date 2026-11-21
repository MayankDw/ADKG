"""
Export the AlzKG Explorer page (ad_kg_explorer.html, written by visualize_ad_kg.py) to a
static Word document: pipeline summary tables, a screenshot of the core subgraph, and
full node and edge tables for that subgraph.

    python knowledge_graph/export_explorer_docx.py [explorer.html] [output.docx]

Defaults: datasets/data/kg/ad_kg_explorer.html -> datasets/data/kg/AlzKG Explorer.docx.
Needs python-docx and Google Chrome or Microsoft Edge (headless, for the graph image;
the page loads vis-network from cdnjs, so an internet connection is required).
"""
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = Path(__file__).resolve().parents[1]
KG_DIR = ROOT / "datasets" / "data" / "kg"
html = Path(sys.argv[1]) if len(sys.argv) > 1 else KG_DIR / "ad_kg_explorer.html"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else html.with_name("AlzKG Explorer.docx")

BROWSERS = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]
# Show only the graph canvas, sized to fill the browser window.
GRAPH_ONLY_CSS = """<style>
header, section[aria-labelledby="how"], .bench-head, #filters, .card, .foot { display: none !important }
body, .wrap { margin: 0 !important; padding: 0 !important; max-width: none !important; background: #fff !important }
.bench { display: block !important; border: 0 !important; border-radius: 0 !important }
#graph { height: 900px !important; width: 1200px !important; background: #fff !important }
</style>"""

s = html.read_text(encoding="utf8")
D = json.loads(re.search(r"const D = (\{.*?\});\n", s).group(1))
S = D["summary"]
TYPES = {"disease": "Disease", "gene/protein": "Gene / protein", "variant": "Variant",
         "drug/chemical/compound": "Drug / compound", "clinical_intervention": "Trial intervention",
         "effect/phenotype": "Phenotype"}
fmt = lambda n: f"{n:,}"


def screenshot_graph(tmp: Path) -> Path:
    browser = next((b for b in BROWSERS if Path(b).exists()), None) or shutil.which("chrome") or shutil.which("msedge")
    if not browser:
        sys.exit("Chrome or Edge is needed to render the graph image.")
    page = tmp / "graph_only.html"
    page.write_text(s.replace("</style>", "</style>\n" + GRAPH_ONLY_CSS, 1), encoding="utf8")
    png = tmp / "graph.png"
    subprocess.run([browser, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                    "--force-device-scale-factor=2", "--window-size=1200,900",
                    "--virtual-time-budget=15000", f"--user-data-dir={tmp / 'profile'}",
                    f"--screenshot={png}", page.as_uri()],
                   check=True, capture_output=True, timeout=120)
    if not png.exists():
        sys.exit("Browser did not produce a screenshot.")
    return png


tmpdir = tempfile.TemporaryDirectory()
graph_png = str(screenshot_graph(Path(tmpdir.name)))

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Inches(8.5), Inches(11)
for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
    setattr(sec, m, Inches(1))

st = doc.styles
st["Normal"].font.name = "Calibri"; st["Normal"].font.size = Pt(11)
for name, size in (("Title", 26), ("Heading 1", 16), ("Heading 2", 13)):
    f = st[name].font; f.name = "Georgia"; f.size = Pt(size); f.bold = True
    f.color.rgb = RGBColor(0x18, 0x21, 0x2C)
    rf = st[name].element.get_or_add_rPr().get_or_add_rFonts()
    for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
        rf.attrib.pop(qn(a), None)
    rf.set(qn("w:ascii"), "Georgia"); rf.set(qn("w:hAnsi"), "Georgia")
    st[name].paragraph_format.space_before = Pt(14 if name != "Title" else 0)
    st[name].paragraph_format.space_after = Pt(6)
    st[name].paragraph_format.keep_with_next = True

def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd"); sh.set(qn("w:val"), "clear"); sh.set(qn("w:color"), "auto"); sh.set(qn("w:fill"), hexcolor)
    tcPr.append(sh)

def table(headers, rows, widths, mono_cols=(), size=9.5):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"; t.alignment = WD_TABLE_ALIGNMENT.CENTER; t.autofit = False
    for i, h in enumerate(headers):
        c = t.rows[0].cells[i]; c.text = ""
        r = c.paragraphs[0].add_run(h); r.bold = True; r.font.size = Pt(size)
        shade(c, "ECEFF3")
    # repeat header row on each page
    trPr = t.rows[0]._tr.get_or_add_trPr(); th = OxmlElement("w:tblHeader"); th.set(qn("w:val"), "true"); trPr.append(th)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run(str(v)); r.font.size = Pt(size)
            if i in mono_cols: r.font.name = "Consolas"
    for row in t.rows:
        for i, w in enumerate(widths):
            row.cells[i].width = Inches(w)
    if len(rows) <= 12:  # keep short tables on one page
        for row in t.rows[:-1]:
            for c in row.cells:
                for pp in c.paragraphs: pp.paragraph_format.keep_with_next = True
    doc.add_paragraph()
    return t

def para(text, italic=False, color=None, size=None):
    p = doc.add_paragraph(); r = p.add_run(text); r.italic = italic
    if color: r.font.color.rgb = RGBColor.from_string(color)
    if size: r.font.size = Pt(size)
    return p

# ---- Header ----
p = para("ALZHEIMER'S DISEASE KNOWLEDGE GRAPH", color="5A6574", size=9)
doc.add_paragraph("AlzKG Explorer", style="Title")
p = doc.add_paragraph()
p.add_run("Seven AD-specific datasets are filtered, mapped onto shared identifiers, and joined into one graph anchored on Alzheimer's disease (")
r = p.add_run("MONDO:0004975"); r.font.name = "Consolas"
p.add_run("). The figure below shows its core: the most strongly supported genes, variants, subtypes and drugs.")

# ---- How the graph is built ----
doc.add_heading("How the graph is built", level=1)

doc.add_heading("1. Pull AD evidence from each source", level=2)
src = sorted(S["per_source"].items(), key=lambda kv: -kv[1])
table(["Source", "AD edges retrieved"], [(k, fmt(v) if v else "manual download (not in this build)") for k, v in src], [3.0, 3.0])

doc.add_heading("2. Harmonize identifiers", level=2)
para("Each source names the same thing differently. Every node is mapped to one canonical id.")
table(["Entity", "Mapping"], [
    ("Genes", "APOE · 348 · PA55 → NCBIGene:348"),
    ("Variants", "RS#:429358 · ClinVar · GWAS → rs429358"),
    ("AD", "MESH:D000544 · PA443319 · \"Probable Alzheimer's Disease\" → MONDO:0004975"),
    ("Drugs", "merged by normalized name"),
], [1.3, 5.2], mono_cols=(1,))

doc.add_heading("3. Restrict to the AD neighborhood", level=2)
para("Keep edges within two hops of an AD node: direct evidence plus the second hop that matters for repurposing (gene → drug).")
p = doc.add_paragraph(); r = p.add_run(fmt(S["ad_degree"])); r.bold = True; r.font.size = Pt(18)
p.add_run("  edges touch AD directly")

doc.add_heading("4. Indexed graph", level=2)
p = doc.add_paragraph()
for n, lab in ((S["final_nodes"], "nodes"), (S["final_edges"], "edges")):
    r = p.add_run(fmt(n)); r.bold = True; r.font.size = Pt(18); p.add_run(f"  {lab}     ")
table(["Node type", "Count"], [(TYPES.get(t, t), fmt(n)) for t, n in S["node_types"].items()], [3.0, 3.0])
doc.add_paragraph("Edges by source", style="Heading 2")
table(["Edge source", "Edges"], [(k, fmt(v)) for k, v in S["edge_sources"].items()], [3.0, 3.0])

# ---- Core subgraph ----
doc.add_heading("Core subgraph", level=1).paragraph_format.page_break_before = True
para(f"{S['subgraph_nodes']} of {fmt(S['final_nodes'])} nodes and {S['subgraph_links']} of {fmt(S['final_edges'])} edges, chosen for multi-source support and clinical relevance.")
doc.add_picture(graph_png, width=Inches(6.5))
doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
doc.paragraphs[-1].paragraph_format.keep_with_next = True
para("Figure 1. Core subgraph centered on Alzheimer's disease. Node colors: red = disease, teal = gene/protein, "
     "amber = variant, indigo = drug/compound. Thicker edges are asserted by more datasets.", italic=True, color="5A6574", size=9.5)
para("Edge colors are neutral; each link lists every dataset that asserts it (see Table B). "
     "DisGeNET and NIAGADS need manual downloads and are not in this build.", color="5A6574", size=10)

byId = {n["id"]: n for n in D["nodes"]}
deg = {n["id"]: 0 for n in D["nodes"]}
for l in D["links"]:
    deg[l["from"]] += 1; deg[l["to"]] += 1

doc.add_heading("Table A. Nodes in the core subgraph", level=1)
order = list(TYPES)
nodes = sorted(D["nodes"], key=lambda n: (order.index(n["type"]) if n["type"] in order else 99, -n["degree"]))
table(["Name", "Type", "Identifier", "Vocabulary", "Edges (full graph)"],
      [(n["name"], TYPES.get(n["type"], n["type"]), n["cid"], n["vocab"], fmt(n["degree"])) for n in nodes],
      [2.2, 1.1, 1.3, 1.1, 0.8], mono_cols=(2,), size=8.5)

doc.add_heading("Table B. Edges in the core subgraph", level=1)
links = sorted(D["links"], key=lambda l: (-len(l["sources"]), byId[l["from"]]["name"]))
table(["From", "To", "Relation", "Sources"],
      [(byId[l["from"]]["name"], byId[l["to"]]["name"], "; ".join(l["labels"]) or l["relation"], ", ".join(l["sources"])) for l in links],
      [2.0, 2.0, 1.3, 1.2], size=8)

doc.save(out)
tmpdir.cleanup()
print("saved", out)
