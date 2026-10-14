"""Fig. SX y Table SX (a: sitios; b: resumen por sitio, altitud y campaña) — Registros de Humicolopsis en la rizósfera de Nothofagus pumilio en Chile.
Fuente: Saona et al. 2025, ISME Commun 5:ycaf010 (Table S4: ASV × muestra; Table S1: coordenadas y altitudes).
  (a) mapa de los cuatro sitios (treeline y 200 m bajo el treeline), con detección de Humicolopsis
  (b) frecuencia de detección y abundancia relativa media por sitio, altitud y campaña de muestreo
Uso: python figSX_humicolopsis_chile.py"""
import json, math
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap, Normalize

plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
HUM = "#2a78d6"                                                            # Helotiales (paleta de COL en treeplot.py)
SEQ = LinearSegmentedColormap.from_list("seq", ["#dbe9f9", "#2a78d6", "#0d3b75"])   # secuencial de un solo tono

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "07_extra"
S4 = ROOT / "241213_tabla_s4_ycaf010.xlsx"
GEO = OUT / "datos/ne_50m_southern_cone.geojson"


def dms(d, m, s): return -(d + m / 60 + s / 3600)


# Table S1 de Saona et al. 2025 (coordenadas GPS y altitud, m s.n.m.)
SITES = pd.DataFrame([
    ("Chillán", "Nevados de Chillán", "Treeline", dms(36, 54, 22.086), dms(71, 23, 38.184), 2094),
    ("Chillán", "Nevados de Chillán", "Below treeline", dms(36, 54, 19.724), dms(71, 24, 24.602), 1821),
    ("Villarrica", "Villarrica National Park", "Treeline", dms(39, 23, 36.931), dms(71, 57, 50.767), 1399),
    ("Villarrica", "Villarrica National Park", "Below treeline", dms(39, 22, 52.41), dms(71, 58, 14.36), 1234),
    ("Antillanca", "Puyehue National Park (Antillanca)", "Treeline", dms(40, 47, 10.774), dms(72, 11, 33.644), 1332),
    ("Antillanca", "Puyehue National Park (Antillanca)", "Below treeline", dms(40, 46, 37.531), dms(72, 11, 55.366), 1132),
    ("Coyhaique", "Coyhaique National Reserve", "Treeline", dms(45, 30, 58.525), dms(72, 2, 38.72), 1211),
    ("Coyhaique", "Coyhaique National Reserve", "Below treeline", dms(45, 31, 14.632), dms(72, 3, 3.207), 998),
], columns=["Location", "Site", "Altitude", "Latitude", "Longitude", "Elevation_m"])
REGION = {"Chillán": "Ñuble", "Villarrica": "La Araucanía", "Antillanca": "Los Lagos", "Coyhaique": "Aysén"}
SEASONS = [("Fall", "Apr 2022"), ("Winter", "Jun 2022"), ("Spring", "Oct 2022"), ("Summer", "Jan 2023")]
DATE = dict(SEASONS)
LOCS = ["Chillán", "Villarrica", "Antillanca", "Coyhaique"]
ALTS = ["Treeline", "Below treeline"]

# ---------- datos ----------
d = pd.read_excel(S4, header=1)
tot = d.groupby("Sample_ID")["Abundance"].sum()
meta = d.drop_duplicates("Sample_ID").set_index("Sample_ID")[["Name", "Location", "Altitude", "Season", "Replicate"]]
h = d[d["Genus"] == "Humicolopsis"].copy()
h["Total_reads"] = h["Sample_ID"].map(tot)
h["RA"] = 100 * h["Abundance"] / h["Total_reads"]

# nivel género (suma de ASVs) por muestra
gen = h.groupby("Sample_ID").agg(reads=("Abundance", "sum")).join(meta)
gen["Total_reads"] = gen.index.map(tot)
gen["RA"] = 100 * gen["reads"] / gen["Total_reads"]
gen["pos"] = gen["reads"] > 0

summ = (gen.groupby(["Location", "Altitude", "Season"])
        .agg(n=("pos", "size"), pos=("pos", "sum"), reads=("reads", "sum"), mean_RA=("RA", "mean"), max_RA=("RA", "max"))
        .reset_index())
for a in ("ASV18", "ASV3239"):
    x = h[h["OTU_ID"] == a].assign(p=lambda t: t["Abundance"] > 0).groupby(["Location", "Altitude", "Season"])["p"].sum()
    summ = summ.merge(x.rename(f"pos_{a}").reset_index(), on=["Location", "Altitude", "Season"])

# ---------- Table SX ----------
tA = SITES.copy()
tA.insert(1, "Region", tA["Location"].map(REGION))
tA["Latitude"] = tA["Latitude"].round(5); tA["Longitude"] = tA["Longitude"].round(5)
st = gen.groupby(["Location", "Altitude"]).agg(Samples=("pos", "size"), Positive=("pos", "sum"),
                                               Reads=("reads", "sum")).reset_index()
tA = tA.merge(st, on=["Location", "Altitude"])
tA["Detection_%"] = (100 * tA["Positive"] / tA["Samples"]).round(1)
tA["Seasons_detected"] = tA.apply(lambda r: ", ".join(f"{s} ({DATE[s]})" for s, _ in SEASONS
    if ((gen.Location == r.Location) & (gen.Altitude == r.Altitude) & (gen.Season == s) & gen.pos).any()) or "Not detected", axis=1)
tA = tA.rename(columns={"Latitude": "Latitude_DD", "Longitude": "Longitude_DD", "Samples": "Samples_sequenced",
                        "Positive": "Samples_with_Humicolopsis", "Reads": "Humicolopsis_reads"})

tB = summ.copy()
tB["Sampling_date"] = tB["Season"].map(DATE)
tB["Detection_%"] = (100 * tB["pos"] / tB["n"]).round(1)
tB["mean_RA"] = tB["mean_RA"].round(3); tB["max_RA"] = tB["max_RA"].round(3)
tB["Location"] = pd.Categorical(tB["Location"], LOCS); tB["Altitude"] = pd.Categorical(tB["Altitude"], ALTS)
tB["Season"] = pd.Categorical(tB["Season"], [s for s, _ in SEASONS])
tB = tB.sort_values(["Location", "Altitude", "Season"])[
    ["Location", "Altitude", "Season", "Sampling_date", "n", "pos", "pos_ASV18", "pos_ASV3239", "Detection_%", "reads", "mean_RA", "max_RA"]]
tB.columns = ["Location", "Altitude", "Season", "Sampling_date", "Samples_sequenced", "Samples_with_Humicolopsis",
              "Samples_with_ASV18", "Samples_with_ASV3239", "Detection_%", "Humicolopsis_reads",
              "Mean_relative_abundance_%", "Max_relative_abundance_%"]

xy = SITES.set_index(["Location", "Altitude"])
with pd.ExcelWriter(OUT / "Table_SX_Humicolopsis_Chile.xlsx") as w:
    tA.to_excel(w, sheet_name="a_sites", index=False)
    tB.to_excel(w, sheet_name="b_site_season_summary", index=False)
tA.to_csv(OUT / "Table_SXa_sites.tsv", sep="\t", index=False)
tB.to_csv(OUT / "Table_SXb_summary.tsv", sep="\t", index=False)

# ---------- figura ----------
fig = plt.figure(figsize=(7.2, 4.3))
ax = fig.add_axes([0.03, 0.06, 0.30, 0.86])
geo = json.load(open(GEO))
for f in geo["features"]:
    g = f["geometry"]; polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
    chile = f["properties"]["ADMIN"] == "Chile"
    for p in polys:
        ax.add_patch(Polygon(p[0], closed=True, fc="#f3f2ee" if chile else "#ffffff", ec=MUTED if chile else "#c9c8c3",
                             lw=0.5, zorder=1))
ax.set_xlim(-76.5, -66.5); ax.set_ylim(-47.5, -33.5)
ax.set_aspect(1 / math.cos(math.radians(40.5)))
ax.set_xticks([-75, -72, -69]); ax.set_xticklabels(["75°W", "72°W", "69°W"])
ax.set_yticks([-46, -43, -40, -37, -34]); ax.set_yticklabels(["46°S", "43°S", "40°S", "37°S", "34°S"])
ax.tick_params(labelsize=5.8, length=2, color=MUTED)
for s in ax.spines.values(): s.set_color(MUTED)
ax.text(-74.3, -36.2, "Pacific\nOcean", fontsize=5.6, color=MUTED, style="italic", ha="center")
ax.text(-68.3, -42.3, "ARGENTINA", fontsize=5.6, color=MUTED, ha="center")
ax.text(-70.4, -34.6, "CHILE", fontsize=5.6, color=INK2, ha="center")
site = gen.groupby("Location").agg(pos=("pos", "sum"), n=("pos", "size")).reindex(LOCS)
for loc in LOCS:
    la, lo = SITES[SITES.Location == loc][["Latitude", "Longitude"]].mean()
    p, n = site.loc[loc, "pos"], site.loc[loc, "n"]
    ax.scatter(lo, la, s=40, marker="o", fc=HUM if p else "white", ec=INK if p else INK2, lw=0.8, zorder=3)
    ax.text(lo + 0.55, la, f"{loc}\n{p}/{n} samples", fontsize=5.8, va="center", color=INK, zorder=4)
ax.scatter([], [], s=30, fc=HUM, ec=INK, lw=0.8, label="Detected")
ax.scatter([], [], s=30, fc="white", ec=INK2, lw=0.8, label="Not detected")
ax.legend(loc="lower right", fontsize=5.6, frameon=False, handletextpad=0.3, borderaxespad=0.3)
fig.text(0.01, 0.955, "a", fontsize=11, fontweight="bold")

# (b) matriz sitio × campaña: tamaño = abundancia relativa media; color = % de muestras con detección
bx = fig.add_axes([0.52, 0.17, 0.40, 0.72])
rows = [(l, a) for l in LOCS for a in ALTS]
norm = Normalize(0, 100)
smax = summ["mean_RA"].max()
for i, (l, a) in enumerate(rows):
    for j, (s, _) in enumerate(SEASONS):
        r = summ[(summ.Location == l) & (summ.Altitude == a) & (summ.Season == s)]
        if r.empty:
            bx.text(j, i, "not\nsampled", ha="center", va="center", fontsize=4.8, color=MUTED); continue
        r = r.iloc[0]
        if r.pos == 0:
            bx.scatter(j, i, s=6, marker="o", fc="white", ec=MUTED, lw=0.6)
        else:
            bx.scatter(j, i, s=12 + 260 * r.mean_RA / smax, fc=SEQ(norm(100 * r.pos / r.n)), ec="white", lw=1.0, zorder=3)
        bx.text(j + 0.40, i + 0.30, f"{r.pos}/{r.n}", ha="right", va="center", fontsize=4.8, color=INK2)
bx.set_xticks(range(4)); bx.set_xticklabels([f"{s}\n{dt}" for s, dt in SEASONS], fontsize=5.8)
bx.xaxis.tick_top()
bx.set_yticks(range(len(rows))); bx.set_yticklabels([("T" if a == "Treeline" else "BT") + f"  ({int(xy.loc[(l, a), 'Elevation_m'])} m)"
                                                    for l, a in rows], fontsize=5.6)
for k, l in enumerate(LOCS):
    bx.annotate(l, xy=(0, 2 * k + 0.5), xycoords=("axes fraction", "data"), xytext=(-62, 0), textcoords="offset points",
                rotation=90, ha="center", va="center", fontsize=6.2, fontweight="bold")
    if k: bx.axhline(2 * k - 0.5, color=GRID, lw=0.8)
bx.set_xlim(-0.5, 3.5); bx.set_ylim(len(rows) - 0.5, -0.5)
bx.tick_params(length=0, pad=2)
for s in bx.spines.values(): s.set_visible(False)
fig.text(0.37, 0.955, "b", fontsize=11, fontweight="bold")

# leyendas de (b)
cax = fig.add_axes([0.925, 0.50, 0.012, 0.30])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SEQ), cax=cax)
cb.set_label("Samples with\n$\\it{Humicolopsis}$ (%)", fontsize=5.6); cb.ax.tick_params(labelsize=5.4, length=2); cb.outline.set_linewidth(0.4)
lx = fig.add_axes([0.915, 0.17, 0.08, 0.26]); lx.axis("off")
for k, v in enumerate([0.5, 1.5, 3.0]):
    lx.scatter(0.25, 2 - k, s=12 + 260 * v / smax, fc="#9fc1ea", ec="white", lw=1.0)
    lx.text(0.62, 2 - k, f"{v:g}%", va="center", fontsize=5.4, color=INK2)
lx.set_xlim(0, 1.2); lx.set_ylim(-0.6, 3.0)
lx.text(0.0, 2.75, "Mean relative\nabundance", fontsize=5.6, va="center")
fig.text(0.50, 0.06, "T, treeline; BT, ~200 m below treeline. Numbers: samples with $\\it{Humicolopsis}$ / samples sequenced.",
         fontsize=5.4, color=INK2)

for ext in ("pdf", "png"):
    fig.savefig(OUT / f"FigSX_Humicolopsis_Chile.{ext}", dpi=400 if ext == "png" else None, bbox_inches="tight")

# cifras para el texto
asv = h[h.Abundance > 0].groupby("OTU_ID")["Sample_ID"].nunique()
print("muestras secuenciadas:", gen.shape[0], "| con Humicolopsis:", int(gen.pos.sum()), "| por ASV:", asv.to_dict())
print("reads Humicolopsis:", int(gen.reads.sum()), "de", int(tot.sum()))
print(tA[["Location", "Altitude", "Elevation_m", "Samples_sequenced", "Samples_with_Humicolopsis", "Detection_%", "Seasons_detected"]].to_string())
print(gen[gen.pos].sort_values("RA", ascending=False).head(3)[["Location", "Altitude", "Season", "RA"]])
print("RA media en positivas:", round(gen[gen.pos].RA.mean(), 2), "| mediana:", round(gen[gen.pos].RA.median(), 2))
for c in ("Location", "Altitude", "Season"):
    print(gen.groupby(c).agg(pos=("pos", "sum"), n=("pos", "size"), mRA=("RA", "mean")).round(2).to_string())
