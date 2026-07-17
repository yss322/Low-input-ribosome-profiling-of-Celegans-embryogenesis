---
title: About
toc: false
---

# About

An interactive tool for exploring translational regulation in the early *C. elegans* embryo (1- to 8-cell stage).

---

## Data Sources

| Dataset | Source | What it measures |
|---------|--------|-----------------|
| TE, Ribo, RNA | [Shukla et al. 2025, *Cell Reports*](https://www.cell.com/cell-reports/fulltext/S2211-1247(25)00549-2) | Ribosome occupancy & mRNA abundance per stage (CLR-normalized, batch-corrected) |
| Spatial mRNA | [Tintori et al. 2016, *Developmental Cell*](https://www.cell.com/developmental-cell/fulltext/S1534-5807(16)30518-4) | Per-cell transcript abundance at blastomere resolution |

**Gene set:** 4,905 genes passing CPM filters in both Ribo-seq and RNA-seq.

---

## How to Read the Plots

### Top Panels — TE, Ribosome Occupancy & mRNA

```js
display(html`<img src="${await FileAttachment("data/img/example_plots.png").url()}" style="width:100%; border:1px solid #ddd; border-radius:6px; margin:8px 0" />`)
```

**Left — Translational Efficiency (TE)**

- **TE = Ribosome CLR − mRNA CLR** (log-ratio of ribosome footprints to mRNA)
- Positive = translated more efficiently than average. Negative = translationally repressed.
- Grey dashed lines are reference genes: **lem-2** (high TE), **gpd-4** (housekeeping), **nos-2** (low TE, repressed)

**Right — Ribosome Occupancy & mRNA Abundance**

- Orange = ribosome occupancy. Blue = mRNA abundance.
- Lines diverging → translational regulation. Lines converging → mRNA-driven change.

### Bottom Panel — Spatial mRNA (Embryo Pictogram)

```js
display(html`<img src="${await FileAttachment("data/img/example_embryo.png").url()}" style="width:100%; border:1px solid #ddd; border-radius:6px; margin:8px 0" />`)
```

- Darker = more mRNA in that cell. Asymmetric coloring = mRNA is localized.
- Cells at each stage: P0 → AB/P1 → ABa/ABp/EMS/P2 → ABal/ABar/ABpl/ABpr/MS/E/C/P3
- This shows mRNA only, not where translation occurs.

---

## Gene Names

The TE data uses common names (`par-3`, `nos-2`). The Tintori data uses systematic names (`F54E7.3`). This tool maps between them automatically via WormBase aliases. If the spatial panel shows "not found," no alias exists for that gene.

---

## Citation

> Shukla et al. (2025). Landscape and regulation of mRNA translation in the early *C. elegans* embryo. *Cell Reports*, 44(4).

> Tintori et al. (2016). A Transcriptional Lineage of the Early *C. elegans* Embryo. *Developmental Cell*, 38(4).
