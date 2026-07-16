---
title: C. elegans Translational Efficiency Explorer
toc: false
---

<h1 style="white-space:nowrap"><em>C. elegans</em> — Translational Efficiency Explorer</h1>

<div class="grid grid-cols-2">
<div class="card" style="grid-column: span 2; padding: 12px 16px;">
<p style="margin:0; font-size:14px; max-width:100%">Select a gene to view its translational efficiency, ribosome occupancy, and mRNA abundance across early embryonic cell stages. TE and ribo/RNA from CLR-transformed, batch-corrected ribosome profiling (Shukla et al.). Spatial mRNA from Tintori et al. 2016. All panels: 1- to 8-cell stage.</p>
</div>
</div>

```js
const raw = await FileAttachment("data/averaged_clr_all.csv").csv({typed: true});
```

```js
const geneNames = raw.map(d => d.gene_name).sort();
const gene = view(Inputs.select(geneNames, {value: "par-3", label: "Gene"}));
```

```js
const stages = ["1-cell", "2-cell", "4-cell", "8-cell"];

function reshapeTE(row) {
  return stages.map((s, i) => ({
    gene: row.gene_name,
    stage: s,
    stageIdx: i,
    TE: [row.WT_one_cell_TE, row.WT_two_cell_TE, row.WT_four_cell_TE, row.WT_eight_cell_TE][i]
  }));
}

function reshapeRiboRNA(row) {
  const ribo = [row["one_cell.ribo"], row["two_cell.ribo"], row["four_cell.ribo"], row["eight_cell.ribo"]];
  const rna = [row["one_cell.rna"], row["two_cell.rna"], row["four_cell.rna"], row["eight_cell.rna"]];
  return stages.flatMap((s, i) => [
    {stage: s, stageIdx: i, value: ribo[i], type: "Ribosome occupancy"},
    {stage: s, stageIdx: i, value: rna[i], type: "mRNA abundance"}
  ]);
}

const refGeneNames = ["lem-2", "gpd-4", "nos-2"];
const refLabels = {"lem-2": "lem-2 (High TE)", "gpd-4": "gpd-4 (Housekeeping)", "nos-2": "nos-2 (Low TE)"};
const refRows = raw.filter(d => refGeneNames.includes(d.gene_name));
const refTE = refRows.flatMap(reshapeTE);

const selectedRow = raw.find(d => d.gene_name === gene);
const targetTE = selectedRow ? reshapeTE(selectedRow) : [];
const riboRnaData = selectedRow ? reshapeRiboRNA(selectedRow) : [];
```

<div class="grid grid-cols-2">
<div class="card">

```js
display(resize((width) => Plot.plot({
  width,
  height: 450,
  marginLeft: 140,
  marginRight: 80,
  x: {type: "point", domain: stages, label: null, padding: 0.2},
  y: {domain: [-5, 5], ticks: [-5, -2.5, 0, 2.5, 5], label: "Normalized TE"},
  marks: [
    Plot.ruleY([0], {stroke: "#ddd"}),
    Plot.lineY(refTE, {x: "stage", y: "TE", z: "gene", stroke: "#999", strokeWidth: 2.5, strokeDasharray: "10 6"}),
    Plot.dot(refTE, {x: "stage", y: "TE", stroke: "#999", fill: "white", r: 5, strokeWidth: 1.5}),
    Plot.text(refTE.filter(d => d.stage === "1-cell"), {x: "stage", y: "TE", text: d => refLabels[d.gene], dx: -25, textAnchor: "end", fontSize: 11, fill: "#444"}),
    Plot.lineY(targetTE, {x: "stage", y: "TE", z: "gene", stroke: "red", strokeWidth: 3}),
    Plot.dot(targetTE, {x: "stage", y: "TE", fill: "red", r: 5}),
    Plot.text(targetTE.filter(d => d.stage === "8-cell"), {x: "stage", y: "TE", text: "gene", dx: 18, textAnchor: "start", fontSize: 13, fill: "red", fontWeight: "bold"})
  ]
})))
```

</div>
<div class="card">

```js
// Offset labels if values are close
const endPoints = riboRnaData.filter(d => d.stage === "8-cell");
const labelData = endPoints.map(d => ({...d}));
if (labelData.length === 2 && Math.abs(labelData[0].value - labelData[1].value) < 0.8) {
  const mid = (labelData[0].value + labelData[1].value) / 2;
  labelData[0].value = mid + 0.5;
  labelData[1].value = mid - 0.5;
}

display(resize((width) => Plot.plot({
  width,
  height: 450,
  marginLeft: 60,
  marginRight: 150,
  x: {type: "point", domain: stages, label: null, padding: 0.2},
  y: {domain: [-5, 5], ticks: [-5, -2.5, 0, 2.5, 5], label: "Normalized Values"},
  color: {domain: ["Ribosome occupancy", "mRNA abundance"], range: ["rgb(228,88,10)", "rgb(55,135,192)"], legend: false},
  marks: [
    Plot.ruleY([0], {stroke: "#ddd"}),
    Plot.lineY(riboRnaData, {x: "stage", y: "value", stroke: "type", z: "type", strokeWidth: 2.5}),
    Plot.dot(riboRnaData, {x: "stage", y: "value", fill: "type", r: 5}),
    Plot.text(labelData, {x: "stage", y: "value", text: "type", dx: 18, textAnchor: "start", fontSize: 12, fill: "type"})
  ]
})))
```

</div>
</div>

## Spatial mRNA Expression (Tintori et al. 2016)

```js
const embryos = await FileAttachment("data/embryos.json").json();
const tintori = await FileAttachment("data/tintori_expression.csv").csv({typed: true});
const aliasRaw = await FileAttachment("data/tintori_aliases.csv").text();
```

```js
// Build alias lookup: any alias → primary Tintori ID (first column)
const aliasToId = new Map();
aliasRaw.split("\n").forEach(line => {
  const parts = line.split(",").map(s => s.trim()).filter(Boolean);
  if (parts.length === 0) return;
  const primary = parts[0];
  parts.forEach(name => aliasToId.set(name, primary));
});

const cellNameMap = embryos.cell_name_map;
const tintoriId = aliasToId.get(gene) || gene;
const tintoriRow = tintori.find(d => d[""] === tintoriId);
```

```js
function renderEmbryoPictogram(geneRow, maxVal) {
  if (!geneRow) {
    const p = document.createElement("p");
    p.style.color = "#888";
    p.textContent = "Gene not found in Tintori et al. dataset";
    return p;
  }

  const colorScale = d3.scaleSequential(t => d3.interpolateBlues(t * 0.85 + 0.15)).domain([0, maxVal]);
  const stageLabels = ["1-cell", "2-cell", "4-cell", "8-cell"];

  const container = document.createElement("div");
  container.style.cssText = "display:flex; flex-direction:column; gap:12px; padding:10px 0";

  // Legend — horizontal gradient bar
  const legendWrap = document.createElement("div");
  legendWrap.style.cssText = "display:flex; align-items:center; gap:8px; font-size:12px; color:#555";
  const legendLabel = document.createElement("span");
  legendLabel.textContent = "TPM:";
  legendWrap.append(legendLabel);

  const legendBar = document.createElement("div");
  legendBar.style.cssText = "display:flex; height:16px; border:1px solid #ccc; border-radius:3px; overflow:hidden";
  const numSteps = 40;
  for (let i = 0; i < numSteps; i++) {
    const v = (maxVal / numSteps) * i;
    const seg = document.createElement("div");
    seg.style.cssText = `width:4px;height:100%;background:${colorScale(v)}`;
    legendBar.append(seg);
  }
  legendWrap.append(legendBar);

  const legendMin = document.createElement("span");
  legendMin.textContent = "0";
  const legendMax = document.createElement("span");
  legendMax.textContent = Math.round(maxVal);
  // Insert min before bar, max after
  legendWrap.insertBefore(legendMin, legendBar);
  legendWrap.append(legendMax);
  container.append(legendWrap);

  // Embryo stages row
  const stagesRow = document.createElement("div");
  stagesRow.style.cssText = "display:flex; align-items:flex-end; gap:24px; justify-content:center";

  // Embryo stages
  embryos.stages.slice(0, 4).forEach((stage, idx) => {
    const vb = stage.svg_attrs.viewBox;
    const [, , w, h] = vb.split(" ").map(Number);

    const wrapper = document.createElement("div");
    wrapper.style.textAlign = "center";

    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("viewBox", vb);
    svg.setAttribute("width", w * 1.1);
    svg.setAttribute("height", h * 1.1);
    svg.style.overflow = "visible";

    const g1 = document.createElementNS("http://www.w3.org/2000/svg", "g");
    if (stage.g1_attrs?.transform) g1.setAttribute("transform", stage.g1_attrs.transform);
    const g2 = document.createElementNS("http://www.w3.org/2000/svg", "g");
    if (stage.g2_attrs?.transform) g2.setAttribute("transform", stage.g2_attrs.transform);

    stage.cells.forEach(cell => {
      const colName = cellNameMap[cell.cell_name] || cell.cell_name;
      const val = geneRow[colName] || 0;
      const fill = colorScale(val);

      let el;
      if (cell.type === "ellipse") {
        el = document.createElementNS("http://www.w3.org/2000/svg", "ellipse");
        el.setAttribute("cx", cell.attrs.cx);
        el.setAttribute("cy", cell.attrs.cy);
        el.setAttribute("rx", cell.attrs.rx);
        el.setAttribute("ry", cell.attrs.ry);
      } else if (cell.type === "path") {
        el = document.createElementNS("http://www.w3.org/2000/svg", "path");
        el.setAttribute("d", cell.attrs.d);
      }
      if (el) {
        el.setAttribute("fill", fill);
        el.setAttribute("stroke", "#555");
        el.setAttribute("stroke-width", "0.6");
        g2.append(el);
      }
    });

    g1.append(g2);
    svg.append(g1);
    wrapper.append(svg);

    const label = document.createElement("div");
    label.style.cssText = "font-size:12px; color:#555; margin-top:6px; font-weight:500";
    label.textContent = stageLabels[idx];
    wrapper.append(label);

    stagesRow.append(wrapper);
  });

  container.append(stagesRow);
  return container;
}
```

```js
const perCellCols = ["P0","AB","P1","ABa","ABp","EMS","P2","ABal","ABar","ABpl","ABpr","MS","E","C","P3"];
const maxExpr = tintoriRow ? Math.max(...perCellCols.map(c => tintoriRow[c] || 0), 1) : 500;

const card = document.createElement("div");
card.className = "card";
const heading = document.createElement("h3");
heading.style.cssText = "margin:0 0 8px";
heading.innerHTML = `<strong>${gene}</strong> — mRNA abundance per cell <span style="font-weight:normal;color:#888">(Tintori et al. 2016)</span>`;
card.append(heading);
card.append(renderEmbryoPictogram(tintoriRow, maxExpr));
display(card);
```

---

