const fs = require('node:fs');
const path = require('node:path');

const root = __dirname;
const W = 595;
const H = 842;
const ink = [0.12, 0.17, 0.20];
const muted = [0.32, 0.39, 0.42];
const teal = [0.09, 0.42, 0.40];
const palette = [
  [0.10, 0.40, 0.72], [0.92, 0.39, 0.12], [0.10, 0.62, 0.34],
  [0.80, 0.18, 0.20], [0.49, 0.31, 0.68], [0.48, 0.34, 0.24]
];

function csvRows(file, header = false) {
  const lines = fs.readFileSync(path.join(root, file), 'utf8').trim().split(/\r?\n/);
  return lines.slice(header ? 1 : 0).map(line => line.split(',').map(Number));
}

function esc(value) {
  return String(value).replace(/[^\x20-\x7e]/g, '?').replace(/([\\()])/g, '\\$1');
}

function color(rgb) { return `${rgb.map(v => v.toFixed(3)).join(' ')} rg\n`; }
function strokeColor(rgb) { return `${rgb.map(v => v.toFixed(3)).join(' ')} RG\n`; }

class Page {
  constructor() { this.ops = ''; }
  text(x, y, value, size = 10, bold = false, rgb = ink) {
    this.ops += color(rgb) + `BT /${bold ? 'F2' : 'F1'} ${size} Tf 1 0 0 1 ${x.toFixed(2)} ${y.toFixed(2)} Tm (${esc(value)}) Tj ET\n`;
  }
  line(x1, y1, x2, y2, rgb = ink, width = 0.7) {
    this.ops += strokeColor(rgb) + `${width} w ${x1.toFixed(2)} ${y1.toFixed(2)} m ${x2.toFixed(2)} ${y2.toFixed(2)} l S\n`;
  }
  rect(x, y, w, h, rgb, fill = true) {
    this.ops += (fill ? color(rgb) : strokeColor(rgb)) + `${x} ${y} ${w} ${h} re ${fill ? 'f' : 'S'}\n`;
  }
  dot(x, y, size = 2.4, rgb = teal) {
    this.rect(x - size / 2, y - size / 2, size, size, rgb, true);
  }
  paragraph(value, x, y, width = 510, size = 10, lineHeight = 14, bold = false, rgb = ink) {
    const maxChars = Math.max(18, Math.floor(width / (size * 0.52)));
    const words = String(value).split(/\s+/);
    let current = '';
    for (const word of words) {
      const next = current ? `${current} ${word}` : word;
      if (next.length > maxChars && current) {
        this.text(x, y, current, size, bold, rgb);
        y -= lineHeight;
        current = word;
      } else current = next;
    }
    if (current) { this.text(x, y, current, size, bold, rgb); y -= lineHeight; }
    return y;
  }
  heading(value, y) {
    this.text(42, y, value, 16, true, teal);
    this.line(42, y - 6, 553, y - 6, [0.78, 0.83, 0.84], 0.65);
    return y - 22;
  }
  frame(pageNumber) {
    this.rect(0, 806, W, 36, [0.14, 0.24, 0.27]);
    this.text(42, 819, 'MULTIOBJECTIVE OPTIMIZATION  /  LAB 7', 9, true, [1, 1, 1]);
    this.line(42, 32, 553, 32, [0.78, 0.83, 0.84], 0.6);
    this.text(42, 18, 'Results report | source outputs under moo_lab7/results', 8, false, muted);
    this.text(535, 18, String(pageNumber), 8, false, muted);
  }
  chart(data, options) {
    const { x, y, w, h, xMin, xMax, yMin, yMax, xLabel, yLabel, logY = false } = options;
    const left = x + 38, right = x + w - 10, bottom = y + 28, top = y + h - 18;
    const px = value => left + (value - xMin) / (xMax - xMin) * (right - left);
    const py = value => {
      const a = logY ? Math.log10(Math.max(value, yMin)) : value;
      const lo = logY ? Math.log10(yMin) : yMin;
      const hi = logY ? Math.log10(yMax) : yMax;
      return bottom + (a - lo) / (hi - lo) * (top - bottom);
    };
    this.line(left, bottom, right, bottom, ink, 0.8);
    this.line(left, bottom, left, top, ink, 0.8);
    const xTicks = options.xTicks || 4;
    for (let i = 0; i <= xTicks; i++) {
      const v = xMin + (xMax - xMin) * i / xTicks;
      const xx = px(v);
      this.line(xx, bottom, xx, bottom - 3, ink, 0.6);
      this.text(xx - 8, bottom - 14, options.xFormat ? options.xFormat(v) : v.toFixed(1), 7, false, muted);
    }
    const yTickValues = options.yTickValues || Array.from({ length: 5 }, (_, i) => yMin + (yMax - yMin) * i / 4);
    for (const v of yTickValues) {
      const yy = py(v);
      this.line(left - 3, yy, left, yy, ink, 0.6);
      this.text(x + 1, yy - 2, options.yFormat ? options.yFormat(v) : v.toFixed(1), 7, false, muted);
    }
    this.text((left + right) / 2 - 18, y + 2, xLabel, 8, false, ink);
    this.text(x + 1, top + 6, yLabel, 8, false, ink);
    return { px, py, left, right, bottom, top };
  }
  seriesLine(rows, map, xi, yi, rgb, width = 1) {
    const points = rows.map(row => [map.px(row[xi]), map.py(row[yi])]);
    if (!points.length) return;
    this.ops += strokeColor(rgb) + `${width} w ${points[0][0].toFixed(2)} ${points[0][1].toFixed(2)} m `;
    for (const point of points.slice(1)) this.ops += `${point[0].toFixed(2)} ${point[1].toFixed(2)} l `;
    this.ops += 'S\n';
  }
}

function addTable(page, x, y, columns, rows, rowHeight = 39, fontSize = 8) {
  const totalWidth = columns.reduce((sum, col) => sum + col.width, 0);
  page.rect(x, y - 20, totalWidth, 20, [0.14, 0.24, 0.27]);
  let cx = x;
  for (const col of columns) {
    page.text(cx + 5, y - 13, col.title, 8, true, [1, 1, 1]);
    cx += col.width;
  }
  let rowY = y - 20;
  rows.forEach((row, index) => {
    page.rect(x, rowY - rowHeight, totalWidth, rowHeight, index % 2 ? [0.96, 0.97, 0.97] : [1, 1, 1]);
    page.line(x, rowY - rowHeight, x + totalWidth, rowY - rowHeight, [0.84, 0.87, 0.88], 0.45);
    cx = x;
    row.forEach((value, i) => {
      page.paragraph(value, cx + 5, rowY - 12, columns[i].width - 9, fontSize, 10, false, ink);
      cx += columns[i].width;
    });
    rowY -= rowHeight;
  });
  return rowY;
}

const pages = [];
function newPage(number) {
  const page = new Page();
  page.frame(number);
  pages.push(page);
  return page;
}

// Page 1: summary, metric conventions, and Task D comparison table.
{
  const p = newPage(1);
  p.text(42, 774, 'Study cards, comparisons, and trace explanations', 20, true, ink);
  p.text(42, 752, '37 catalogued methods | 10,000-evaluation budget where specified | fixed seeds', 9, false, muted);
  let y = p.heading('Executive summary', 720);
  y = p.paragraph('Performance depends on the problem structure. On P3, MOBA and PESA-II lead with nearly tied hypervolume and low IGD. On P5, many methods reach almost the same HV ceiling, so spacing and IGD are needed to distinguish their fronts. For DTLZ2 with five objectives, HypE has the highest HV while RVEA has lower IGD. On the 15-city P8 TSP, MOACO leads the compared algorithms by mean HV.', 42, y, 510, 10, 14);
  y = p.heading('How to read the metrics', y - 5);
  y = p.paragraph('All objectives are minimised internally. HV uses fixed per-problem normalisation and reference point 1.1; larger is better. IGD measures distance to a reference front and Schott spacing measures regularity; lower is better. Task A3 is one seeded run; Task D aggregates 10 seeds. P5 IGD uses the pooled non-dominated runs, not an analytic global front. Five-objective HV is estimated by Monte Carlo.', 42, y, 510, 9, 12);
  y = p.heading('Task D representative results', y - 5);
  addTable(p, 42, y, [
    { title: 'Problem', width: 100 }, { title: 'Top by HV', width: 108 },
    { title: 'HV mean +/- SD', width: 110 }, { title: 'IGD mean +/- SD', width: 110 },
    { title: 'Interpretation', width: 83 }
  ], [
    ['P3, 2 objectives', 'MOBA (29)', '0.53874 +/- 0.00003', '0.00379 +/- 0.00003', 'PESA-II close second'],
    ['P5, constrained', 'SMS-EMOA (33)', '1.19580 +/- 0.00000', '0.00509 +/- 0.00022', 'MOBA has lower IGD'],
    ['DTLZ2, M=5', 'HypE (34)', '1.2324 +/- 0.0066', '0.2622 +/- 0.0048', 'RVEA IGD is 0.1735'],
    ['P8, 15-city TSP', 'MOACO (25)', '0.7832 +/- 0.0219', '0.0108 +/- 0.0182', 'NSGA-II lowest spacing']
  ], 42, 7.3);
  p.paragraph('A rank by HV is not an all-metric win. On P5, SMS-EMOA ranks first in HV, while MOBA has lower mean spacing and IGD. Values are from taskD_compare.txt; SD is sample SD over 10 seeds.', 42, 270, 510, 8, 11, false, muted);
}

// Page 2: Task C hand traces and sensitivity.
{
  const p = newPage(2);
  let y = p.heading('Task C: trace explanations', 774);
  p.text(42, y, 'Algorithm 8 | TOPSIS on P12', 11, true, ink); y -= 17;
  y = p.paragraph('Vector-normalise each criterion, apply weights, and compute distances D+ and D- from the positive and negative ideal. Closeness C = D- / (D+ + D-); larger C ranks better. With weights (0.3, 0.3, 0.2, 0.2), S3 ranks first (0.6298), followed closely by S1 (0.6214), then S5, S2, and S4. Shifting 0.1 weight from criterion 2 to criterion 1 changes the order to S1, S2, S3, S5, S4. The chosen ranking depends on decision-maker preference.', 42, y, 510, 9, 12);
  y -= 5; p.text(42, y, 'Algorithm 32 | IBEA on the small P1 trace', 11, true, ink); y -= 17;
  y = p.paragraph('The additive epsilon indicator compares normalised objective vectors: Ieps(i,j) = max_k(fn_i,k - fn_j,k). IBEA combines pairwise indicators into exponential fitness. At kappa=0.05 the first trace fitness spans -2.6930 to approximately 0; environmental selection repeatedly removes the least-fit member from the parent-offspring union until four remain. Kappa controls the sharpness of the fitness response.', 42, y, 510, 9, 12);
  y = p.paragraph('P3 sensitivity (10 seeds): kappa=0.05 gives HV 0.5370 +/- 0.0001 and IGD 0.0085 +/- 0.0004. kappa=0.01 has slightly lower HV (0.5356 +/- 0.0077) but better mean IGD (0.0067 +/- 0.0027). kappa=0.25 is worse on both (HV 0.5165, IGD 0.0453). The preferred setting depends on metric; 0.05 is more stable in HV.', 42, y, 510, 9, 12);
  y -= 5; p.text(42, y, 'Algorithm 20 | MOEA/D-ACO on the hand-checkable TSP', 11, true, ink); y -= 17;
  y = p.paragraph('Each weight-vector subproblem constructs a tour using its pheromone matrix and a distance/time heuristic. A tour replaces neighbouring incumbents when its normalised Tchebycheff value improves; pheromone evaporates by rho and is reinforced along neighbour tours. In iteration 1 the two tours have objectives (15,25) and (13,21); in iteration 2 the reported tours are (13,21) and (13,19). This five-city trace is distinct from the 15-city Task D benchmark.', 42, y, 510, 9, 12);
  y = p.paragraph('P8 rho sensitivity: rho=0.02 gives the highest tested mean HV (0.7854 +/- 0.0113); rho=0.10 gives 0.7730 +/- 0.0232; rho=0.50 falls to 0.7203 +/- 0.0523. In this experiment, high evaporation weakens useful pheromone memory.', 42, y, 510, 9, 12);
  p.paragraph('Source: handtrace_taskC.txt and taskC_sensitivity.txt. Trace values illustrate algorithm steps; they are not the multi-seed benchmark means.', 42, y - 3, 510, 8, 11, false, muted);
}

// Page 3: P3 and P8 figures redrawn from stored fronts.
{
  const p = newPage(3);
  let y = p.heading('Figures: continuous and discrete trade-offs', 774);
  p.text(42, y, 'Figure 1. NSGA-II on P3 with analytic true front', 10, true, ink);
  const chart1 = p.chart([], { x: 45, y: 470, w: 500, h: 255, xMin: 0, xMax: 1, yMin: 0, yMax: 1, xLabel: 'f1', yLabel: 'f2', xTicks: 5, yTickValues: [0, 0.2, 0.4, 0.6, 0.8, 1], yFormat: v => v.toFixed(1) });
  const trueP3 = csvRows('results/true_P3.csv');
  const frontP3 = csvRows('results/fronts/alg14_P3.csv', true);
  p.seriesLine(trueP3, chart1, 0, 1, [0.08, 0.09, 0.10], 1.25);
  frontP3.forEach(row => p.dot(chart1.px(row[0]), chart1.py(row[1]), 2.6, palette[0]));
  p.line(385, 704, 408, 704, [0.08, 0.09, 0.10], 1.2); p.text(413, 701, 'analytic front', 7, false, muted);
  p.dot(396, 690, 3, palette[0]); p.text(413, 687, 'NSGA-II points', 7, false, muted);
  p.paragraph('The points closely follow the analytic curve, consistent with the single-run IGD of 0.00484. The visual overlap makes small density differences difficult to assess.', 45, 462, 500, 8, 11, false, muted);

  p.text(42, 431, 'Figure 2. MOEA/D-ACO on the P8 distance-time plane', 10, true, ink);
  const p8 = csvRows('results/fronts/alg20_P8.csv');
  const chart2 = p.chart([], { x: 45, y: 164, w: 500, h: 250, xMin: 300, xMax: 385, yMin: 350, yMax: 440, xLabel: 'total distance', yLabel: 'total time', xTicks: 5, yTicks: 4, xFormat: v => `${Math.round(v)}`, yFormat: v => `${Math.round(v)}` });
  p8.forEach(row => p.dot(chart2.px(row[0]), chart2.py(row[1]), 3, palette[0]));
  p.paragraph('Lower distance tends to accompany higher travel time. These are single-run tours; no reference front is drawn. Use the multi-seed table for algorithm comparisons.', 45, 155, 500, 8, 11, false, muted);
}

// Page 4: P13 trajectories, interpretation, constraints, and limitations.
{
  const p = newPage(4);
  let y = p.heading('Figure 3 and study interpretation', 774);
  p.text(42, y, 'P13(a) task losses over 300 iterations (log scale)', 10, true, ink);
  const trajectories = csvRows('results/P13_trajectories.csv', true);
  const chart = p.chart([], { x: 45, y: 435, w: 500, h: 300, xMin: 0, xMax: 300, yMin: 0.4, yMax: 100, logY: true, xLabel: 'iteration', yLabel: 'task loss', xTicks: 6, xFormat: v => `${Math.round(v)}`, yTickValues: [1, 10, 100], yFormat: v => `${v}` });
  const series = [
    ['fixed_L1', 'fixed L1'], ['fixed_L2', 'fixed L2'], ['mgda_L1', 'MGDA L1'],
    ['mgda_L2', 'MGDA L2'], ['gradnorm_L1', 'GradNorm L1'], ['gradnorm_L2', 'GradNorm L2']
  ];
  series.forEach(([field], i) => {
    const index = ['iter', 'fixed_L1', 'fixed_L2', 'mgda_L1', 'mgda_L2', 'gradnorm_L1', 'gradnorm_L2'].indexOf(field);
    p.seriesLine(trajectories, chart, 0, index, palette[i], 1.15);
    p.line(330, 712 - i * 13, 346, 712 - i * 13, palette[i], 1.4);
    p.text(351, 709 - i * 13, series[i][1], 7, false, muted);
  });
  p.paragraph('Fixed equal weights approach losses (18,18). MGDA sharply lowers L1 while leaving L2 near 60.5, selecting a different Pareto point; GradNorm balances differently and oscillates. The methods encode different trade-offs, not one universally best trajectory.', 45, 425, 500, 8, 11, false, muted);

  y = p.heading('Other results and limitations', 388);
  y = p.paragraph('P12: TOPSIS selects S3 while VIKOR at v=0.5 selects S1; weights and v alter rankings. P11: goal programming meets the profit goal but misses the 600-unit target by about 100 units. Lexicographic optimisation selects (300,200), reaching 500 units while preserving profit 2400. P14: the requested 1.8 kg / 12 h target lies behind the front; the best runtime at 1.8 kg is 14.04 h.', 42, y, 510, 8.5, 11);
  y = p.paragraph('Constraint handling matters: A3 plain NSGA-II on P5 returns 53% feasible points, while constrained NSGA-II returns 100% feasible points. A high HV does not certify feasibility.', 42, y, 510, 8.5, 11);
  y = p.paragraph('P9 was not run because its instructor dataset and 5-fold 5-NN evaluation are unavailable; P13(b) is not implemented. P8 is a seeded synthetic 15-city TSP. MOEA/D-ACO, MOBA, MOSA, and MO-TLBO are compact educational variants. Claims are limited by budgets, seeds, and reference-front availability.', 42, y, 510, 8.5, 11);
  p.paragraph('Viva prompts: Why can weighted sums miss non-convex fronts? Why can HV and IGD disagree? What changes when feasibility rules are disabled? How do TOPSIS weights and VIKOR v encode preference? How does rho affect pheromone memory?', 42, y - 2, 510, 8, 11, true, teal);
}

function makePdf(allPages) {
  const objects = [];
  const set = (id, value) => { objects[id] = value; };
  set(1, '<< /Type /Catalog /Pages 2 0 R >>');
  const pageIds = allPages.map((_, i) => 5 + i * 2);
  set(2, `<< /Type /Pages /Kids [${pageIds.map(id => `${id} 0 R`).join(' ')}] /Count ${allPages.length} >>`);
  set(3, '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>');
  set(4, '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>');
  allPages.forEach((page, index) => {
    const pageId = pageIds[index];
    const streamId = pageId + 1;
    const stream = page.ops;
    set(pageId, `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${W} ${H}] /Resources << /Font << /F1 3 0 R /F2 4 0 R >> >> /Contents ${streamId} 0 R >>`);
    set(streamId, `<< /Length ${Buffer.byteLength(stream, 'ascii')} >>\nstream\n${stream}endstream`);
  });
  let pdf = '%PDF-1.4\n% Report generated from CSV outputs\n';
  const offsets = [0];
  for (let id = 1; id < objects.length; id++) {
    offsets[id] = Buffer.byteLength(pdf, 'ascii');
    pdf += `${id} 0 obj\n${objects[id]}\nendobj\n`;
  }
  const xrefOffset = Buffer.byteLength(pdf, 'ascii');
  pdf += `xref\n0 ${objects.length}\n0000000000 65535 f \n`;
  for (let id = 1; id < objects.length; id++) pdf += `${String(offsets[id]).padStart(10, '0')} 00000 n \n`;
  pdf += `trailer\n<< /Size ${objects.length} /Root 1 0 R >>\nstartxref\n${xrefOffset}\n%%EOF\n`;
  return Buffer.from(pdf, 'ascii');
}

const output = path.join(root, 'REPORT.pdf');
fs.writeFileSync(output, makePdf(pages));
console.log(`Wrote ${output} (${pages.length} pages, ${fs.statSync(output).size} bytes)`);