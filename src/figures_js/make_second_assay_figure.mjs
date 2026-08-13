import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createCanvas } from "@napi-rs/canvas";
import sharp from "sharp";

const here = path.dirname(fileURLToPath(import.meta.url));
const csvPath = path.join(here, "outputs", "critique2_second_assay_operating.csv");
const outPath = path.join(here, "figures", "fig_second_assay_redesigned.png");

function parseCsv(text) {
  const [header, ...lines] = text.trim().split(/\r?\n/);
  const keys = header.split(",");
  return lines.map((line) => Object.fromEntries(line.split(",").map((value, i) => [keys[i], value])));
}

const rows = parseCsv(await fs.readFile(csvPath, "utf8"));
const datasets = [
  { cv: 0, label: "0% assay CV", dashed: false },
  { cv: 20, label: "20% assay CV", dashed: true },
].map((spec) => ({
  ...spec,
  values: rows
    .filter((row) => Number(row.assay_cv_pct) === spec.cv)
    .map((row) => ({
      rho: Number(row.rho),
      accuracy: Number(row.accuracy_pct),
      sensitivity: Number(row.sensitivity_pct),
      specificity: Number(row.specificity_pct),
      ppv: Number(row.ppv_pct),
      npv: Number(row.npv_pct),
      falseReassurance: Number(row.false_reassurance_pct),
    }))
    .sort((a, b) => a.rho - b.rho),
}));

const width = 2100;
const height = 1240;
const canvas = createCanvas(width, height);
const ctx = canvas.getContext("2d");
ctx.fillStyle = "#FFFFFF";
ctx.fillRect(0, 0, width, height);
ctx.textBaseline = "middle";

const ink = "#18212B";
const grid = "#DDE3EA";
const rule = "#7C8794";
const colors = {
  accuracy: "#205B9A",
  sensitivity: "#3E8D74",
  specificity: "#8A5DA8",
  ppv: "#205B9A",
  npv: "#D07035",
  falseReassurance: "#BE3D3D",
};
const margins = { left: 130, right: 106, top: 190, bottom: 248 };
const gap = 126;
const panelWidth = (width - margins.left - margins.right - gap) / 2;
const panelHeight = height - margins.top - margins.bottom;

function line(points, color, dashed, marker = "circle") {
  ctx.strokeStyle = color;
  ctx.lineWidth = 4;
  ctx.save();
  if (dashed) ctx.setLineDash([10, 8]);
  ctx.beginPath();
  points.forEach((point, index) => {
    if (index === 0) ctx.moveTo(point.x, point.y); else ctx.lineTo(point.x, point.y);
  });
  ctx.stroke();
  ctx.restore();
  ctx.fillStyle = color;
  for (const point of points) {
    if (marker === "square") ctx.fillRect(point.x - 5.5, point.y - 5.5, 11, 11);
    else if (marker === "triangle") {
      ctx.beginPath();
      ctx.moveTo(point.x, point.y - 7);
      ctx.lineTo(point.x - 6, point.y + 5.5);
      ctx.lineTo(point.x + 6, point.y + 5.5);
      ctx.closePath();
      ctx.fill();
    } else {
      ctx.beginPath(); ctx.arc(point.x, point.y, 6, 0, Math.PI * 2); ctx.fill();
    }
  }
}

function swatch(x, y, color, text, marker = "circle") {
  line([{ x, y }, { x: x + 38, y }], color, false, marker);
  ctx.fillStyle = ink;
  ctx.font = "500 25px Arial";
  ctx.textAlign = "left";
  ctx.fillText(text, x + 53, y);
}

function drawAxes({ x, y, minY, maxY, ticks, leftLabel, title, panelLetter, rightAxis }) {
  const base = y + panelHeight;
  const sy = (value) => base - ((value - minY) / (maxY - minY)) * panelHeight;
  const sx = (value) => x + (value / 0.94) * panelWidth;
  const xTicks = [0, 0.5, 0.75, 0.94];
  const xLabels = ["0", "0.50", "0.75", "0.94"];

  ctx.fillStyle = ink;
  ctx.font = "600 34px Arial";
  ctx.textAlign = "left";
  ctx.fillText(`${panelLetter}. ${title}`, x, y - 58);
  for (const tick of ticks) {
    const yy = sy(tick);
    ctx.strokeStyle = grid;
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x, yy); ctx.lineTo(x + panelWidth, yy); ctx.stroke();
    ctx.fillStyle = ink;
    ctx.font = "25px Arial";
    ctx.textAlign = "right";
    ctx.fillText(String(tick), x - 14, yy);
  }
  for (const tick of xTicks) {
    const xx = sx(tick);
    ctx.strokeStyle = grid;
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(xx, y); ctx.lineTo(xx, base); ctx.stroke();
    ctx.fillStyle = ink;
    ctx.font = "25px Arial";
    ctx.textAlign = "center";
    ctx.fillText(xLabels[xTicks.indexOf(tick)], xx, base + 36);
  }
  ctx.strokeStyle = rule;
  ctx.lineWidth = 1.8;
  ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x, base); ctx.lineTo(x + panelWidth, base); ctx.stroke();
  ctx.save();
  ctx.translate(x - 91, y + panelHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.font = "28px Arial";
  ctx.textAlign = "center";
  ctx.fillStyle = ink;
  ctx.fillText(leftLabel, 0, 0);
  ctx.restore();
  ctx.font = "28px Arial";
  ctx.textAlign = "center";
  ctx.fillStyle = ink;
  ctx.fillText("Assumed clearance correlation (\u03c1)", x + panelWidth / 2, base + 87);

  ctx.save();
  ctx.strokeStyle = "#89939D";
  ctx.lineWidth = 1.5;
  ctx.setLineDash([4, 6]);
  ctx.beginPath(); ctx.moveTo(sx(0.94), y); ctx.lineTo(sx(0.94), base); ctx.stroke();
  ctx.restore();
  ctx.fillStyle = "#56616D";
  ctx.font = "22px Arial";
  ctx.textAlign = "right";
  ctx.fillText("source-model \u03c1", sx(0.94) - 7, y + 22);

  if (rightAxis) {
    const rightSy = (value) => base - (value / 15) * panelHeight;
    for (const tick of [0, 5, 10, 15]) {
      const yy = rightSy(tick);
      ctx.strokeStyle = colors.falseReassurance;
      ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(x + panelWidth, yy); ctx.lineTo(x + panelWidth + 8, yy); ctx.stroke();
      ctx.fillStyle = colors.falseReassurance;
      ctx.font = "24px Arial";
      ctx.textAlign = "left";
      ctx.fillText(String(tick), x + panelWidth + 16, yy);
    }
    ctx.save();
    ctx.translate(x + panelWidth + 82, y + panelHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.font = "27px Arial";
    ctx.textAlign = "center";
    ctx.fillStyle = colors.falseReassurance;
    ctx.fillText("False reassurance (%)", 0, 0);
    ctx.restore();
    return { sx, sy, rightSy };
  }
  return { sx, sy };
}

function drawGlobalStyleLegend() {
  const y = 62;
  ctx.fillStyle = ink;
  ctx.font = "500 27px Arial";
  ctx.textAlign = "left";
  ctx.fillText("Simulated assay error:", margins.left, y);
  let x = margins.left + ctx.measureText("Simulated assay error:").width + 36;
  ctx.strokeStyle = "#4B5563";
  ctx.lineWidth = 4;
  ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + 45, y); ctx.stroke();
  ctx.fillText("0% CV", x + 60, y);
  x += 180;
  ctx.save(); ctx.setLineDash([10, 8]); ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + 45, y); ctx.stroke(); ctx.restore();
  ctx.fillText("20% CV", x + 60, y);
}

drawGlobalStyleLegend();

const a = drawAxes({
  x: margins.left,
  y: margins.top,
  minY: 0,
  maxY: 100,
  ticks: [0, 20, 40, 60, 80, 100],
  leftLabel: "Classification performance (%)",
  title: "Classification performance",
  panelLetter: "A",
});
for (const data of datasets) {
  line(data.values.map((row) => ({ x: a.sx(row.rho), y: a.sy(row.accuracy) })), colors.accuracy, data.dashed, "circle");
  line(data.values.map((row) => ({ x: a.sx(row.rho), y: a.sy(row.sensitivity) })), colors.sensitivity, data.dashed, "square");
  line(data.values.map((row) => ({ x: a.sx(row.rho), y: a.sy(row.specificity) })), colors.specificity, data.dashed, "triangle");
}
swatch(margins.left + 18, height - 76, colors.accuracy, "Accuracy", "circle");
swatch(margins.left + 222, height - 76, colors.sensitivity, "Sensitivity", "square");
swatch(margins.left + 465, height - 76, colors.specificity, "Specificity", "triangle");

const bx = margins.left + panelWidth + gap;
const b = drawAxes({
  x: bx,
  y: margins.top,
  minY: 50,
  maxY: 100,
  ticks: [50, 60, 70, 80, 90, 100],
  leftLabel: "Predictive value (%)",
  title: "Predictive values and false reassurance",
  panelLetter: "B",
  rightAxis: true,
});
for (const data of datasets) {
  line(data.values.map((row) => ({ x: b.sx(row.rho), y: b.sy(row.ppv) })), colors.ppv, data.dashed, "circle");
  line(data.values.map((row) => ({ x: b.sx(row.rho), y: b.sy(row.npv) })), colors.npv, data.dashed, "square");
  line(data.values.map((row) => ({ x: b.sx(row.rho), y: b.rightSy(row.falseReassurance) })), colors.falseReassurance, data.dashed, "triangle");
}
swatch(bx + 18, height - 105, colors.ppv, "Positive predictive value", "circle");
swatch(bx + 350, height - 105, colors.npv, "Negative predictive value", "square");
swatch(bx + 18, height - 55, colors.falseReassurance, "False reassurance (right axis)", "triangle");

await sharp(await canvas.encode("png"))
  .withMetadata({ density: 300 })
  .png()
  .toFile(outPath);
console.log(outPath);
