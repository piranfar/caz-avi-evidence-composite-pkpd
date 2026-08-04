import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createCanvas } from "@napi-rs/canvas";
import sharp from "sharp";

const here = path.dirname(fileURLToPath(import.meta.url));
const csvPath = path.join(here, "outputs", "prescriptive_decision_grid.csv");
const outPath = path.join(here, "figures", "fig_individualised_dose_redesigned.png");

function parseCsv(text) {
  const [header, ...lines] = text.trim().split(/\r?\n/);
  const keys = header.split(",");
  return lines.map((line) => Object.fromEntries(line.split(",").map((value, i) => [keys[i], value])));
}

const rows = parseCsv(await fs.readFile(csvPath, "utf8"));
const classes = [
  { key: "0\u201330", label: "0\u201330", color: "#1F4E85" },
  { key: "31\u201360", label: "31\u201360", color: "#28758C" },
  { key: "61\u201390", label: "61\u201390", color: "#4A987E" },
  { key: "91\u2013120", label: "91\u2013120", color: "#C86438" },
  { key: "121\u2013150", label: "121\u2013150", color: "#8067A5" },
];
const grouped = new Map(classes.map((item) => [item.key, []]));
for (const row of rows) {
  if (grouped.has(row.ekfc_class)) grouped.get(row.ekfc_class).push({
    mic: Number(row.mic_mg_l),
    dose: Number(row.median_placing_dose_g_day),
    within: Number(row.placing_dose_within_licensed_pct),
  });
}
for (const values of grouped.values()) values.sort((a, b) => a.mic - b.mic);

const width = 2100;
const height = 1120;
const canvas = createCanvas(width, height);
const ctx = canvas.getContext("2d");
ctx.fillStyle = "#FFFFFF";
ctx.fillRect(0, 0, width, height);
ctx.textBaseline = "middle";

const ink = "#18212B";
const grid = "#DDE3EA";
const rule = "#7C8794";
const alert = "#BE3D3D";
const margins = { left: 132, right: 52, top: 168, bottom: 128 };
const gap = 102;
const panelWidth = (width - margins.left - margins.right - gap) / 2;
const panelHeight = height - margins.top - margins.bottom;

function drawLegend() {
  const y = 64;
  ctx.fillStyle = ink;
  ctx.font = "500 27px Arial";
  ctx.textAlign = "left";
  ctx.fillText("EKFC class (mL/min/1.73 m\u00b2)", margins.left, y);
  let x = margins.left + ctx.measureText("EKFC class (mL/min/1.73 m\u00b2)").width + 40;
  for (const item of classes) {
    ctx.strokeStyle = item.color;
    ctx.lineWidth = 5;
    ctx.beginPath();
    ctx.moveTo(x, y);
    ctx.lineTo(x + 34, y);
    ctx.stroke();
    ctx.fillStyle = item.color;
    ctx.beginPath();
    ctx.arc(x + 17, y, 6.5, 0, Math.PI * 2);
    ctx.fill();
    ctx.fillStyle = ink;
    ctx.font = "500 27px Arial";
    ctx.fillText(item.label, x + 48, y);
    x += 48 + ctx.measureText(item.label).width + 50;
  }
}

function drawLinePanel() {
  const x = margins.left;
  const y = margins.top;
  const base = y + panelHeight;
  const minLogX = Math.log2(0.0625);
  const maxLogX = Math.log2(64);
  const minLogY = Math.log10(0.5);
  const maxLogY = Math.log10(80);
  const sx = (value) => x + ((Math.log2(value) - minLogX) / (maxLogX - minLogX)) * panelWidth;
  const sy = (value) => base - ((Math.log10(value) - minLogY) / (maxLogY - minLogY)) * panelHeight;
  const xTicks = [0.0625, 0.25, 1, 4, 16, 64];
  const xLabels = ["0.06", "0.25", "1", "4", "16", "64"];
  const yTicks = [0.5, 1, 2, 5, 10, 20, 50];

  ctx.fillStyle = ink;
  ctx.font = "600 35px Arial";
  ctx.textAlign = "left";
  ctx.fillText("A. Median daily product dose required", x, y - 61);

  for (const tick of yTicks) {
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
  }

  ctx.strokeStyle = alert;
  ctx.lineWidth = 2;
  ctx.save();
  ctx.setLineDash([8, 6]);
  ctx.beginPath(); ctx.moveTo(x, sy(10)); ctx.lineTo(x + panelWidth, sy(10)); ctx.stroke();
  ctx.restore();
  ctx.fillStyle = alert;
  ctx.font = "24px Arial";
  ctx.textAlign = "left";
  ctx.fillText("10-g/day licensed maximum", x + 12, sy(10) - 19);

  for (const item of classes) {
    const values = grouped.get(item.key);
    ctx.strokeStyle = item.color;
    ctx.lineWidth = 4;
    ctx.beginPath();
    values.forEach((point, index) => {
      const px = sx(point.mic); const py = sy(point.dose);
      if (index === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    });
    ctx.stroke();
    ctx.fillStyle = item.color;
    for (const point of values) {
      ctx.beginPath(); ctx.arc(sx(point.mic), sy(point.dose), 5.5, 0, Math.PI * 2); ctx.fill();
    }
  }

  ctx.strokeStyle = rule;
  ctx.lineWidth = 1.8;
  ctx.beginPath(); ctx.moveTo(x, base); ctx.lineTo(x + panelWidth, base); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x, base); ctx.stroke();
  ctx.fillStyle = ink;
  ctx.font = "25px Arial";
  ctx.textAlign = "center";
  xTicks.forEach((tick, i) => ctx.fillText(xLabels[i], sx(tick), base + 36));
  ctx.font = "28px Arial";
  ctx.fillText("MIC (mg/L)", x + panelWidth / 2, base + 87);
  ctx.save();
  ctx.translate(x - 94, y + panelHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText("Median daily product dose (g/day)", 0, 0);
  ctx.restore();
}

function drawBarPanel() {
  const x = margins.left + panelWidth + gap;
  const y = margins.top;
  const base = y + panelHeight;
  const mics = [4, 8, 16];
  const groupStep = panelWidth / mics.length;
  const barWidth = 42;
  const barGap = 9;
  const blockWidth = classes.length * barWidth + (classes.length - 1) * barGap;
  const sy = (value) => base - (value / 100) * panelHeight;

  ctx.fillStyle = ink;
  ctx.font = "600 35px Arial";
  ctx.textAlign = "left";
  ctx.fillText("B. Subjects within the licensed dose range", x, y - 61);

  for (let tick = 0; tick <= 100; tick += 20) {
    const yy = sy(tick);
    ctx.strokeStyle = tick === 0 ? rule : grid;
    ctx.lineWidth = tick === 0 ? 1.8 : 1;
    ctx.beginPath(); ctx.moveTo(x, yy); ctx.lineTo(x + panelWidth, yy); ctx.stroke();
    ctx.fillStyle = ink;
    ctx.font = "25px Arial";
    ctx.textAlign = "right";
    ctx.fillText(String(tick), x - 14, yy);
  }
  for (let i = 0; i < mics.length; i += 1) {
    const center = x + groupStep * (i + 0.5);
    const start = center - blockWidth / 2;
    for (let j = 0; j < classes.length; j += 1) {
      const item = classes[j];
      const point = grouped.get(item.key).find((row) => row.mic === mics[i]);
      const bx = start + j * (barWidth + barGap);
      const by = sy(point.within);
      ctx.fillStyle = item.color;
      ctx.fillRect(bx, by, barWidth, base - by);
    }
    ctx.fillStyle = ink;
    ctx.font = "600 28px Arial";
    ctx.textAlign = "center";
    ctx.fillText(`MIC ${mics[i]}`, center, base + 43);
  }
  ctx.save();
  ctx.translate(x - 91, y + panelHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = ink;
  ctx.font = "28px Arial";
  ctx.textAlign = "center";
  ctx.fillText("Subjects within licensed range (%)", 0, 0);
  ctx.restore();
}

drawLegend();
drawLinePanel();
drawBarPanel();

await sharp(await canvas.encode("png"))
  .withMetadata({ density: 300 })
  .png()
  .toFile(outPath);
console.log(outPath);
