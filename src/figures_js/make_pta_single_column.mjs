import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createCanvas } from "@napi-rs/canvas";
import sharp from "sharp";

const here = path.dirname(fileURLToPath(import.meta.url));
const csvPath = path.join(here, "outputs", "primary_pta_results.csv");
const outPath = path.join(here, "figures", "fig_pta_vs_mic_single_column.png");

function parseCsv(text) {
  const [header, ...lines] = text.trim().split(/\r?\n/);
  const keys = header.split(",");
  return lines.map((line) => Object.fromEntries(line.split(",").map((value, i) => [keys[i], value])));
}

const rows = parseCsv(await fs.readFile(csvPath, "utf8"));
const selected = [
  ["R1", "0–30"], ["R8", "31–60"], ["R10", "61–90"],
  ["R12", "91–120"], ["R13", "121–150"],
];
const byRegimen = new Map(selected.map(([regimen]) => [regimen, []]));
for (const row of rows) {
  if (byRegimen.has(row.regimen)) byRegimen.get(row.regimen).push({
    mic: Number(row.mic_mg_l),
    caz: Number(row.caz_pta_pct),
    avi: Number(row.avi_attainment_pct),
    joint: Number(row.joint_pta_pct),
  });
}
for (const values of byRegimen.values()) values.sort((a, b) => a.mic - b.mic);

const width = 1080;
const height = 1600;
const canvas = createCanvas(width, height);
const ctx = canvas.getContext("2d");
ctx.fillStyle = "#FFFFFF";
ctx.fillRect(0, 0, width, height);
ctx.textBaseline = "middle";

const palette = { caz: "#2673C9", avi: "#C84D45", joint: "#1A1A1A" };
const margins = { left: 125, right: 36, top: 70, bottom: 100 };
const gapX = 54;
const gapY = 82;
const plotWidth = (width - margins.left - margins.right - gapX) / 2;
const plotHeight = (height - margins.top - margins.bottom - 2 * gapY) / 3;
const minLog = Math.log2(0.0625);
const maxLog = Math.log2(64);
const xTicks = [0.0625, 0.25, 1, 4, 16, 64];
const xTickLabels = ["0.06", "0.25", "1", "4", "16", "64"];

function panelPosition(index) {
  const col = index % 2;
  const row = Math.floor(index / 2);
  return {
    x: margins.left + col * (plotWidth + gapX),
    y: margins.top + row * (plotHeight + gapY),
  };
}

function sx(value, x) {
  return x + ((Math.log2(value) - minLog) / (maxLog - minLog)) * plotWidth;
}

function sy(value, y) {
  return y + plotHeight - (value / 100) * plotHeight;
}

function line(points, color, marker) {
  ctx.strokeStyle = color;
  ctx.lineWidth = marker === "triangle" ? 5 : 4;
  ctx.beginPath();
  points.forEach((point, index) => {
    if (index === 0) ctx.moveTo(point.x, point.y);
    else ctx.lineTo(point.x, point.y);
  });
  ctx.stroke();
  ctx.fillStyle = color;
  for (const point of points) {
    if (marker === "circle") {
      ctx.beginPath(); ctx.arc(point.x, point.y, 5, 0, Math.PI * 2); ctx.fill();
    } else if (marker === "square") {
      ctx.fillRect(point.x - 4.5, point.y - 4.5, 9, 9);
    } else {
      ctx.beginPath();
      ctx.moveTo(point.x, point.y - 6);
      ctx.lineTo(point.x - 5.5, point.y + 5);
      ctx.lineTo(point.x + 5.5, point.y + 5);
      ctx.closePath(); ctx.fill();
    }
  }
}

function drawPanel(index, regimen, ekfc) {
  const { x, y } = panelPosition(index);
  const values = byRegimen.get(regimen);
  ctx.strokeStyle = "#222222";
  ctx.lineWidth = 1.2;
  ctx.strokeRect(x, y, plotWidth, plotHeight);

  ctx.strokeStyle = "#E1E1E1";
  ctx.lineWidth = 1;
  for (const value of [0, 50, 100]) {
    ctx.beginPath(); ctx.moveTo(x, sy(value, y)); ctx.lineTo(x + plotWidth, sy(value, y)); ctx.stroke();
  }
  for (const tick of xTicks) {
    ctx.beginPath(); ctx.moveTo(sx(tick, x), y); ctx.lineTo(sx(tick, x), y + plotHeight); ctx.stroke();
  }

  ctx.setLineDash([8, 7]);
  ctx.strokeStyle = "#777777";
  ctx.lineWidth = 1.6;
  ctx.beginPath(); ctx.moveTo(sx(8, x), y); ctx.lineTo(sx(8, x), y + plotHeight); ctx.stroke();
  ctx.setLineDash([3, 6]);
  ctx.strokeStyle = "#8B8B8B";
  ctx.beginPath(); ctx.moveTo(x, sy(90, y)); ctx.lineTo(x + plotWidth, sy(90, y)); ctx.stroke();
  ctx.setLineDash([]);

  line(values.map((v) => ({ x: sx(v.mic, x), y: sy(v.caz, y) })), palette.caz, "circle");
  line(values.map((v) => ({ x: sx(v.mic, x), y: sy(v.avi, y) })), palette.avi, "square");
  line(values.map((v) => ({ x: sx(v.mic, x), y: sy(v.joint, y) })), palette.joint, "triangle");

  ctx.fillStyle = "#111111";
  ctx.font = "600 28px Arial";
  ctx.textAlign = "center";
  ctx.fillText(`EKFC ${ekfc}  •  ${regimen}`, x + plotWidth / 2, y - 22);
  ctx.font = "24px Arial";
  ctx.textAlign = "center";
  xTicks.forEach((tick, i) => ctx.fillText(xTickLabels[i], sx(tick, x), y + plotHeight + 30));
  if (index % 2 === 0) {
    ctx.textAlign = "right";
    [0, 50, 100].forEach((value) => ctx.fillText(String(value), x - 12, sy(value, y)));
  }
  if (index % 2 === 0) {
    ctx.save();
    ctx.translate(x - 84, y + plotHeight / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.font = "25px Arial";
    ctx.fillText("Target attainment (%)", 0, 0);
    ctx.restore();
  }
  if (index >= 3) {
    ctx.textAlign = "center";
    ctx.font = "25px Arial";
    ctx.fillText("MIC (mg/L)", x + plotWidth / 2, y + plotHeight + 66);
  }
}

selected.forEach(([regimen, ekfc], index) => drawPanel(index, regimen, ekfc));

const legend = panelPosition(5);
ctx.fillStyle = "#111111";
ctx.font = "600 29px Arial";
ctx.textAlign = "left";
ctx.fillText("Key", legend.x, legend.y + 24);
const legendItems = [
  [palette.caz, "circle", "Ceftazidime"],
  [palette.avi, "square", "Avibactam"],
  [palette.joint, "triangle", "Joint"],
];
legendItems.forEach(([color, marker, label], index) => {
  const y = legend.y + 74 + index * 52;
  line([{ x: legend.x, y }, { x: legend.x + 58, y }], color, marker);
  ctx.fillStyle = "#111111";
  ctx.font = "25px Arial";
  ctx.fillText(label, legend.x + 76, y);
});
ctx.setLineDash([8, 7]);
ctx.strokeStyle = "#777777";
ctx.lineWidth = 1.6;
ctx.beginPath(); ctx.moveTo(legend.x, legend.y + 252); ctx.lineTo(legend.x + 58, legend.y + 252); ctx.stroke();
ctx.setLineDash([3, 6]);
ctx.strokeStyle = "#8B8B8B";
ctx.beginPath(); ctx.moveTo(legend.x, legend.y + 302); ctx.lineTo(legend.x + 58, legend.y + 302); ctx.stroke();
ctx.setLineDash([]);
ctx.fillStyle = "#111111";
ctx.font = "24px Arial";
ctx.fillText("EUCAST breakpoint", legend.x + 76, legend.y + 252);
ctx.fillText("90% PTA", legend.x + 76, legend.y + 302);

await sharp(await canvas.encode("png")).withMetadata({ density: 300 }).png().toFile(outPath);
console.log(outPath);
