import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createCanvas } from "@napi-rs/canvas";
import sharp from "sharp";

const here = path.dirname(fileURLToPath(import.meta.url));
const csvPath = path.join(here, "outputs", "lung_penetration.csv");
const outPath = path.join(here, "figures", "fig_lung_penetration_redesigned.png");

function parseCsv(text) {
  const [header, ...lines] = text.trim().split(/\r?\n/);
  const keys = header.split(",");
  return lines.map((line) => Object.fromEntries(line.split(",").map((value, i) => [keys[i], value])));
}

const rows = parseCsv(await fs.readFile(csvPath, "utf8"));
const regimens = [
  { id: "R1", label: "R1", ekfc: "0\u201330" },
  { id: "R8", label: "R8", ekfc: "31\u201360" },
  { id: "R10", label: "R10", ekfc: "61\u201390" },
  { id: "R12", label: "R12", ekfc: "91\u2013120" },
  { id: "R13", label: "R13", ekfc: "121\u2013150" },
];
const scenarios = [
  { key: "plasma", label: "Plasma", color: "#2B6CB0" },
  { key: "central estimate", label: "ELF central (CAZ/AVI 0.52/0.42)", color: "#D97706" },
  { key: "conservative", label: "ELF conservative (CAZ/AVI 0.30/0.30)", color: "#8C3B20" },
];
const byCell = new Map(rows.map((row) => [`${row.compartment}|${row.regimen}`, row]));

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
const margins = { left: 122, right: 52, top: 172, bottom: 128 };
const gap = 92;
const panelWidth = (width - margins.left - margins.right - gap) / 2;
const panelHeight = height - margins.top - margins.bottom;

function drawLegend() {
  ctx.font = "500 29px Arial";
  let x = margins.left;
  const y = 66;
  for (const item of scenarios) {
    ctx.fillStyle = item.color;
    ctx.fillRect(x, y - 10, 30, 20);
    ctx.fillStyle = ink;
    ctx.textAlign = "left";
    ctx.fillText(item.label, x + 43, y);
    x += 43 + ctx.measureText(item.label).width + 58;
  }
}

function drawPanel(index, heading, valueKey, yLabel) {
  const x = margins.left + index * (panelWidth + gap);
  const y = margins.top;
  const baseline = y + panelHeight;
  const groupStep = panelWidth / regimens.length;
  const barWidth = 38;
  const barGap = 10;
  const groupBarWidth = scenarios.length * barWidth + (scenarios.length - 1) * barGap;
  const sy = (value) => baseline - (value / 100) * panelHeight;

  ctx.fillStyle = ink;
  ctx.font = "600 35px Arial";
  ctx.textAlign = "left";
  ctx.fillText(`${index === 0 ? "A" : "B"}. ${heading}`, x, y - 61);

  for (let value = 0; value <= 100; value += 20) {
    const yy = sy(value);
    ctx.strokeStyle = value === 0 ? rule : grid;
    ctx.lineWidth = value === 0 ? 1.8 : 1;
    ctx.beginPath();
    ctx.moveTo(x, yy);
    ctx.lineTo(x + panelWidth, yy);
    ctx.stroke();
    ctx.fillStyle = ink;
    ctx.font = "25px Arial";
    ctx.textAlign = "right";
    ctx.fillText(String(value), x - 14, yy);
  }

  const ninetyY = sy(90);
  ctx.save();
  ctx.strokeStyle = "#7E8791";
  ctx.lineWidth = 1.5;
  ctx.setLineDash([5, 6]);
  ctx.beginPath();
  ctx.moveTo(x, ninetyY);
  ctx.lineTo(x + panelWidth, ninetyY);
  ctx.stroke();
  ctx.restore();
  ctx.fillStyle = "#56616D";
  ctx.font = "22px Arial";
  ctx.textAlign = "right";
  ctx.fillText("90%", x + panelWidth - 2, ninetyY - 18);

  for (let g = 0; g < regimens.length; g += 1) {
    const groupCenter = x + groupStep * (g + 0.5);
    const start = groupCenter - groupBarWidth / 2;
    for (let s = 0; s < scenarios.length; s += 1) {
      const scenario = scenarios[s];
      const row = byCell.get(`${scenario.key}|${regimens[g].id}`);
      const value = Number(row[valueKey]);
      const bx = start + s * (barWidth + barGap);
      const by = sy(value);
      ctx.fillStyle = scenario.color;
      ctx.fillRect(bx, by, barWidth, baseline - by);
    }

    ctx.fillStyle = ink;
    ctx.font = "600 26px Arial";
    ctx.textAlign = "center";
    ctx.fillText(regimens[g].label, groupCenter, baseline + 35);
    ctx.font = "23px Arial";
    ctx.fillText(`EKFC ${regimens[g].ekfc}`, groupCenter, baseline + 69);
  }

  ctx.save();
  ctx.translate(x - 86, y + panelHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = ink;
  ctx.font = "28px Arial";
  ctx.textAlign = "center";
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();
}

drawLegend();
drawPanel(0, "MIC-weighted joint CFR", "joint_cfr_pct", "Joint CFR (%)");
drawPanel(1, "Joint PTA at MIC 8 mg/L", "joint_pta_mic8_pct", "Joint PTA (%)");

await sharp(await canvas.encode("png"))
  .withMetadata({ density: 300 })
  .png()
  .toFile(outPath);
console.log(outPath);
