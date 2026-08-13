import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createCanvas } from "@napi-rs/canvas";
import sharp from "sharp";

const here = path.dirname(fileURLToPath(import.meta.url));
const csvPath = path.join(here, "outputs", "lung_penetration_icu_trial_detail.csv");
const outPath = path.join(here, "figures", "fig_lung_penetration_icu_updated.png");

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const split = (line) => {
    const out = []; let field = ""; let quoted = false;
    for (let i = 0; i < line.length; i += 1) {
      const char = line[i];
      if (char === '"' && line[i + 1] === '"') { field += '"'; i += 1; }
      else if (char === '"') quoted = !quoted;
      else if (char === "," && !quoted) { out.push(field); field = ""; }
      else field += char;
    }
    out.push(field); return out;
  };
  const header = split(lines[0]);
  return lines.slice(1).map((line) => Object.fromEntries(split(line).map((value, i) => [header[i], value])));
}

const rows = parseCsv(await fs.readFile(csvPath, "utf8"));
const regimens = [
  { id: "R1", ekfc: "0–30" }, { id: "R8", ekfc: "31–60" },
  { id: "R10", ekfc: "61–90" }, { id: "R12", ekfc: "91–120" },
  { id: "R13", ekfc: "121–150" },
];
const scenarios = [
  { key: "plasma", label: "Plasma", color: "#1F5B99" },
  { key: "icu_trial", label: "ELF, ICU trial (0.41 / 0.44)", color: "#2F855A" },
  { key: "healthy_volunteer", label: "ELF, healthy volunteers (0.52 / 0.42)", color: "#C87320" },
  { key: "conservative", label: "ELF, conservative (0.30 / 0.30)", color: "#8E3B2D" },
];
const cell = new Map(rows.map((row) => [`${row.scenario}|${row.regimen}`, row]));

const width = 2400, height = 1220;
const canvas = createCanvas(width, height);
const ctx = canvas.getContext("2d");
ctx.fillStyle = "#FFFFFF"; ctx.fillRect(0, 0, width, height); ctx.textBaseline = "middle";
const ink = "#18212B", grid = "#D9E0E7", rule = "#6C7783";
const margins = { left: 145, right: 55, top: 218, bottom: 140 }, gap = 120;
const panelWidth = (width - margins.left - margins.right - gap) / 2;
const panelHeight = height - margins.top - margins.bottom;

function drawLegend() {
  ctx.font = "500 29px Arial"; ctx.textAlign = "left";
  let x = margins.left, y = 65;
  for (const item of scenarios) {
    ctx.fillStyle = item.color; ctx.fillRect(x, y - 10, 32, 20);
    ctx.fillStyle = ink; ctx.fillText(item.label, x + 45, y);
    x += 45 + ctx.measureText(item.label).width + 56;
  }
  ctx.fillStyle = "#4C5966"; ctx.font = "25px Arial";
  ctx.fillText("ELF/plasma ratios: ceftazidime / avibactam", margins.left, 122);
}

function drawPanel(index, heading, key, axisLabel) {
  const x = margins.left + index * (panelWidth + gap), y = margins.top, base = y + panelHeight;
  const groupStep = panelWidth / regimens.length, barWidth = 31, barGap = 8;
  const groupWidth = scenarios.length * barWidth + (scenarios.length - 1) * barGap;
  const sy = (value) => base - (value / 100) * panelHeight;
  ctx.fillStyle = ink; ctx.font = "600 35px Arial"; ctx.textAlign = "left";
  ctx.fillText(`${index === 0 ? "A" : "B"}. ${heading}`, x, y - 60);
  for (let tick = 0; tick <= 100; tick += 20) {
    const yy = sy(tick); ctx.strokeStyle = tick === 0 ? rule : grid; ctx.lineWidth = tick === 0 ? 1.8 : 1;
    ctx.beginPath(); ctx.moveTo(x, yy); ctx.lineTo(x + panelWidth, yy); ctx.stroke();
    ctx.fillStyle = ink; ctx.font = "25px Arial"; ctx.textAlign = "right"; ctx.fillText(String(tick), x - 15, yy);
  }
  const ninety = sy(90); ctx.save(); ctx.strokeStyle = "#7E8791"; ctx.lineWidth = 1.4; ctx.setLineDash([5, 6]);
  ctx.beginPath(); ctx.moveTo(x, ninety); ctx.lineTo(x + panelWidth, ninety); ctx.stroke(); ctx.restore();
  ctx.fillStyle = "#56616D"; ctx.font = "22px Arial"; ctx.textAlign = "right"; ctx.fillText("90%", x + panelWidth - 1, ninety - 18);
  for (let g = 0; g < regimens.length; g += 1) {
    const group = regimens[g], centre = x + groupStep * (g + 0.5), start = centre - groupWidth / 2;
    for (let s = 0; s < scenarios.length; s += 1) {
      const scenario = scenarios[s], row = cell.get(`${scenario.key}|${group.id}`), value = Number(row[key]);
      const bx = start + s * (barWidth + barGap), by = sy(value);
      ctx.fillStyle = scenario.color; ctx.fillRect(bx, by, barWidth, base - by);
    }
    ctx.fillStyle = ink; ctx.font = "600 26px Arial"; ctx.textAlign = "center"; ctx.fillText(group.id, centre, base + 34);
    ctx.font = "23px Arial"; ctx.fillText(`EKFC ${group.ekfc}`, centre, base + 69);
  }
  ctx.save(); ctx.translate(x - 95, y + panelHeight / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = ink; ctx.font = "28px Arial"; ctx.textAlign = "center"; ctx.fillText(axisLabel, 0, 0); ctx.restore();
}

drawLegend();
drawPanel(0, "MIC-weighted joint CFR", "joint_cfr_pct", "Joint CFR (%)");
drawPanel(1, "Joint PTA at MIC 8 mg/L", "joint_pta_mic8_pct", "Joint PTA (%)");
await sharp(await canvas.encode("png")).withMetadata({ density: 300 }).png().toFile(outPath);
console.log(outPath);
