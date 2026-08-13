import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const workbookPath = process.argv[2];
const mode = process.argv[3] ?? "inspect";
const previewPath = process.argv[4] ?? "supplementary-preview.png";

const input = await FileBlob.load(workbookPath);
const workbook = await SpreadsheetFile.importXlsx(input);

if (mode === "inspect") {
  const summary = await workbook.inspect({
    kind: "workbook,sheet,table,computedStyle",
    sheetId: "S1 Label Taxonomy",
    range: "A1:D12",
    maxChars: 7000,
    tableMaxRows: 8,
    tableMaxCols: 8,
  });
  console.log(summary.ndjson);
  const preview = await workbook.render({
    sheetName: "S1 Label Taxonomy",
    range: "A1:D20",
    autoCrop: "all",
    scale: 1.5,
    format: "png",
  });
  await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));
  process.exit(0);
}

const name = "S2a Evidence Sources";
const existing = workbook.worksheets.getItemOrNullObject(name);
if (!existing.isNullObject) {
  existing.delete();
}
const sheet = workbook.worksheets.add(name);
sheet.showGridLines = false;

sheet.getRange("A1:G1").merge();
sheet.getRange("A1").values = [["Supplementary Table S2a. Evidence-source roles in the primary analysis"]];
sheet.getRange("A2:G2").merge();
sheet.getRange("A2").values = [["Sources are classified by their role in the prespecified simulation framework. Contextual sources did not contribute numerical inputs to the primary model."]];

sheet.getRange("A4:G4").values = [["A. SOURCE ROLES", null, null, null, null, null, null]];
sheet.getRange("A5:G5").values = [[
  "Source", "Population", "Component contributed", "Analytical role",
  "Used in primary model", "Used in sensitivity analyses", "Used for calibration"
]];
sheet.getRange("A6:G13").values = [
  ["Cojutti et al.", "ICU, non-RRT", "CAZ/AVI clearance equations, IIV, correlation", "Primary input", "Yes", "Yes", "Yes"],
  ["Gatti et al.", "CVVHDF", "Exposure and circuit context", "Donor sensitivity input", "No", "Yes", "No"],
  ["Curtiaud et al.", "ECMO", "Trough concentrations", "Contextual PK anchor", "No", "No", "No"],
  ["O'Jeanson et al.", "CVVHDF simulation", "PK and circuit removal", "Donor sensitivity input", "No", "Yes", "No"],
  ["Wölfl-Duchek et al.", "Healthy volunteers/CSF", "Target-site PK", "Contextual evidence", "No", "No", "No"],
  ["Aubry et al.", "In vitro", "Killing and inoculum-effect model", "Contextual evidence", "No", "No", "No"],
  ["Lee et al.", "Surveillance", "KPC/OXA-48-like MIC distributions", "Primary MIC input", "Yes", "Yes", "No"],
  ["Bakthavatchalam et al.", "Surveillance", "Alternative OXA-48-like MIC distribution", "Donor sensitivity input", "No", "Yes", "No"],
];

sheet.getRange("A1:G1").format = {
  font: { bold: true, size: 12, name: "Arial" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
sheet.getRange("A2:G2").format = {
  font: { italic: true, color: "#555555", size: 10, name: "Arial" },
  wrapText: true,
  verticalAlignment: "center",
};
sheet.getRange("A4:G4").format = {
  font: { bold: true, size: 10, name: "Arial" },
  horizontalAlignment: "left",
  verticalAlignment: "center",
};
sheet.getRange("A5:G5").format = {
  fill: "#2F5496",
  font: { bold: true, color: "#FFFFFF", size: 10, name: "Arial" },
  wrapText: true,
  horizontalAlignment: "center",
  verticalAlignment: "center",
  borders: { preset: "all", style: "thin", color: "#000000" },
};
sheet.getRange("A6:G13").format = {
  font: { size: 10, name: "Arial" },
  wrapText: true,
  verticalAlignment: "center",
  borders: { preset: "all", style: "thin", color: "#000000" },
};
sheet.getRange("A7:G7").format.fill = "#D9E6F2";
sheet.getRange("A9:G9").format.fill = "#D9E6F2";
sheet.getRange("A11:G11").format.fill = "#D9E6F2";
sheet.getRange("A13:G13").format.fill = "#D9E6F2";
sheet.getRange("E6:G13").format.horizontalAlignment = "center";
sheet.getRange("A1:G1").format.rowHeight = 24;
sheet.getRange("A2:G2").format.rowHeight = 26;
sheet.getRange("A4:G4").format.rowHeight = 22;
sheet.getRange("A5:G5").format.rowHeight = 34;
sheet.getRange("A6:G13").format.rowHeight = 34;

const widths = [22, 22, 38, 24, 18, 22, 18];
for (let i = 0; i < widths.length; i += 1) {
  sheet.getRangeByIndexes(0, i, 13, 1).format.columnWidth = widths[i];
}
sheet.freezePanes.freezeRows(5);

const preview = await workbook.render({
  sheetName: name,
  range: "A1:G13",
  autoCrop: "all",
  scale: 1.5,
  format: "png",
});
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(workbookPath);
