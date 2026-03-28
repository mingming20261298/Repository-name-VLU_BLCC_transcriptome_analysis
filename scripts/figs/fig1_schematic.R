if (!requireNamespace("DiagrammeR", quietly=TRUE)) install.packages("DiagrammeR", repos="https://cloud.r-project.org")
if (!requireNamespace("DiagrammeRsvg", quietly=TRUE)) install.packages("DiagrammeRsvg", repos="https://cloud.r-project.org")
if (!requireNamespace("rsvg", quietly=TRUE)) install.packages("rsvg", repos="https://cloud.r-project.org")
library(DiagrammeR)
dir.create("~/JTM_GSE84571/figs", showWarnings=FALSE, recursive=TRUE)
gr <- grViz("
digraph G {
  graph [layout = dot, rankdir = LR, fontsize=12]
  node  [shape = box, style=rounded, fontsize=11]
  edge  [fontsize=10]
  A [label='VLU patients\\nRandomized controlled trial']
  B [label='Control\\nStandard compression']
  C [label='BLCC (Apligraf) +\\nStandard compression']
  D [label='Paired biopsies\\nUlcer edge: week0 & week1\\n(per patient)']
  E [label='Microarray (GPL570)\\nCEL -> RMA + QC']
  F [label='Primary model\\nDiD: arm × time\\n(ΔBLCC − ΔControl)']
  G [label='Downstream\\nDEG / Volcano / Heatmap']
  H [label='Pathways\\nGSVA / GSEA\\n(PROGENy/TF optional)']
  I [label='Microenvironment\\nxCell / MCP-counter']
  J [label='Translational output\\nSignature score\\nLOOCV + external validation']
  A -> B
  A -> C
  B -> D
  C -> D
  D -> E
  E -> F
  F -> G
  F -> H
  F -> I
  F -> J
}
")
svg_file <- "~/JTM_GSE84571/figs/Fig1_study_schematic.svg"
png_file <- "~/JTM_GSE84571/figs/Fig1_study_schematic.png"
pdf_file <- "~/JTM_GSE84571/figs/Fig1_study_schematic.pdf"
svg <- DiagrammeRsvg::export_svg(gr)
writeLines(svg, svg_file)
rsvg::rsvg_png(svg_file, file = png_file, width = 2200)
rsvg::rsvg_pdf(svg_file, file = pdf_file)
message("Saved Fig1 to ~/JTM_GSE84571/figs/")
