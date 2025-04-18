# Load the required libraries
library(pacman)
p_load(DiagrammeR)
p_load(DiagrammeRsvg) # For exporting SVG
p_load(rsvg)          # For converting SVG

# --- Diagram 1: MAPseq Single Run Pipeline ---
mapseq_run_diagram <- grViz("
digraph mapseq_run {

  # Graph aesthetic settings
  graph [layout = dot, rankdir = TB, overlap = false, splines=ortho, bgcolor = transparent]
  node [shape = box, style = filled, fillcolor = LightBlue, fontname = Helvetica, penwidth=2]
  edge [color = DimGray, arrowhead = vee, penwidth=1.5]

  # Define Nodes
  start_run [label = 'Start Single Run\\n(run_mapseq.sh)', fillcolor=PaleGreen]
  config [label = 'Source project_config.txt']
  preflight [label = 'Perform Pre-flight Checks\\n(preflight_checks.sh)', fillcolor=Khaki]
  genfastq [label = 'Demultiplex FASTQs (Optional)\\n(run_genfastq.sh)']

  rna_check [label = 'RNA FASTQs available?', shape=diamond, fillcolor=LightGoldenrodYellow]
  atac_check [label = 'ATAC FASTQs available?', shape=diamond, fillcolor=LightGoldenrodYellow]

  run_rna [label = 'Run Cell Ranger RNA Multi\\n(run_cellranger_RNA.FB.VDJ.sh)']
  run_asap2kite [label = 'Convert ASAP to Kite FASTQs\\n(run_asap_to_kite.sh)']
  run_kite [label = 'Run Kite for ASAP\\n(run_kite.sh)']
  run_atac [label = 'Run Cell Ranger ATAC Count\\n(run_cellranger_ATAC.sh)']

  copy_rna_reports [label = 'Copy RNA Reports']
  copy_atac_reports [label = 'Copy ATAC Reports']

  end_run [label = 'End Single Run', fillcolor=LightPink]


  # Define Edges (Workflow)
  start_run -> config
  config -> preflight
  preflight -> genfastq
  genfastq -> rna_check
  genfastq -> atac_check

  rna_check -> run_rna [label = 'Yes']
  run_rna -> copy_rna_reports
  copy_rna_reports -> end_run [style=dashed]
  rna_check -> end_run [label = 'No', style=dashed]


  atac_check -> run_asap2kite [label = 'Yes']
  run_asap2kite -> run_kite
  run_kite -> run_atac
  run_atac -> copy_atac_reports
  copy_atac_reports -> end_run [style=dashed]
  atac_check -> end_run [label = 'No', style=dashed]

}
")

# --- Diagram 2: MAPseq Aggregation Pipeline ---
mapseq_aggr_diagram <- grViz("
digraph mapseq_aggr {

  # Graph aesthetic settings
  graph [layout = dot, rankdir = TB, overlap = false, splines=ortho, bgcolor = transparent]
  node [shape = box, style = filled, fillcolor = LightSkyBlue, fontname = Helvetica, penwidth=2]
  edge [color = DarkSlateGray, arrowhead = vee, penwidth=1.5]

  # Define Nodes
  start_aggr [label = 'Start Aggregation\\n(run_mapseq_aggr.sh)', fillcolor=PaleGreen]
  aggr_config [label = 'Source project_config.txt']
  aggr_rna [label = 'Run Cell Ranger RNA Aggr (background)\\n(run_cellranger.aggr_rna.sh)']
  aggr_atac [label = 'Run Cell Ranger ATAC Aggr (background)\\n(run_cellranger.aggr_atac.sh)']
  wait [label = 'Wait for Aggr jobs', shape=ellipse, fillcolor=Khaki]
  demux_rna [label = 'Demultiplex RNA with Seurat\\n(run_seurat.demux_rna.R)']
  demux_atac [label = 'Demultiplex ATAC/ASAP with Seurat\\n(run_seurat.demux_atac.R)']
  end_aggr [label = 'End Aggregation\\n(RDS Objects Created)', fillcolor=LightPink]

  # Define Edges (Workflow)
  start_aggr -> aggr_config
  aggr_config -> aggr_rna
  aggr_config -> aggr_atac
  aggr_rna -> wait [style=dashed]
  aggr_atac -> wait [style=dashed]
  wait -> demux_rna
  wait -> demux_atac
  demux_rna -> end_aggr [style=dashed]
  demux_atac -> end_aggr [style=dashed]

}
")


# Print the flowchart objects
print(mapseq_run_diagram)
print(mapseq_aggr_diagram)

# --- Saving the Diagrams ---

# Define the output file names
svg_file_run <- "mapseq_run_pipeline.svg"
png_file_run <- "mapseq_run_pipeline.png"
pdf_file_run <- "mapseq_run_pipeline.pdf"

svg_file_aggr <- "mapseq_aggr_pipeline.svg"
png_file_aggr <- "mapseq_aggr_pipeline.png"
pdf_file_aggr <- "mapseq_aggr_pipeline.pdf"

# Function to save a diagram
save_diagram <- function(diagram_obj, svg_file, png_file, pdf_file) {
  # 1. Save as SVG
  svg_output <- tryCatch({
     export_svg(diagram_obj)
  }, error = function(e) {
    message("Could not export SVG for: ", svg_file)
    print(e)
    return(NULL)
  })

  if (!is.null(svg_output)) {
    tryCatch({
        conn <- file(svg_file, "w")
        writeLines(svg_output, conn)
        close(conn)
        cat("Diagram saved as SVG:", svg_file, "\n")

        # Check if rsvg is available before attempting conversion
        if (requireNamespace("rsvg", quietly = TRUE)) {
          # 2. Convert SVG to PNG
          tryCatch({
            rsvg_png(svg = svg_file, file = png_file)
            cat("Diagram saved as PNG:", png_file, "\n")
          }, error = function(e_png) {
            message("Could not convert SVG to PNG: ", png_file, ". Ensure librsvg is installed on your system.")
            print(e_png)
          })

          # 3. Convert SVG to PDF
          tryCatch({
            rsvg_pdf(svg = svg_file, file = pdf_file)
            cat("Diagram saved as PDF:", pdf_file, "\n")
          }, error = function(e_pdf) {
            message("Could not convert SVG to PDF: ", pdf_file, ". Ensure librsvg is installed on your system.")
            print(e_pdf)
          })
        } else {
          message("Install the 'rsvg' package and ensure librsvg is installed on your system to convert SVG to PNG/PDF.")
        }
    }, error = function(e_write) {
        message("Could not write SVG file: ", svg_file)
        print(e_write)
    }) # End tryCatch for writing/converting
  } # End if !is.null(svg_output)
} # End save_diagram function

# Save the diagrams
save_diagram(mapseq_run_diagram, svg_file_run, png_file_run, pdf_file_run)
save_diagram(mapseq_aggr_diagram, svg_file_aggr, png_file_aggr, pdf_file_aggr)

# Optional: Clean up the intermediate SVG files if desired
# file.remove(svg_file_run)
# file.remove(svg_file_aggr)
