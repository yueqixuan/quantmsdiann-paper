# Build the quantmsdiann MCP manuscript PDF and convert analysis SVGs to PDF.
#
# Requires: pdflatex, bibtex, and one of inkscape / rsvg-convert / cairosvg
# for figure conversion.

MAIN      := main
SUPP      := supplementary
GA        := graphical_abstract
TEXDIR    := .
FIG_SRC   := ../analysis/figures
FIG_PDF   := figures/pdf
LATEXMK   := latexmk
BUILD_DIR := build

SVGS := \
	manuscript/fig1_architecture_scaling.svg \
	performance/parallelism_vs_wallclock.svg \
	performance/runtime_per_step.svg \
	performance/resources_per_step.svg \
	performance/queue_size_sweep.svg \
	quantmsdiann_benchmarks/main_benchmarks_overview.svg \
	quantmsdiann_benchmarks/main_benchmarks_precursors.svg \
	quantmsdiann_benchmarks/main_diann_quantmsdiann_parity.svg \
	quantmsdiann_benchmarks/main_accuracy.svg \
	quantmsdiann_benchmarks/supp_vs_proteobench_matched_min3.svg \
	quantmsdiann_benchmarks/supp_vs_proteobench_min3.svg \
	PXD003539/main_comparison.svg \
	PXD003539/supp_genes_per_condition.svg \
	PXD003539/supp_missing_values_per_run.svg \
	PXD003539/supp_peptides_per_protein.svg \
	PXD003539/supp_venn_protein_accessions.svg \
	PXD004701/main_comparison.svg \
	PXD004701/supp_missing_values_per_run.svg \
	PXD004701/supp_proteins_per_subtype.svg \
	PXD030304/main_comparison.svg \
	PXD030304/supp_missing_values_per_run.svg \
	PXD030304/supp_proteins_per_tissue.svg \
	plexDIA/MSV000093870/main_galatidou_total_pg.svg \
	plexDIA/MSV000093870/main_galatidou_comparison.svg \
	PXD064049/main_comparison.svg \
	combined/atlas_main.svg

PDF_TARGETS := $(addprefix $(FIG_PDF)/,$(SVGS:.svg=.pdf)) \
	$(FIG_PDF)/manuscript/workflow_panel.pdf

WORKFLOW_SVG    := figures/source/quantmsdiann_workflow.svg
WORKFLOW_DPI    := 600

.PHONY: all pdf supplementary graphical-abstract figures clean distclean

all: figures pdf supplementary

pdf: $(BUILD_DIR)/$(MAIN).pdf $(BUILD_DIR)/$(GA).pdf
	cp $(BUILD_DIR)/$(MAIN).pdf .
	cp $(BUILD_DIR)/$(GA).pdf .

supplementary: $(BUILD_DIR)/$(SUPP).pdf
	cp $(BUILD_DIR)/$(SUPP).pdf .

graphical-abstract: $(BUILD_DIR)/$(GA).pdf
	cp $(BUILD_DIR)/$(GA).pdf .

$(BUILD_DIR)/$(MAIN).pdf: $(MAIN).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -outdir=$(BUILD_DIR) $(MAIN).tex

$(BUILD_DIR)/$(SUPP).pdf: $(SUPP).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -outdir=$(BUILD_DIR) $(SUPP).tex

$(BUILD_DIR)/$(GA).pdf: $(GA).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -outdir=$(BUILD_DIR) $(GA).tex

figures: $(PDF_TARGETS)

$(FIG_PDF)/%.pdf: $(FIG_SRC)/%.svg | $(FIG_PDF)
	@mkdir -p $(@D)
	@if command -v rsvg-convert >/dev/null 2>&1; then \
		rsvg-convert -f pdf -o "$@" "$<"; \
	elif command -v inkscape >/dev/null 2>&1; then \
		inkscape "$<" --export-type=pdf --export-filename="$@"; \
	elif python3 -c "import cairosvg" >/dev/null 2>&1; then \
		python3 -c "import cairosvg; cairosvg.svg2pdf(url='$<', write_to='$@')"; \
	elif command -v magick >/dev/null 2>&1; then \
		magick -background white -density 300 "$<" "$@"; \
	else \
		echo "No SVG converter found (rsvg-convert, inkscape, cairosvg, or magick)." >&2; \
		echo "Install one, or open SVGs directly while drafting." >&2; \
		exit 1; \
	fi

$(FIG_PDF)/manuscript/workflow_panel.pdf: $(WORKFLOW_SVG) | $(FIG_PDF)
	mkdir -p $(@D)
	@if command -v rsvg-convert >/dev/null 2>&1; then \
		rsvg-convert -f pdf -o "$@" "$<"; \
	elif command -v inkscape >/dev/null 2>&1; then \
		inkscape "$<" --export-type=pdf --export-filename="$@"; \
	elif command -v magick >/dev/null 2>&1; then \
		magick -background white -density $(WORKFLOW_DPI) "$<" "$@"; \
	else \
		echo "No SVG converter found for workflow panel." >&2; \
		exit 1; \
	fi

$(FIG_PDF):
	mkdir -p $(FIG_PDF)

clean:
	$(LATEXMK) -C -outdir=$(BUILD_DIR) $(MAIN).tex
	$(LATEXMK) -C -outdir=$(BUILD_DIR) $(SUPP).tex
	$(LATEXMK) -C -outdir=$(BUILD_DIR) $(GA).tex
	rm -rf $(BUILD_DIR)
	rm -f $(MAIN).pdf $(SUPP).pdf $(GA).pdf

distclean: clean
	rm -rf $(FIG_PDF)
