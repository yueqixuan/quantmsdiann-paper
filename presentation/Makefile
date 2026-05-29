# Build the quantmsdiann results presentation.
# Figures are pulled from ../paper/figures/pdf/ — build them first with:
#   cd ../paper && make figures

SLIDES    := slides
LATEXMK   := latexmk
BUILD_DIR := build

.PHONY: all clean

all: $(SLIDES).pdf

$(SLIDES).pdf: $(SLIDES).tex
	$(LATEXMK) -pdf -interaction=nonstopmode -outdir=$(BUILD_DIR) $(SLIDES).tex
	cp $(BUILD_DIR)/$(SLIDES).pdf .

clean:
	$(LATEXMK) -C -outdir=$(BUILD_DIR) $(SLIDES).tex
	rm -rf $(BUILD_DIR)
	rm -f $(SLIDES).pdf
