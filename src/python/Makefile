#
# An example Makefile that builds using Pweave and latexmk.
# The tex file are intermediate/temporary files in this setting,
# although some of the commented code could save tex results
# in a separate tex directory (see the R example).
#
# There is also a commented line that compiles Python noweb using
# R's knitr and knitron packages.
#
#

SRC_DIR = $(realpath $(dir $<))
DEST_DIR = $(abspath ../../output)
FIG_DIR = $(abspath ../../figures)
TEX_DIR = $(abspath ../tex)

.PHONY: all clean clean-tex clean-pweave clean-html clean-md
.PRECIOUS: $(TEX_DIR)/%.tex

PNW_FILES = $(wildcard *.pnw) $(wildcard *.Pnw)
TEXW_FILES = $(wildcard *.texw)
TEX_FILES = $(addsuffix .tex, $(basename $(TEXW_FILES))) 
PDF_FILES = $(addsuffix .pdf, $(basename $(TEXW_FILES))) 
HTML_FILES = $(addsuffix .html, $(basename $(PNW_FILES)))
 
PANDOC_EXTS = +tex_math_dollars+yaml_metadata_block+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash+latex_macros+footnotes+inline_notes
#+citations
PANDOC_OPTIONS = --smart --toc --self-contained -s --mathjax 
#--filter pandoc-citeproc

HTML_DEPS = $(wildcard _output.yaml) 
ifneq ($(wildcard include/before_body.html),) 
PANDOC_OPTIONS += -H include/before_body.html 
HTML_DEPS += include/before_body.html
endif

ifneq ($(wildcard include/after_body.html),) 
PANDOC_OPTIONS += -A include/after_body.html 
HTML_DEPS += include/after_body.html
endif

#PDF_DEPS = ""
PDFLATEX = "pdflatex -interaction nonstopmode -shell-escape "

# add citations with:
# --biblio file.bib --csl chicago-author-date.csl 

all: $(HTML_FILES) $(PDF_FILES)

%.md: %.Pnw
	Pweave -f pandoc -c $<

%.md: $(TEX_DIR)/%.tex
	pandoc -f latex -t markdown$(PANDOC_EXTS) $(PANDOC_OPTIONS) -o $@ $<

%.ipynb: %.md
	notedown $< > $@

%.html: %.md $(HTML_DEPS)
	pandoc -f markdown$(PANDOC_EXTS) -t html $(PANDOC_OPTIONS) -o $@ $<

#
# Pweave results in the tex directory
#
$(TEX_DIR)/%.tex: %.texw
	PynowebWeave -c -d -F $(FIG_DIR) -o $(TEX_DIR)/$(addsuffix .tex, $(basename $<)) $<


%.pdf: $(TEX_DIR)/%.tex
	# in the case that the file was knit'ed outside of this Makefile
	# and not moved to the tex directory, move it now.
	#-mv -f $(addsuffix .tex,k $(basename $<)) $(TEX_DIR)
	make -C $(TEX_DIR) -f Makefile $@

clean: clean-pdf clean-texw clean-html clean-md

clean-html:
	-$(RM) *.html

clean-md:
	-$(RM) *.html

clean-pdf:
	-cd $(DEST_DIR) && latexmk -f -C -bibtex ${PDF_FILES}
	-latexmk -f -C -bibtex ${PDF_FILES} 
	#-make -C $(TEX_DIR) -f Makefile clean

clean-texw:
	-cd $(TEX_DIR) && $(RM) ${TEX_FILES}
	-$(RM) -rf _minted-
	-$(RM) -rf cache
	-$(RM) -rf figures
 

