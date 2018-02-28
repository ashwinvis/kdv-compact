
name := KdV_JSC
dir := Latex
path := $(dir)/$(name)

latex := cd $(dir) && pdflatex $(name).tex

latex_files := 

figures :=

.PHONY: all clean cleanall help

all: $(path).pdf

help:
	@cat README.rst

clean:
	rm -f $(dir)/*.log $(dir)/*.aux $(dir)/*.out $(dir)/*.bbl $(dir)/*.blg $(dir)/*.tmp

cleanall: clean
	rm -rf tmp
	rm -f $(path).pdf

$(path).pdf: $(path).log
	@if [ `grep "Package rerunfilecheck Warning: File" $(path).log | wc -l` != 0 ]; then $(latex); fi
	@if [ `grep "Rerun to get cross-references right." $(path).log | wc -l` != 0 ]; then $(latex); fi
	@if [ `grep "Package natbib Warning: Citation(s) may have changed." $(path).log | wc -l` != 0 ]; then $(latex); fi

$(path).log: $(path).tex $(figures) $(latex_files)
	$(latex)

$(path).bbl: $(path).aux $(dir)/KdV-To-cite.bib
	cd $(dir) && bibtex $(name).aux

$(path).aux: $(path).tex
	$(latex)
