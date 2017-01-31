.PHONY: all clean util sample

PREFIX=.

all: util

install:
ifneq ($(PREFIX),.)
	mkdir -p $(PREFIX)
	cp -drT $(PWD) $(PREFIX)
endif

clean:
	-rm -fr sample-results

util:
	git submodule init
	git submodule update --force
	cd util/bwa && make

update:
	git pull
	git submodule update --recursive --force
	cd util/bwa && make clean && make

documentation:
	-emacs --batch -l ~/.emacs README.org -f org-latex-export-to-latex
	latexmk -pdflatex='lualatex -shell-escape -interaction nonstopmode' -pdf -f README.tex
	latexmk -c

sample:
	bin/proovread --sample --pre sample-results -o --threads 1
	@echo "for results of run see: sample-results"

