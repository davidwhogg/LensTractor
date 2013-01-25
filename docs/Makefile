.SUFFIXES: .tex .ps .pdf .png .jpg

all: paper1.pdf

clean:
	rm paper1.pdf *.log *.aux

%.pdf: %.tex
	pdflatex $<
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
