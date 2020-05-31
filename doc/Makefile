TEX = note.tex
DVI = $(TEX:.tex=.dvi)
PDF = $(TEX:.tex=.pdf)

AUX = $(TEX:.tex=.aux)
LOG = $(TEX:.tex=.log)


.SUFFIXES: .tex .dvi .pdf


all: $(PDF)
$(PDF): $(TEX)
	platex $(TEX)
	dvipdfm $(DVI)

clean: ;
	rm -f $(DVI) $(PDF) $(AUX) $(LOG) missfont.log
