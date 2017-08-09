


TARGETS = exat.pdf

LOGS = $(patsubst %.pdf, %.log, $(TARGETS))
AUXS = $(patsubst %.pdf, %.aux, $(TARGETS))

.PHONY: clean
all: $(TARGETS)

%.pdf: %.tex
	pdflatex $< 
	pdflatex $<

clean: 
	rm $(LOGS) $(AUXS)

