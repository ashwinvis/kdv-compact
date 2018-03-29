name:=KdV_JSC

all: pdf clean 

pdf:
	latexmk -view=pdf $(name).tex

clean:
	latexmk -c

cleanall:
	latexmk -C
