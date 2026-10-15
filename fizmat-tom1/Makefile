# Сборка переработанного учебника (нужен XeLaTeX и latexmk).
#   make glava2   — только глава 2
#   make all      — весь том (готовые главы)
#   make clean    — удалить промежуточные файлы
LATEXMK = latexmk -xelatex -interaction=nonstopmode -file-line-error -outdir=build

.PHONY: all glava1 glava2 glava3 clean

all:
	$(LATEXMK) main.tex && cp build/main.pdf pdf/fizmat-tom1.pdf

glava1:
	$(LATEXMK) glava1.tex && cp build/glava1.pdf pdf/glava1-kristallografiya.pdf

glava2:
	$(LATEXMK) glava2.tex && cp build/glava2.pdf pdf/glava2-defekty.pdf

glava3:
	$(LATEXMK) glava3.tex && cp build/glava3.pdf pdf/glava3-fizika-tverdogo-tela.pdf

clean:
	latexmk -C -outdir=build
	rm -rf build
