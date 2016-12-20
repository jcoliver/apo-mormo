#!/bin/bash

INFILE="apo-ms.md"
OUTFILE="Apodemia-MS.docx"

pandoc $INFILE -f markdown -o $OUTFILE --latex-engine=xelatex --reference-docx apo-pandoc-ref.docx
