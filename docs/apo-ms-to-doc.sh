#!/bin/bash

INFILE="apo-ms.md"
OUTFILE="Apodemia-MS.docx"

pandoc $INFILE -f markdown -o $OUTFILE --reference-docx ~/Documents/pandoc-ref.docx
