#!/usr/bin/env bash
awk '
  BEGIN{ OFS="\t" }
  FNR==1{
    cell = FILENAME
    sub(/^.*\//,"",cell)
    sub(/_report\.tsv$/,"",cell)
  }
  /^#count/ && !printed++{
    sub(/^#count/,"count")
    print $0 "\tcell"
    next
  }
  $0 !~ /^#/ { print $0 "\t" cell }
' *report.tsv > combined_trust4.tsv
echo "Wrote combined_trust4.tsv"
