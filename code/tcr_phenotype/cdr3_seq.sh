#!/bin/sh
length_min=$1
length_max=$2

dd=$3
idlist=$(basename -s .tsv "$dd"/*.tsv)

odir=$dd/cdr3_seq.freq/L${length_min}_L${length_max}
mkdir -p "$odir"

ofile=$odir/cdr3_seq.txt
echo "Sample tcr rate" > "$ofile"

for id in $idlist; do
  echo "Processing sample $id..."

  tsvfile="$dd/$id.tsv"
  if [ ! -f "$tsvfile" ]; then
    echo "Warning: File $tsvfile not found, skipping."
    continue
  fi

awk -v min=$length_min -v max=$length_max -v missfile="filter/missing.vgenes" -v id="$id" '
  BEGIN {
    FS="\t"
    while ((getline line < missfile) > 0) {
      miss[line] = 1
    }
    close(missfile)
  }
  NR > 1 {
    cdr3 = $54
    reads = $57
    productive_frequency = $60

    if (productive_frequency == "null" || productive_frequency == "")
      next

    if (cdr3 in miss)
      next

    if (length(cdr3) < min || length(cdr3) > max)
      next
    
    # Count total reads
    count[cdr3] += reads

    # Track total frequency
    total += reads
    }
    END {
     for (b in count)
      print id, b, count[b] / total
  }
' "$tsvfile" >> "$ofile"

done

gzip -f "$ofile"
