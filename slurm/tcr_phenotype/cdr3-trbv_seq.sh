#!/bin/sh
length_min=$1
length_max=$2

dd=$3
idlist=$(basename -s .tsv "$dd"/*.tsv)

odir=$dd/crd3-trbv_seq.freq/L${length_min}_L${length_max}
mkdir -p "$odir"

ofile=$odir/crd3-trbv_seq.txt
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
      # Read missing tcr values into array miss
      while ((getline line < missfile) > 0) {
        miss[line] = 1
      }
      close(missfile)
    }
    NR > 1 {
      productive_frequency = $60
      amino_acid = $54
      v_family = $62
      reads = $57
      
      trbv_cdr3 = v_family"+"amino_acid

      # Skip missing productive_frequency
      if (productive_frequency == "null" || productive_frequency == "")
        next

      # Skip trbv_cdr3 in missing list
      if (trbv_cdr3 in miss)
        next

      # Check amino_acid length
      if (length(amino_acid) < min || length(amino_acid) > max)
        next

      # Track total productive reads
      total += reads 

      # Skip if unresolved TRBV family
      if (v_family == "" || v_family == "unresolved")
        next

    # Count total reads
    count[trbv_cdr3] += reads

}
    END {
     for (b in count)
      print id, b, count[b] / total
  }

  ' "$tsvfile" >> "$ofile"

done

gzip -f "$ofile"
