#!/bin/sh
length=$1
pos=$2
    
dd=$3
idlist=$(basename -s .tsv "$dd"/*.tsv)

odir=$dd/cdr3_aa.freq/L${length}_P${pos}
mkdir -p $odir

ofile=$odir/cdr3_aa.txt
echo "Sample tcr rate" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality and select cdr3-aa length
   cat $dd/$id.tsv |
   cut -f 54,57,60 |
   sed -e "1d" |
   awk -v L=$length 'BEGIN{FS="\t"}{
      amino_acid=$1;
      reads=$2;
      productive_frequency=$3;
      if( productive_frequency != "null" && 
          productive_frequency != "" &&
          length(amino_acid) == L) { 
             print amino_acid, reads
          }
   }' |
   grep -F -w -v -f filter/missing.vgenes |
   awk -v pos=$pos '{
      amino_acid = $1;
      reads=$2;
      aa=substr(amino_acid, pos, 1);
      for(i=1; i <=reads; i++){
         print aa
      }
   }'  > $odir/tmp.$id #target position amino acids

   #calculate template-weighted ratio
   Total=$( cat $odir/tmp.$id | wc -l )

   cat $odir/tmp.$id |
   sort | uniq -c |
   awk -v Total=$Total -v id=$id '{print id, $2, $1/Total}'  >> $ofile
   
   rm -f $odir/tmp.$id
   
done

gzip -f $ofile


