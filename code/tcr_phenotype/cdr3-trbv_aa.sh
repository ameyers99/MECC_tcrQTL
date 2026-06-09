#!/bin/sh
length=$1
pos=$2
    
dd=$3
idlist=$(basename -s .tsv "$dd"/*.tsv)

odir=$dd/cdr3-trbv_aa.freq/L${length}_P${pos}
mkdir -p $odir

ofile=$odir/cdr3-trbv_aa.txt
echo "Sample tcr rate" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality and select cdr3-aa length
   cat $dd/$id.tsv |
   cut -f 54,57,60,62 |
   sed -e "1d" |
   awk -v L=$length 'BEGIN{FS="\t"}{
      amino_acid=$1;
      reads=$2;
      productive_frequency=$3;
      v_family=$4;
      if( productive_frequency != "null" && 
          productive_frequency != "" &&
          v_family != "" &&
          v_family != "unresolved" &&
	  length(amino_acid) == L) { 
             print v_family, amino_acid, reads
          }
   }' |
   grep -F -w -v -f filter/missing.vgenes |
   awk -v pos=$pos '{
      v_family=$1;
      amino_acid=$2;
      reads=$3;
      aa=substr(amino_acid, pos, 1);
      v_aa = v_family"+"aa      
      for(i=1; i <=reads; i++){
         print v_aa
      }
   }'  > $odir/tmp.$id #target position amino acids

   #calculate ratio
   Total=$( cat $odir/tmp.$id | wc -l )

   cat $odir/tmp.$id |
   sort | uniq -c |
   awk -v Total=$Total -v id=$id '{print id, $2, $1/Total}'  >> $ofile
   
   rm -f $odir/tmp.$id
   
done

gzip -f $ofile