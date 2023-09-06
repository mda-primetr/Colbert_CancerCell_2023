#!/usr/bin/env perl


#Run in directory with all the bacterial reconstruction results
find | grep "quant.fa$" | parallel -I {} 'lbzip2 -n1 {}' 
find ./* | grep "quant.fa.bz2$" | parallel -I {} 'echo "scaffoldSplit.pl {} 0 | lbzip2 -n1 -c > {}.split.fa.bz2"' | parallel 
find ./* | grep split.fa.bz2 | xargs readlink -e > ForFilterList
bacMapFilter.pl ForFilterList > Filter.out 2>Filter.err 
find | grep "quant.fa.bz2$" | parallel -I {} 'echo "perl ~/bacMapGetFilterCounts.pl Filter.out {} > {}.counts.txt"' | parallel
find ./* | grep quant.fa.bz2$ | parallel -I {} 'bzcat {} | grep "^>" | sed -re "s:^:{}\t:g"' > OverallFinal.txt
cat OverallFinal.txt| grep -v -f Filter.out > Filtered.final.txt

