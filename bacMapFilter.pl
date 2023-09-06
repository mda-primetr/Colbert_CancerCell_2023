#!/usr/bin/env perl
#
#
#
use warnings;
use strict;
use FAlite;
use File::Temp qw/tempdir/;




my $good = {};
my $bad = {};
my $kmers = {};
my $kmerLength = 25;
my $readLetters = {};
my $used = {};
my $seenIn = {};


my $tmpdir = tempdir( CLEANUP => 1 );
open IN, "$ARGV[0]";



my $samples = 0;
while (my $line = <IN>) { 
	$samples++;
	chomp $line;
#	system("scaffoldSplit.pl $line 0 > $tmpdir/test.fa");
	my $cmd = "cat";
	if ($line =~ m/bz2$/) {
		$cmd = "lbzcat";
	}
	open TEMP, "$cmd $line |";
	my $fasta_file = new FAlite(\*TEMP); # or any other filehandle
	while(my $entry = $fasta_file->nextEntry) {	
		my $def = $entry->def;
		$def =~ s/;coords=.*//g;
		$def =~ m/taxId=(\d+)/;
		my $taxId = $1;
		$seenIn->{$taxId}->{$line} = 1;
		$def =~ m/readLetters=(\d+)/;
		my $readLettersLocal = $1;
		unless ($taxId) {
			next;
		}
		unless ($readLetters) {
			next;
		}
		unless (exists $used->{$line}->{$def}) {
			$readLetters->{$taxId} += $readLettersLocal;
			$used->{$line}->{$def} = 1;
		}
		my $seq = uc($entry->seq);
		for (my $i = 0; $i < length($seq) - $kmerLength; $i++) {
			my $kmer = substr($seq, $i, $kmerLength);
			$kmers->{$taxId}->{$kmer} = 1;
		}

	}
	close TEMP;
}
close IN;
print STDERR "Taxa\tReadLetters\tKmersUsed\tSamplesIn\tRatio\tFiltered\n";
foreach my $taxa (keys %$kmers) {
	my $totalReadLetters = $readLetters->{$taxa};
	my $kmersUsed = scalar(keys %{$kmers->{$taxa}});
	my $ratio = $totalReadLetters / $kmersUsed;
	my $inSamples = scalar(keys %{$seenIn->{$taxa}});
	my $filtered = "";
	if (not ($kmersUsed > 50000) and ($ratio > 10 or $kmersUsed < 10000 or ($inSamples < int(0 * $samples) and $kmersUsed < 50000))) {
		print "taxId=$taxa;\n";
		$filtered = "x";
	}
	print STDERR  "$taxa\t$totalReadLetters\t$kmersUsed\t$inSamples\t$ratio\t$filtered\n";
		

}



