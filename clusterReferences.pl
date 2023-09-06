#!/usr/bin/env perl


use warnings;
use strict;
use DataBrowser qw(browseErr browse);

my $in = $ARGV[0];
my $data = {};

my $prot2NucMap = {};
my $nuc2FullName = {};
my $nucHits = {};
open IN, "lbzip2 -d -c -n3 $ARGV[1] | cut -f1,3 | grep -v ^@ |";

while (<IN>) {
	my $start = index($_, "\t");
	my $read = substr($_, 0, $start);
	my $genome = substr($_, $start + 1, -1);
	$data->{$genome}->{$read} = 1;
	$nucHits->{$genome}++;
	
	$start = index($genome, 'b');
	my $end = index($genome, '|', $start + 1);
	$start = index($genome, '|', $end + 1);
	
	my $acc = substr($genome, $end + 1, ($start - ($end + 1)));
	$nuc2FullName->{$acc} = $genome;
	
	
}
close IN;


my $full2parts = {};
open IN, "lbzip2 -d -c -n3 $in | cut -f1,3 | grep -v ^@ |";
while (<IN>) {
	my $start = index($_, "\t");
	my $read = substr($_, 0, $start);
	my $geneHit = substr($_, $start + 1, -1);
	$start = index($geneHit, '|');
	my $end = index($geneHit, '|', $start + 1);
	$start = index($geneHit, '|', $end + 1);	
	my $nucAcc = substr($geneHit, $end + 1, ($start - ($end + 1)));

	unless (exists $nuc2FullName->{$nucAcc}) {
		$start = rindex($geneHit, "|");
		$end = "gi|NoNucHits|gb|$nucAcc|" . substr($geneHit, $start + 1, (index($geneHit, ";") - ($start + 1))) . ";taxId=" . substr($geneHit, rindex($geneHit, ";taxId=") + 7);
		$nuc2FullName->{$nucAcc} = $end;
		$nucHits->{$end} = 0;
	} else {
		$end = $nuc2FullName->{$nucAcc};
	}
	unless (exists $full2parts->{$end}->{$geneHit}) {
		$full2parts->{$end}->{$geneHit} = 1;
	}
	$data->{$end}->{$read} = 1;
}
close IN;

my @references = sort {scalar(keys %{$data->{$b}}) <=> scalar(keys %{$data->{$a}}) || $nucHits->{$b} <=> $nucHits->{$a} || $a cmp $b} keys %{$data};
my %good;
my $clusterRadius = .96;
my $blackList = {};
my $chimeras = {};
my $whiteList = {};
my $e = 2.7182818284590452353602874713527;
my $gr = 1.6180339887498948482045868343656;
my $pi = 3.1415926535897932384626433832795;
my $ratio = $e;
foreach my $reference (@references) {
	my $nonInto = {};
	my $blackListCount = 0;
	my $count = 0;
	my $intersect = 0;
	foreach my $read (keys %{$data->{$reference}}) {
		if (exists $blackList->{$read}) {
			$blackListCount++;
			next;
		}
		$count++;
		if (exists $whiteList->{$read}) {
			$intersect++;
		} else {
			$nonInto->{$read} = 1;
		}
	}
	unless ($count) {
		next;
	}
	my $goodRefReadAmount = scalar(keys %$whiteList);
	unless ($goodRefReadAmount) {
		$goodRefReadAmount = $count;
	}
	my $adjustment = ((log($goodRefReadAmount/$count)/log($ratio)) / 100) * $ratio;
	my $intersectionPerc = ($intersect/$count);
	my $overallReadCount = scalar(keys %{$data->{$reference}});
	my $blackListIntersection = ($blackListCount/$overallReadCount);
	if ($blackListIntersection > $intersectionPerc) {
		$intersectionPerc = $blackListIntersection;
	}
	my $explainPerc = $intersect / $overallReadCount;
	my $radius = $clusterRadius - $adjustment;
	if ($intersectionPerc > $radius) {
		print STDERR "Absorbed because ($intersectionPerc > ($radius) [adjustment: $adjustment] or ($blackListIntersection) ($blackListCount/$overallReadCount) > ($radius))\n";
		print STDERR "Rejected: $reference\n$intersect explained by more than one previously accepted centroid; OverallReads: $overallReadCount; perc: $explainPerc\n\n";
		foreach my $read (keys %$nonInto) {
			$blackList->{$read} = 1;
		}
		next;
	} else {
		print STDERR "Not Absorbed because ($intersectionPerc > ($radius) [adjustment: $adjustment] or ($blackListIntersection) ($blackListCount/$overallReadCount) > ($radius))\n";
		print STDERR "Accepted: $reference\n$intersect explained by more than one previously accepted centroid; OverallReads: $overallReadCount; perc: $explainPerc\n\n";
		print STDERR "Added $reference as a new centroid\n";
		$good{$reference} = 1;
		foreach my $read (keys %{$data->{$reference}}) {
			$whiteList->{$read} = 1;
		}
	}
}
foreach my $reference (sort keys %good) {
	print "Parent:$reference\n";
	if (exists $full2parts->{$reference}) {
		foreach my $acc (sort keys %{$full2parts->{$reference}}) {
			print "Child:$acc\n";
		}	
	}
}	
