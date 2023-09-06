#!/usr/bin/env perl
#
#
use warnings;
use strict;
use DataBrowser;
use FAlite;
use threads;
use threads::shared;
use Thread::Queue qw( );
use Compress::Zstd qw(compress decompress compress_mt);
use File::Temp qw/tempdir/;
use Sereal;



my $tmpdir = tempdir(PERMS => 0650, DIR => "/dev/shm", CLEANUP => 1);
my $pid = $$;
my $lock :shared;

system("touch $tmpdir/pid.$pid");
my $THREADS = $ARGV[2];
unless (defined $THREADS) {
	$THREADS = 8;
}

my $q = Thread::Queue->new();
my $startQ = Thread::Queue->new();
my $encoder = Sereal::Encoder->new();

my $cleanup = threads->create(\&cleaner, $tmpdir);

my @threads;
for (my $i = 0; $i < $THREADS; $i++) {
	$threads[$i] = threads->create(\&worker, $tmpdir)
}


my $reads = {};
my $limits = {};
open IN, "lbzcat $ARGV[0] | cut -f1,12 |";
#while (my $line = <IN>) {
while(<IN>) {
	#chomp $line;
	#my @parts = split /\t/, $line;
	my $start = index($_, "\t");
	my $read = substr($_, 0, $start);
	my $score = substr($_, $start + 1, -1);
	
	#print "$read\n";
	#print "$genome\n";
	if (not exists $reads->{$read} or (exists $reads->{$read} and $limits->{$read} < $score)) { 
		$reads->{$read} = 1;
#	if (not exists $limits->{$read} or $limits->{$read} < $score) {
		$limits->{$read} = $score;
	}

#	if (not exists $reads->{$parts[0]} or (exists $reads->{$parts[0]} and $limits->{$parts[0]} < $parts[1])) {
#		$reads->{$parts[0]} = 1;
#		$limits->{$parts[0]} = $parts[1];
#	}

}
close IN;
my $data = $encoder->encode($limits);
open OUT, ">$tmpdir/limits";
print OUT "$data";
close OUT;

open IN, "$ARGV[1]";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
while(my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	if (exists $reads->{$def}) {
		$reads->{$def} = length($entry->seq);
	}

}
close IN;

my $readData = $encoder->encode($reads);
open OUT, ">$tmpdir/reads";
print OUT "$readData";
close OUT;


for (my $i = 0; $i < $THREADS; $i++) {
	$startQ->enqueue("GO");
}

open IN, "lbzcat $ARGV[0] |";
my $count = 0;
my @lines;
my $jobCount = 0;
open OUT, ">$tmpdir/$jobCount";
while (<IN>) {
	$count++;
	print OUT "$_";
	if ($count == 1000000) {
		close OUT;
		$q->enqueue($jobCount);
		$jobCount++;
		$count = 0;
		open OUT, ">$tmpdir/$jobCount";
		while (1) {
			if ($q->pending > 10) {
				sleep (1);
			} else {
				last;
			}
		}
	}
}
close IN;
close OUT;
$q->enqueue($jobCount);
for (my $i = 0; $i < $THREADS; $i++) {
	$q->enqueue(undef);
}

foreach my $thr (@threads) {
	$thr->join();
}
system("touch $tmpdir/overallFinished");

$cleanup->join();

sub cleaner {
	my $tmpdir = $_[0];
	while (1) {
		sleep(1);
		#print STDERR "checking dir\n";
		opendir my $dh, "$tmpdir";
		my @files = readdir $dh;
		closedir $dh;
		my $noFinished = 1;
		my $overallFinished = 0;
		foreach my $file (@files) {
			if ($file =~ m/overallFinished/) {
				$overallFinished = 1;
			}	
			if ($file =~ m/(\d+)\.finished$/) {
				print STDERR "got $file\n";
				$noFinished = 0;
				open my $fh, "$tmpdir/$1.result";
				while (<$fh>) {
					print "$_";
				}
				close $fh;
				system("rm $tmpdir/$1.result; rm $tmpdir/$1.finished");
			}
		}
		if ($noFinished and $overallFinished) {
			last;
		}
			

		
	}
	return();

}



sub worker {
	my $tmpdir = $_[0];
	my $tid = threads->tid();
	$startQ->dequeue();
	my $decoder = Sereal::Decoder->new();
	my $limitData;
	my $limitfh;
	open $limitfh, "$tmpdir/limits";
	binmode($limitfh);
	{
		local $/;
		undef $/;

		$limitData = <$limitfh>;
	}
	my $limits = $decoder->decode($limitData);
	close $limitfh;
	my $readfh;
	open $readfh, "$tmpdir/reads";
	my $readData;
	binmode($readfh);
	{
		local $/;
		undef $/;
		$readData = <$readfh>;
	}
	close $readfh;
	my $reads = $decoder->decode($readData);
	while (1) {
		my $job = $q->dequeue();
		if (not defined $job) {
			my $peek = $q->peek();
			if (not defined $peek) {
				last;
			} else {
				redo
			}
		}
		my $infh;
		open $infh, "$tmpdir/$job";
		my $outfh;
		open $outfh, ">$tmpdir/$job.result";
		#my @linesOut;
		while (my $line = <$infh>) {
			chomp $line;
			my @parts = split /\t/, $line;
			my $bitScore = $parts[11];
			if ($bitScore < .85 * $limits->{$parts[0]}) {
				next;
			}
			my $btop = $parts[12];
			#print "$btop\n";
			my @btop = $btop =~ /\d+|\D+/g;
			#browse(\@btop);
			my $cigar = btopToSam(\@btop);
			my $qStart = $parts[6];
			my $qEnd = $parts[7];
			my $length = $reads->{$parts[0]};
			if ($qEnd < $length) {
				my $endPad = $length - $qEnd;
				push @$cigar, $endPad . "S";
			}
			if ($qStart > 1) {
				my $startPad = $qStart - 1;
				unshift @$cigar, $startPad . "S";
			}
			my $refStart = $parts[8];
			my $refEnd = $parts[9];
			my $errors = $parts[4] + $parts[5];
			my $reverse = 0;
			if ($refEnd < $refStart) {
				my @cigarPrint = reverse(@$cigar);
				my $cigarPrint = join "", @cigarPrint;
				print $outfh "$parts[0]\t16\t$parts[1]\t$refEnd\t*\t$cigarPrint\tNM:i:$errors\n";
				#push @linesOut, "$parts[0]\t16\t$parts[1]\t$refEnd\t*\t$cigarPrint\tNM:i:$errors";
			} else {
				my $cigarPrint = join "", @$cigar;
				print $outfh "$parts[0]\t0\t$parts[1]\t$refStart\t*\t$cigarPrint\tNM:i:$errors\n";
				#push @linesOut, "$parts[0]\t0\t$parts[1]\t$refStart\t*\t$cigarPrint\tNM:i:$errors";
			}
		}
		close $infh;
		close $outfh;
		system("touch $tmpdir/$job.finished");
		system("rm $tmpdir/$job");
	}	
	return();
}
sub btopToSam {
	my @in = @{$_[0]};
	my @cigarCount;
	my @cigarElem;
	foreach my $elem (@in) {
		if ($elem =~ m/^(\d+)$/) {
			push @cigarCount, $1;
			push @cigarElem, "M";
		} elsif (not $elem =~ m/\-/) {
			push @cigarCount, length($elem) / 2;
			push @cigarElem, "M";
		} else {
			my @subchunks = ( $elem =~ m/../g );
	#		print "subchunks\n";
	#		browse(\@subchunks);
			foreach my $subelem (@subchunks) {
				if (not $subelem =~ m/\-/) {
					push @cigarCount, "1";
					push @cigarElem, "M";
				} elsif ($subelem =~ m/^\-/) {
					push @cigarCount, "1";
					push @cigarElem, "D";
				} elsif ($subelem =~ m/\-$/) {
					push @cigarCount, "1";
					push @cigarElem, "I";
				}
			}
		}
	}
	#browse(\@cigarCount);
	#browse(\@cigarElem);
	my @cigar;
	my $currentCount = 0;
	my $currentElem = "";
	for (my $i = 0; $i < scalar(@cigarCount); $i++) {
		unless ($currentElem) {
			$currentCount = $cigarCount[$i];
			$currentElem = $cigarElem[$i];
			next;
		}
		if ($cigarElem[$i] eq $currentElem) {
			$currentCount += $cigarCount[$i];
		} else {
			my $combined = $currentCount . $currentElem;
			push @cigar, $combined;
			$currentCount = $cigarCount[$i];
			$currentElem = $cigarElem[$i];
		}
	}
	my $combined = $currentCount . $currentElem;
	push @cigar, $combined;
	#browse(\@cigar);
	return (\@cigar);
	
}
