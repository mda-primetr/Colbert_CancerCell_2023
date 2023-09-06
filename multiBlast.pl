#!/usr/bin/env perl
#
#
use warnings;
use strict;
use threads;
use File::Temp qw/tempdir/;
use Getopt::Long qw/GetOptionsFromString GetOptions/;
use Time::HiRes qw(time usleep);
use DataBrowser qw(browseErr);

my $pid = $$;

my $devShmTmp = tempdir(PERMS => 0650, DIR => "/dev/shm/", CLEANUP => 1);
system("touch $devShmTmp/pid.$pid");
my $in;
my $out;
my $dbPrefix;
my $threads;
my $params;
my $useZstd = 0;
my $perThread = 0;
GetOptions ("in=s" => \$in, "out=s" => \$out, "dbPrefix=s" => \$dbPrefix, "threads=s" => \$threads, "params=s" => \$params, "useZstd" => \$useZstd, "perThread" => \$perThread);

my $blastxThreads = int($threads / 2) + 1;
my $paramCont = {
	"virmapBlastn" => "blastn -task blastn -word_size 17 -max_target_seqs 100 -outfmt 6 -num_threads 2 -max_hsps 2 -evalue 0.1",
	"virmapMegablast" => "blastn -task megablast -word_size 45 -max_target_seqs 20 -num_threads 2 -outfmt 6 -max_hsps 1 -evalue 0.1 -ungapped",
	"megablast" => "blastn -task megablast -word_size 29 -max_target_seqs 100 -num_threads 1 -outfmt 6 -max_hsps 1 -evalue 0.01 -ungapped",
	"megablastLong" => "blastn -task megablast -word_size 31 -max_target_seqs 10000 -num_threads 1 -outfmt 6 -max_hsps 1 -evalue 0.01",
	"megablastSensitive" => "blastn -task megablast -word_size 21 -max_target_seqs 10000 -num_threads 1 -outfmt 6 -max_hsps 1 -evalue 0.01",
	"megablastMAG" => "blastn -task megablast -word_size 21 -max_target_seqs 10000 -num_threads 2 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop score\" -perc_identity 97",
	"blastn" => "blastn -task blastn -max_target_seqs 1000 -outfmt 6 -word_size 17 -num_threads 1 -max_hsps 10 -evalue 0.01",
	"bacmap" => "blastn -task megablast -word_size 29 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop\" -num_threads 1 -evalue .00001 -max_target_seqs 500 -perc_identity 95 -max_hsps 1 -qcov_hsp_perc 95",
	"bacmapLoose" => "blastn -task megablast -word_size 21 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop\" -num_threads 1 -evalue .00001 -max_target_seqs 500 -perc_identity 85 -max_hsps 1 -qcov_hsp_perc 95",
	"bacmapPri" => "-blastn task megablast -word_size 21 -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop\" -num_threads 1 -evalue .00001 -max_target_seqs 500 -perc_identity 85 -max_hsps 1 -qcov_hsp_perc 95",
	"16S" => "blastn -task megablast -word_size 29 -outfmt 6 -num_threads 1 -evalue 1e-10 -max_target_seqs 10 -perc_identity 97 -max_hsps 1 -qcov_hsp_perc 95",
	"mirna" => "blastn -task blastn-short -word_size 16 -outfmt 6 -num_threads 1 -evalue 0.01 -perc_identity 93.5",
	"probes" => "blastn -task megablast -word_size 15 -outfmt 6 -max_hsps 1 -ungapped -qcov_hsp_perc 90 -perc_identity 85 -num_threads 2",
	"probesLoose" => "blastn -task dc-megablast -word_size 12 -template_length 21 -template_type coding -ungapped -max_hsps 1 -outfmt 6 -qcov_hsp_perc 90 -perc_identity 85 -num_threads 1",
	"probesLoose3" => "blastn -task dc-megablast -word_size 12 -template_length 21 -template_type optimal -ungapped -max_hsps 1 -outfmt 6 -qcov_hsp_perc 90 -perc_identity 85 -num_threads 1",
	"probesLoose2" => "blastn -task blastn-short -word_size 13 -ungapped -max_hsps 1 -outfmt 6 -qcov_hsp_perc 90 -perc_identity 85 -num_threads 1",
	"probesStage2" => "blastn -task megablast -word_size 11 -outfmt 6 -qcov_hsp_perc 90 -perc_identity 85 -num_threads 2",
	"probesStrict" => "blastn -task megablast -word_size 15 -ungapped -max_hsps 1 -outfmt 6 -qcov_hsp_perc 100 -perc_identity 93 -num_threads 2",
	"probesShort" => "blastn -task megablast -word_size 7 -ungapped -max_hsps 1 -outfmt 6 -qcov_hsp_perc 90 -perc_identity 85 -num_threads 2",
	"probesUltraShort" => "blastn -task blastn -word_size 11 -ungapped -max_hsps 1 -num_threads 2 -outfmt 6 -perc_identity 85",
	"probesHost" => "blastn -task blastn -word_size 11 -ungapped -max_hsps 1 -num_threads 2 -outfmt 6 -perc_identity 100",
	"probesBlastx" => "blastx -task blastx -outfmt 6 -qcov_hsp_perc 90 -num_threads $blastxThreads"
	
};

if (-e "$out") {
	system("rm $out");
}
unless (exists $paramCont->{$params}) {
	die "param not defined\n";
}

my $suffix = "ndb";
if ($paramCont->{$params} =~ m/^blastx/) {
	$suffix = "pdb";
	$threads = 2;
	
}
my $fileList = `ls $dbPrefix.*.$suffix 2>/dev/null`;
my @files = split /\n/, $fileList;
my @dbList;
foreach my $file (@files) {
	my $prefix = $file;
	$prefix =~ s/.$suffix$//g;
	push @dbList, $prefix;
}
unless (scalar @dbList) {
	if (-e "$dbPrefix.$suffix") {
		push @dbList, "$dbPrefix";
	}
}
unless (scalar (@dbList)) {
	die "no database found\n";
}
#my $splitBytes = 2097152;
#my $splitBytes = 3145728;
my $splitBytes = 8388608;
#my $splitBytes = 524288;

my $splitThresh = 1.1 * $threads * $splitBytes;
my $inBase = $in;
$inBase =~ s/.*\///g;
#system("cat $tmpPrefix.forBlastnFilter.split.fa | paste - - | shuf | tr \"\\t\" \"\\n\" > $tmpPrefix.forBlastnFilter.split.shuf.fa");
my $tmpdir = tempdir(PERMS => 0650, CLEANUP => 1);
system("chmod g+s $tmpdir");

system("cp $in $tmpdir");

system("cat $tmpdir/$inBase | perl ~/fastaClean.pl | paste - - | shuf | tr \"\\t\" \"\\n\" > $tmpdir/forSplit.fa");

my $splitCount = $threads;
my $inputSize = -s "$tmpdir/forSplit.fa";
if ($perThread or ($inputSize < $splitThresh)) {
	system("fastaSplitter.pl $tmpdir/forSplit.fa $tmpdir/split $threads");
} else {
	my $threadCount = $threads;
	my $customSplit;
	my $i = 1;
	while (1) {
		my $split1 = $inputSize / $threadCount;
		if ($split1 < $splitBytes) {
			$customSplit = int($inputSize / $threadCount) + 1;
			last;
		} else {
			$threadCount += ($threads);	
			$i++;
			if ($i % 2 == 0) {
				$threadCount--;
			}
		}
	}
	#my $split1 = $inputSize / $threads;
	#my $split2 = int($split1 / $splitBytes) + 1;
	#my $parts = $threads * $split2;
	#my $customSplit = int($inputSize / $parts) + 1024;
	print STDERR "split into $customSplit bytes per chunk\n";
	system("cat $tmpdir/forSplit.fa | fastaSplitStream.pl $tmpdir/split $customSplit");
	my $splitList = `ls $tmpdir/split.*.fa`;
	my @splitList = split /\n/, $splitList;
	$splitCount = scalar(@splitList);
	print STDERR "split into $splitCount parts\n";
}
my $splitTimes = {};
for (my $i = 0; $i < $splitCount; $i++) {
	$splitTimes->{$i} = 0;
}
my $count = 0;

foreach my $shard (@dbList) {
	my $t1 = time();
	my $splitLocal = {};
	my $threadToSplit = {};
	my $threadCont = {};
	system("cp $shard.* $devShmTmp");
	my $shardBase = $shard;
	$shardBase =~ s/.*\///g;
	my $threadCount = 0;
	my $resSize = 0;
	#browseErr($splitTimes);
	#my @order;
	my @order = sort {$splitTimes->{$b} <=> $splitTimes->{$a}} keys %$splitTimes;
	#my $initTotal = scalar(@init);
	#print STDERR "total in init: $initTotal\n";
	#while (scalar(@init)) {
	#	push @order, shift @init;
	#	if (scalar (@init)) {
	#		push @order, pop @init;
	#	}
	#}
	my $orderTotal = scalar(@order);
	print STDERR "total in order: $orderTotal\n";
	foreach my $i (@order) {
		if (-e "$tmpdir/split.$i.fa" and -s "$tmpdir/split.$i.fa") {
			#print STDERR "expected out res.split.$count.$i.out\n";
			my $thread = async {
				if ($useZstd) {
					system("$paramCont->{$params} -query $tmpdir/split.$i.fa -db $devShmTmp/$shardBase | zstd -cq > $tmpdir/res.split.$count.$i.out");
				} else {
					system("$paramCont->{$params} -query $tmpdir/split.$i.fa -out $tmpdir/res.split.$count.$i.out -db $devShmTmp/$shardBase");
				}
			};
			$threadCont->{$threadCount} = $thread;
			$splitLocal->{$i} = time();
			$threadToSplit->{$threadCount} = $i;
			$threadCount++;
			while (scalar(keys %$threadCont) == $threads) {
				my @threadKeys = keys %$threadCont;
				foreach my $threadIdx (@threadKeys) {
					my $thread = $threadCont->{$threadIdx};
					if ($thread->is_joinable()) {
						$thread->join();
						my $idx = $threadToSplit->{$threadIdx};
						my $resChunk = -s "$tmpdir/res.split.$count.$idx.out";
						$resSize += $resChunk;
						system("cat $tmpdir/res.split.$count.$idx.out >> $out");
						system("rm $tmpdir/res.split.$count.$idx.out");
						my $finishTime = time();
						my $localElapsed = $finishTime - $splitLocal->{$idx};
						$splitTimes->{$idx} += $localElapsed;
						print STDERR "chunk $idx took $localElapsed for shard $count\n";
						delete($threadCont->{$threadIdx});
					}
				}
				usleep(10000);
			}

		}
	}
	while (scalar(keys %$threadCont)) {
		my @threadKeys = keys %$threadCont;
		foreach my $threadIdx (@threadKeys) {
			my $thread = $threadCont->{$threadIdx};
			if ($thread->is_joinable()) {
				$thread->join();
				my $idx = $threadToSplit->{$threadIdx};
				system("cat $tmpdir/res.split.$count.$idx.out >> $out");
				system("rm $tmpdir/res.split.$count.$idx.out");
				my $finishTime = time();
				my $localElapsed = $finishTime - $splitLocal->{$idx};
				$splitTimes->{$idx} += $localElapsed;
				print STDERR "chunk $idx took $localElapsed for shard $count\n";
				delete($threadCont->{$threadIdx});
			}	
		}
		usleep(10000);
	}
	system("rm $devShmTmp/$shardBase*");
	#system("cat $tmpdir/res.split.$count.*.out >> $out");
	#system("rm $tmpdir/res.split.$count.*.out");
	my $t2 = time();
	my $elapsed = $t2 - $t1;
	print STDERR "finished shard $count in $elapsed seconds\n";
	print STDERR "total results size = $resSize\n";
	$count++;
}

