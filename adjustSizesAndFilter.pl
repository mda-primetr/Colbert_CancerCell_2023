#!/usr/bin/env perl


use warnings;
use strict;
use Time::HiRes qw(time usleep);
use threads;
use Thread::Queue qw( );
use threads::shared;
use FAlite;
use File::Temp qw/tempdir/;
use RocksDB;
use Digest::MD5 qw(md5_hex);
use Compress::Zstd qw(compress decompress compress_mt);
use DataBrowser qw(browseErr browseNumSortValErr browseNumSortErr browse browseErr);
#use Cpanel::JSON::XS;
use Sereal;


my $q = Thread::Queue->new();
my $q2 = Thread::Queue->new();
my $goQ = Thread::Queue->new();
my $goQ2 = Thread::Queue->new();
my $keysQ = Thread::Queue->new();
my $finished :shared;
$finished = 0;
my $doneEnqueue :shared;
$doneEnqueue = 0;
my $donePileups :shared;
$donePileups = &share({});

my $function = shift @ARGV;
my $quant = 0;
unless ($function eq "build" or $function eq "quant" or $function eq "normCount") {
	die;
}
if ($function eq "quant") {
	$quant = 1;
}
if ($function eq "normCount") {
	$quant = 2;
}
my $pad = shift @ARGV;
unless ($pad =~ m/^\d+$/) {
	die "pad is not a number\n";
}

my $t1 = time();
my $startTime = $t1;
my $THREADS = $ARGV[3];
my $tmpdir = tempdir(PERMS => 0650, DIR => "/dev/shm", CLEANUP => 1);
my $pid = $$;
system("touch $tmpdir/pid.$pid");
my @threads;
my @threads2;
for (my $i = 0; $i < $THREADS; $i++) {
	push @threads, threads->create(\&worker, $quant, $tmpdir, $pad);
}
for (my $i = 0; $i < $THREADS; $i++) {
	push @threads2, threads->create(\&worker2);
}


my $encoder = Sereal::Encoder->new();
my $decoder = Sereal::Decoder->new();
my $dbMASK = RocksDB->new("$tmpdir/MASK", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true'});
my $dbREAD = RocksDB->new("$tmpdir/READ", { create_if_missing => 1, db_log_dir => "/dev/null", keep_log_file_num => 1, max_log_file_size => 1, db_stats_log_interval => 9999999, max_bytes_for_level_base => 100000000, target_file_size_base => 20000000, allow_mmap_reads => 'true', IncreaseParallelism => 1, max_background_compactions => 4, max_background_flushes => 4 });
print STDERR "Done Sereal\n";
my $prevMd5sums = {};
if (-e $ARGV[2]) {
	open IN, "$ARGV[2]";
	while (my $line = <IN>) {
		chomp $line;
		my @parts = split /\t/, $line;
		$prevMd5sums->{$parts[0]} = $parts[1];
	}
	close IN;
}
my $skip = {};
open IN, "$ARGV[1]";
my $fasta_file = new FAlite(\*IN); # or any other filehandle
my $skipped = 0;
my $changed = 0;
my $noExist = 0;
open OUT, ">$ARGV[0].md5sums";
while (my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	my $seq = $entry->seq;
	my $md5sum = md5_hex(substr($seq,$pad, (length($seq) - (-1 * $pad))));
	if (exists $prevMd5sums->{$def} and $md5sum eq $prevMd5sums->{$def}) {
		print ">$def\n";
		print "$seq\n";
		$skip->{$def} = 1;
		$skipped++;
		print OUT "$def\t$md5sum\n";
		next;
	} elsif (exists $prevMd5sums->{$def}) {
		$changed++;
	} else {
		$noExist++;
	}
	#my $length = length($seq);
	#print STDERR "seqlength $length\n";
	my $compressed = compress($seq);
	$dbMASK->put("$def", $compressed);

}
$dbMASK->compact_range;
close IN;
close OUT;
print STDERR "done mask compact\n";
my $sortThreads = $THREADS;
#if ($sortThreads > 8) {
#	$sortThreads = 8;
#}
my $command = "lbzip2 -dc -n$THREADS";
if ($ARGV[0] =~ m/\.zst$/ or $ARGV[0] =~ m/\.zstd/) {
	$command = "zstd -d -c -q"
}
open IN, "$command $ARGV[0] | cut -f1,3,4,6,10 | sort -k2,2 -k3,3n --buffer-size=40G --parallel=$sortThreads --compress-program=zstd |";

print STDERR "done open sam\n";

my $noHits = qr/noHits=1/;
my $lines = {};
my $lineAmount = {};
my $keysAmount :shared;
my $printed = {};
my $skippedLines = 0;
my $queuedLines = 0;
my $referenceParts = {};
my $currentReference = {};
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line, 3;
	my $reference = $parts[1];
	if (not defined $reference or $reference eq "*") {
		next;
	}
	if (exists $skip->{$reference}) {
		$skippedLines++;
		#print STDERR "skipped $reference $line\n";
		next;
	}
	#print STDERR "pushed $reference line $line\n";
	

	unless (exists $currentReference->{$reference}) {
		$currentReference->{$reference} = 0;
		my $localKey = "$reference" . ".$currentReference->{$reference}";
		$referenceParts->{$reference}->{$localKey} = 1;
	}
	
	$queuedLines++;
	push @{$lines->{$reference}}, $line;
	$lineAmount->{$reference}++;
	
	if ($lineAmount->{$reference} % 25000 == 0) {
		my $localKey = "$reference" . ".$currentReference->{$reference}";
		$referenceParts->{$reference}->{$localKey} = 1;
		$dbREAD->put("$localKey", compress_mt($encoder->encode($lines->{$reference}), $THREADS));
		#$dbREAD->put("$localKey", $encoder->encode($lines->{$reference}));
		$lines->{$reference} = [];		
		$currentReference->{$reference}++;

		
	}
	
	



}
foreach my $reference (keys %$lines) {
	if (scalar(@{$lines->{$reference}})) {
		my $localKey = "$reference" . ".$currentReference->{$reference}";
		$referenceParts->{$reference}->{$localKey} = 1;
		#$dbREAD->put("$localKey", $encoder->encode($lines->{$reference}));
		$dbREAD->put("$localKey", compress_mt($encoder->encode($lines->{$reference}), $THREADS));
	}

}

print STDERR "done packing dbRead\n";

$keysAmount = scalar(keys %$lines);
close IN;
for (my $i = 0; $i < $THREADS; $i++) {
	$keysQ->enqueue($keysAmount);
}
close IN;
my $t2 = time();
my $elapsed = $t2 - $t1;
print STDERR "finished reading SAM after $elapsed seconds, skipped $skipped references, $changed references changed content, $noExist references changed names/lengths. $skippedLines skipped lines, $queuedLines queued\n";
$t2 = $t1;

my $expectedEnqueue = 0;
foreach my $hit (sort {$lineAmount->{$b} <=> $lineAmount->{$a}} keys %$lines) {
	my $job = [$hit, $referenceParts->{$hit}];
	$q->enqueue($job);
	$expectedEnqueue++;
	#print STDERR "enqueued $hit\n";
}
#lockPrintStderr("MAIN:done enqueue jobs");
undef($lines);
undef($dbMASK);
undef($dbREAD);

sub worker2 {
	select(STDERR);
	$| = 1;
	my $tid = threads->tid();

	my $start = $goQ2->dequeue();
	
	my $encoder = Sereal::Encoder->new();
	my $decoder = Sereal::Decoder->new();
	my $dbREAD;
	{lock($finished);
	$dbREAD = RocksDB->new("$tmpdir/READ", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	}
	

	
	while (1) {
		#lockPrintStderr("top of loop worker2");

		my $job = $q2->dequeue();

		#lockPrintStderr("after dequeue worker2");
		unless (defined $job) {
			my $peek = $q2->peek();
			if (not defined $peek) {
				last;
			} else {
				redo;
			}
		}

		#lockPrintStderr("worker2 thread: $tid got job!");
		#incoming
		my $jobName = $job->[0];
		my $isContig = $job->[1];
		my $lastPos = $job->[2];
		my $firstPos = $job->[3];
		my $quant = $job->[4];
		#lockPrintStderr("Thread $tid for job $jobName");
		#my $readsUsed;
	
		#local
		my $beginStrip = qr/^(\d+)S/;
		my $sizeQr = qr/size=(\d+)/;
		my $stripQr = qr/S\.(.)/;
		
	
		#stuff to return
		my $count = 0;
		my $readLetters = 0;
		my $pileups = {};
		my $alternate = {};
		my $alternateAdjust = {};
		
	
	
	
		my $lineCont = $decoder->decode(decompress($dbREAD->get("$jobName")));
		#my $lineCont = $decoder->decode($dbREAD->get("$jobName"));
		foreach my $line (@{$lineCont}) {
			my @parts = split /\t/, $line;
			my $cigar = $parts[3];
			my $seq = $parts[4];
			my $start = $parts[2];
			#$readsUsed->{$parts[0]} = 1;
			$parts[0] =~ m/$sizeQr/;
			my $size = $1;
			unless ($size) {
				$size = 1;
			}
			$count += $size;
			if ($quant) {
				$readLetters += length($seq) * $size;
				next;
			}
			#my $beforeParse = time;
			my ($cigarReturn, $Mtotal, $Dtotal) = parseCigar($cigar);
			#my $afterParse = time;
			#$parseCigarTime += ($afterParse - $beforeParse);
			my $length = length($seq);
			if ($Mtotal < $length * .5) {
				next;
			}
			my $negative = 0;
			if ($cigar =~ m/$beginStrip/) {
				$negative = $1;
			}
			$start = $start - $negative;
			my @letters = split //, $seq;
			my @directives = @{$cigarReturn};
			#browse(\@directives);
			my $directivePos = 0;
			my $dist = 0;
			my $Dsize = $size;
			if ($Dtotal > .2 * $length) {
				$Dsize = 0;
			}
			my $consect = 0;
			my $insertAddPos = 0;
			my $insertStr = "";
			for (my $i = 0; $i < scalar(@letters); $i++) {
				my $insertPos = $start + $dist;
				if (not $isContig and $insertPos > $lastPos) {
					$letters[$i] = lc($letters[$i]);
				} 
				if (not $isContig and $insertPos < $firstPos) {
					$letters[$i] = lc($letters[$i]);
				}
				if ($directives[$directivePos] eq "M" or (($directives[$directivePos] eq "S" or $directives[$directivePos] eq "H") and $insertPos < $firstPos)) {
					$pileups->{$insertPos}->{$letters[$i]} += $size;
					$directivePos++;
					$dist++;
					next;
				}
				if (($directives[$directivePos] eq "S" or $directives[$directivePos] eq "H") and $insertPos >= 1) {
					$pileups->{$insertPos}->{"S.$letters[$i]"} += $size;
					$dist++;
					$directivePos++;
					next;
				}
				if ($consect and $directives[$directivePos] eq "I") {
					$insertStr .= $letters[$i];
					$consect++;
					$directivePos++;
					next;
				}
				if ($consect and $directives[$directivePos] ne "I") {
					$alternate->{$insertAddPos}->{$insertStr} += $size;
					$alternateAdjust->{$insertAddPos} += $size;
					#$dist += length($insertStr);		
					#$insertPos += length($insertStr);
					$insertStr = "";
					$consect = 0;
						
				}
				if ($directives[$directivePos] eq "D") {
					if ($Dsize) {
						$pileups->{$insertPos}->{"D"} += $Dsize;
					}
					$directivePos++;
					$dist++;
					redo;
				}
				if (not $consect and $directives[$directivePos] eq "I") {
					$directivePos++;
					$insertAddPos = $insertPos;
					$consect = 1;
					$insertStr = "$letters[$i]";
					next;
				}
				if ($directives[$directivePos] eq "X.M") {
					$dist++;
					$directivePos++;
					next;
				}
				if ($directives[$directivePos] eq "X.D") {
					$dist++;
					$directivePos++;
					redo;
				}
				if ($directives[$directivePos] eq "X.I") {
					$directivePos++;
					next;
				}
			
			}
			
		}
		#lockPrintStderr("Threads $tid on job $jobName finshed");
		
		my $return = {};
		$return->{'count'} = $count;
		$return->{'readLetters'} = $readLetters;
		$return->{'pileups'} = $pileups;
		$return->{'alternate'} = $alternate;
		$return->{'alternateAdjust'} = $alternateAdjust; 
		#my $returnData = $encoder->encode($return);
		lock($donePileups);		
		$donePileups->{$jobName} = $encoder->encode($return);
		#lockPrintStderr("Threads $tid on job $jobName pushed results back up");
	
	}
	
}	
	
	
sub worker {
	select(STDERR);
	$| = 1;
	my $tid = threads->tid();	
	my $quat = $_[0];
	my $tmpdir = $_[1];
	my $padded = $_[2];
	my $readRocks = $_[3];
	my $encoder = Sereal::Encoder->new();
	my $decoder = Sereal::Decoder->new();
	my $start = $goQ->dequeue();
	my $keysAmount = $keysQ->dequeue();
	my $dbMASK;
	
	my $dbCHECK;
	{lock($finished);
	$dbMASK = RocksDB->new("$tmpdir/MASK", { read_only => 1, db_log_dir => '/dev/null', keep_log_file_num => 1, allow_mmap_reads => 'true'}) or die "can't open DB\n";
	}
	my $beginStrip = qr/^(\d+)S/;
	my $sizeQr = qr/size=(\d+)/;
	my $stripQr = qr/S\.(.)/;
	my $return = {};
	while (1) {
		#lockPrintStderr("worker1 top of loop");
		my $isContig = 0;
		my $job = $q->dequeue();
		#lockPrintStderr("worker1 after dequeue");
		unless (defined $job) {
			my $peek = $q->peek();
			if (not defined $peek) {
				last;
			} else {
				redo;
			}
		}
		#lockPrintStderr("worker1 thread: $tid got job!");
		my $t1 = time;
		
		my $reference = $job->[0];
		my $referenceParts = $job->[1];
		my $afterLineContTime = time;
		my $lineContDecompressTime = $afterLineContTime - $t1;
		if ($reference =~ m/^contig/) {
			$isContig = 1;
		}
		my $pileups = {};
		my $count = 0;
		my $readLetters = 0;
		my $alternate = {};
		my $alternateAdjust = {};
		#unless (scalar(@{$lineCont})) {
		#	lockPrintStderr("$reference has no reads");
		#	next;
		#}
		my $printed = 0;
	#	print STDERR "got $reference\n";
		my $origScaffold = decompress($dbMASK->get($reference));
		my $scaffold = "a" . $origScaffold;
		my @scafLetters = split //, $scaffold;
		my $caseMask = [];
		foreach my $letter (@scafLetters) {
			if ($isContig) {
				push @$caseMask, 0;
				next;
			}
			if ($letter =~ m/[acgt]/) {
				push @$caseMask, 1;
			} elsif ($letter eq "n") {
				push @$caseMask, 2;
			} else {
				push @$caseMask, 0;
			}
		}
		unless(scalar(@{$caseMask})) {
			browseErr($caseMask);
			die "casemask is empty\n";
		}
		my $lastPos = scalar(@{$caseMask}) - $padded;
		my $firstPos = $padded + 1;
		my $time1 = time;
		my $parseCigarTime = 0;
		#my $readsUsed = {};
		foreach my $subPart (keys %$referenceParts) {
			my $job2 = [$subPart, $isContig, $lastPos, $firstPos, $quant];
			lock($q2);			
			$q2->enqueue($job2);
		}
			
		#lockPrintStderr("worker thread $tid: before infinite loop");	
		while (1) {
			foreach my $subPart (keys %$referenceParts) {
				my $returnHash;
				{lock($donePileups);
				if (exists $donePileups->{$subPart}) {
					$returnHash = $decoder->decode($donePileups->{$subPart});		
					delete($donePileups->{$subPart});
				}
				}
				if ($returnHash) {
					delete($referenceParts->{$subPart});
					$count += $returnHash->{'count'};
					$readLetters += $returnHash->{'readLetters'};
					foreach my $pos (keys %{$returnHash->{'pileups'}}) {
						foreach my $event (keys %{$returnHash->{'pileups'}->{$pos}}) {
							$pileups->{$pos}->{$event} += $returnHash->{'pileups'}->{$pos}->{$event};
						}
					}
					foreach my $pos (keys %{$returnHash->{'alternate'}}) {
						foreach my $event (keys %{$returnHash->{'alternate'}->{$pos}}) {
							$alternate->{$pos}->{$event} += $returnHash->{'alternate'}->{$pos}->{$event};
						}
					}
					foreach my $pos (keys %{$returnHash->{'alternateAdjust'}}) {
						$alternateAdjust->{$pos} += $returnHash->{'alternateAdjust'}->{$pos};
					}		
				}			
			}
			if (scalar keys %$referenceParts == 0) {
				last;
			} else {
				#lockPrintStderr("working1 waiting for results refparts not empty");
				unless($quant) {
					usleep(100);
				}
			}
		}



		#get results from threading
		#also enqueue


		my $time2 = time;
		my $pileTime = $time2 - $time1;
		$time1 = $time2;
		#undef($lineCont);
		my @sequence;
		if (not scalar(keys %$pileups) and not scalar(keys %{$alternate}) and not $quant) {
			#lockPrintStderr("$reference has no valid hits");
			next;
		}
		if ($quant) {
			$pileups->{0} = 1;
		}
		my @positions = sort {$a <=> $b} keys %{$pileups};
		my $start = $positions[0];
		my $stop = $positions[$#positions];
		if (scalar(keys %$alternate) and not $quant) {
			my @altPos = sort {$a <=> $b} keys %{$alternate};
			if ($altPos[0] < $start) {
				$start = $altPos[0];
			}
			if ($altPos[$#altPos] > $stop) {
				$stop = $altPos[$#altPos]
			}
		}
		for (my $i = $start; $i <= $stop; $i++) {
			if ($quant) {
				last
			}
			if (not exists $pileups->{$i} and not exists $alternate->{$i}) {
				if ($caseMask->[$i] == 2) {
					push @sequence, "n";
				} else {
					push @sequence, "N";
				}
			} else {
				my @order;
				my $noOrig = 0;
				if (exists $pileups->{$i}) {
					my @keys = sort {$pileups->{$i}->{$b} <=> $pileups->{$i}->{$a}} keys %{$pileups->{$i}};
					if (scalar(@keys) > 1) {
						my @toDelete;
						foreach my $letter (@keys) {
							if ($letter =~ m/S\./) {
								push @toDelete, $letter;
							}
						}
						if (scalar(@toDelete) >= 1 and scalar(@toDelete) != scalar(@keys)) {
							foreach my $key (@toDelete) {
								delete $pileups->{$i}->{$key};
							}
							@keys = keys %{$pileups->{$i}};
						} 
					}
					@order = sort {$pileups->{$i}->{$b} <=> $pileups->{$i}->{$a}} @keys;
				} else {
					$noOrig = 1;
					@order = qw(X);
					$pileups->{$i}->{$order[0]} = -1000000;
				}
				if (exists $alternate->{$i}) {
					my @altOrder = sort {$alternate->{$i}->{$b} <=> $alternate->{$i}->{$a} || $b cmp $a} keys %{$alternate->{$i}};
					my $altInsert = $altOrder[0];
					if ($alternate->{$i}->{$altInsert} >= ($pileups->{$i}->{$order[0]} - $alternateAdjust->{$i})) {
						if (defined $altOrder[1] and $alternate->{$i}->{$altInsert} == $alternate->{$i}->{$altOrder[1]}) {
							$altInsert = "n" x length($altInsert);;
						}
						if ($caseMask->[$i] == 1) {
							$altInsert = lc($altInsert);
						}
						push @sequence, $altInsert;
					}
					if ($noOrig) {
						delete $pileups->{$i};
					}
				}
				if (exists $pileups->{$i}) {
					if (defined $order[1] and $pileups->{$i}->{$order[0]} == $pileups->{$i}->{$order[1]}) {
						push @sequence, "n";
					} else {
						unless ($order[0] eq "D") {
							if ($order[0] =~ m/$stripQr/) {
								$order[0] = $1;
							}
							unless($order[0]) {
								lockPrintStderr("WARNING: putting nothing: $order[0] in sequence for $reference at pos $i");
							}
							if (not defined $caseMask->[$i] or $caseMask->[$i] == 1) {
								$order[0] = lc($order[0]);
							}
							push @sequence, $order[0];
						}
					}
				}
			}	
		}
		$time2 = time;
		my $rebuildTime = $time2 - $time1;
		my $outSeq = "";
		unless ($quant) {
			$outSeq = join "", @sequence;
		} else {
			$outSeq = $origScaffold;
		}
		$pileups = {};
		$alternate = {};
		my $sizeDiff = 0;
		my $lengthDiff = 0;
		my $sizeErrStr = "";
		my $beforeLength;
		my $newLength = length($outSeq);
		if ($reference =~ m/length=(\d+)/) {
			$beforeLength = $1;
		} else {
			$beforeLength = $newLength;
		}
		if ($reference =~ m/length=\d+/) {
			$reference =~ s/length=\d+/length=$newLength/g;
		} else {
			$reference .= ";length=$newLength";
		}
		$lengthDiff = $newLength / $beforeLength;
		$sizeErrStr = "(lengthDiff = $lengthDiff)";
		if ($quant) {
			my $countName = "size";
			my $letterName = "readLetters";
			if ($quant == 2) {
				$countName = "normSize";
				$letterName = "normLetters";
			}
			my $beforeSize;
			if ($reference =~ m/;$countName=(\d+)/g) {
				$beforeSize = $1;
			}
			unless ($beforeSize) {
				$beforeSize = $count;
			}
			$sizeDiff = $count / $beforeSize;
			if ($reference =~ m/$countName=\d+/) {
				$reference =~ s/$countName=\d+/$countName=$count/g;
			} else {
				$reference .= ";$countName=$count";
			}
			if ($reference =~ m/$letterName=\d+/) {
				$reference =~ s/$letterName=\d+/$letterName=$readLetters/g;
			} else {
				$reference .= ";$letterName=$readLetters";
			}
			$sizeErrStr = "(sizeDiff = $sizeDiff; lengthDiff = $lengthDiff)";
		}
		$return->{$reference} = $outSeq;
		{lock($finished);
		$finished++;
		my $t2 = time;
		my $elapsed = $t2 - $t1;
		#if ($elapsed > 5) {
			lockPrintStderr("finished on thread $tid ($finished / $keysAmount) time: $elapsed seconds $reference $sizeErrStr, pileTime = $pileTime, rebuildTime = $rebuildTime, parseCigarTime = $parseCigarTime, lineContDecompressTime = $lineContDecompressTime");;
			#select->flush();
		#}
		}
	}
	return($return);
}
for (my $i = 0; $i < $THREADS; $i++) {
	$q->enqueue(undef);
}
for (my $i = 0; $i < $THREADS; $i++) {	
	$goQ->enqueue("GO");
}

for (my $i = 0; $i < $THREADS; $i++) {
	$goQ2->enqueue("GO");
}






$t2 = time();
$elapsed = $t2 - $t1;
$t1 = $t2;
my $seqs = {};
#lockPrintStderr("MAIN:before worker join");
foreach my $thr (@threads) {
	my $tid = $thr->tid();
	my $return = $thr->join();
	foreach my $def (keys %$return) {
		$seqs->{$def} = $return->{$def};
	}
}
#lockPrintStderr("MAIN:after join worker");

for (my $i = 0; $i < $THREADS; $i++) {
	$q2->enqueue(undef);
}
#lockPrintStderr("after worker2 enqueue undef");
foreach my $thr (@threads2) {
	$thr->join();
}
#lockPrintStderr("after worker2 join");
open OUT, ">$ARGV[0].md5sums";
foreach my $def (sort keys %$seqs) {
	unless (exists $seqs->{$def} and defined $seqs->{$def}) {
		next;
	}
	print ">$def\n";
	print "$seqs->{$def}\n";
	my $md5sum = md5_hex($seqs->{$def});
	print OUT "$def\t$md5sum\n";
}


my $endTime = time();
my $overallTime = $endTime - $startTime;
#print STDERR "Overall pileup time $overallTime seconds\n";

sub parseCigar {
	my $cigar = $_[0];
	my @directives;
	my @cigNum = split /[A-Z]/, $cigar;
	my @cigLet = split /\d+/, $cigar;
	shift @cigLet;
	my $Mtotal = 0;
	my $Dtotal = 0;
	my $largestD = 0;
	my $largestDi = 0;
	for (my $i = 0; $i < scalar(@cigNum); $i++) {
		if ($cigLet[$i] eq "M") {
			$Mtotal += $cigNum[$i];
		}
		if ($cigLet[$i] eq "D") {
			$Dtotal += $cigNum[$i];
			if ($cigNum[$i] > $largestD) {
				$largestD = $cigNum[$i];
				$largestDi = $i;
			}
		}
		for (my $j = 0; $j < $cigNum[$i]; $j++) {
			push @directives, $cigLet[$i];
		}
	} 
	if ($largestD > .5 * $Mtotal) {
		@directives = ();
		$Mtotal = 0;
		$Dtotal = 0;
		my $leftM = 0;
		my $rightM = 0;
		for (my $i = 0; $i < scalar(@cigNum); $i++) {
			if ($i < $largestDi and $cigLet[$i] eq "M") {
				$leftM += $cigNum[$i];
			}
			if ($i > $largestDi and $cigLet[$i] eq "M") {
				$rightM += $cigNum[$i];
			}
		}
		for (my $i = 0; $i < scalar(@cigNum); $i++) {
			my $letter = $cigLet[$i];
			if ($leftM > $rightM and $i > $largestDi) {
				$letter = "X.$letter";
			} elsif ($rightM > $leftM and $i < $largestDi) {
				$letter = "X.$letter";
			} elsif ($leftM == $rightM or $i == $largestDi) {
				$letter = "X.$letter";
			}
			for (my $j = 0; $j < $cigNum[$i]; $j++) {
				if ($letter eq "M") {
					$Mtotal++;
				}
				if ($letter eq "D") {
					$Dtotal++;
				}
				push @directives, $letter;
			}
		}
	}
	my $directive = \@directives;
	return($directive, $Mtotal, $Dtotal);


}
sub lockPrintStderr {
	lock($finished);
	my $tid = threads->tid();
	my $printout = $_[0];
	print STDERR "THREAD $tid - $printout\n";
}
