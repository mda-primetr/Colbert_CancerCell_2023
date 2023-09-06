#!/usr/bin/env perl
#
#
#
use warnings;
use strict;
use DataBrowser qw(browse);
use FAlite;
use Clone 'clone';
#use Sereal;
#use Compress::Zstd qw(compress decompress compress_mt);

my $bacteria = {};
my $hits = {};
my $highest = {};


#my $encoder = Sereal::Encoder->new();
#my $decoder = Sereal::Decoder->new();


#my ($parents, $names, $ranks, $children) = makeTaxonomy($ARGV[3]);

my $origRead = {};
open IN, "lbzcat $ARGV[0] |";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my $read = $parts[0];
#	push @{$origRead->{$read}}, $line;
	$read =~ s/.kraken.*;/;/g;
	my $ref = $parts[1];
	my $evalue = $parts[10];
	if ($evalue > 1e-5) {
		next;
	}
	my $id = $parts[2];
	#unless ($id > .95) {
	#	next;
	#}
	my $bitScore = $parts[11];
	my $rStart = $parts[8];
	my $rEnd = $parts[9];
	if ($rEnd < $rStart) {
		my $temp = $rStart;
		$rStart = $rEnd;
		$rEnd = $temp;
	}
	#$bacteria->{$read}->{$bitScore}->{$ref}->{$rStart} = $rEnd;
	$bacteria->{$read} = 1;
	#$hits->{$ref}->{$read} = 1;
	if (not exists $highest->{$read} or $highest->{$read} < $bitScore) {
		$highest->{$read} = $bitScore;
	}
	

}

close IN;

my $start = scalar(keys %$highest);

open IN, "lbzcat $ARGV[1] |";

my $force = 0;

if (defined $ARGV[2] and $ARGV[2] eq "force") {
	$force = 1;
}


my $removed = {};

while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my $read = $parts[0];
	unless (exists $bacteria->{$read}) {
	#	printLines($origRead, $read);
		next;
	}
	my $bitScore = $parts[11];
	if ($bitScore > $highest->{$read} or $force) {
		delete $highest->{$read};
	#	$removed->{$read} = clone($bacteria->{$read});
	#	$removed->{$read}->{"HIT"} = $line;
	#	foreach my $bit (keys %{$bacteria->{$read}}) {
	#		foreach my $ref (keys %{$bacteria->{$read}->{$bit}}) {
	#			delete $hits->{$ref}->{$read};
	#		}
	#	}
		delete $bacteria->{$read};
	}
	#} else {
	#	printLines($origRead, $read);
	#}
	
}
close IN;



open IN, "lbzcat $ARGV[0] | ";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	if (exists ($bacteria->{$parts[0]})) {
		print "$line\n";
	}
}
close IN;
#foreach my $read (keys %$bacteria) {
#	printLines($origRead, $read);
#}


sub printLines {
	my $cont = $_[0];
	my $read = $_[1];
	foreach my $line (@{$cont->{$read}}) {
		print "$line\n";
	}
	
}



__END__
my $remain = scalar(keys %$highest);
foreach my $read (keys %$highest) {
	#print "read = $read\n";
}

#print "start = $start, remain = $remain\n";

open IN, "$ARGV[2]";

my $fasta_file = new FAlite(\*IN); # or any other filehandle
while(my $entry = $fasta_file->nextEntry) {
	my $def = $entry->def;
	$def =~ s/^>//g;
	if (exists $bacteria->{$def}) {
	#	$bacteria->{$def}->{"SEQ"} = $entry->seq;
	}
	if (exists $removed->{$def}) {
		$removed->{$def}->{"SEQ"} = $entry->seq;
	}

}
close IN;
#browse($bacteria);
#browse($removed);

my $taxaCounts = {};

my $wantedRanks = {
	'superkingdom' => "k__",
	'phylum' => "p__",
	'class' => "c__",
	'order' => "o__",
	'family' => "f__",
	'genus' => "g__",
	'species' => "s__",

};

my @rankOrder = qw(superkingdom phylum class order family genus species);
my $used = {};
foreach my $read (sort keys %{$bacteria}) {
	#print "$read\n";
	$read =~ m/size=(\d+)/;
	my $size = $1;
	my @scores = sort {$b <=> $a} keys %{$bacteria->{$read}};
	my $topScore = $scores[0];
		
	my $goodTaxas = {};
	my $totalRefs = 0;
	foreach my $score (keys %{$bacteria->{$read}}) {
		if ($score < .9 * $topScore) {		
			next;
		}
		foreach my $ref (keys %{$bacteria->{$read}->{$score}}) {
			$ref =~ m/taxId=(\d+)/;
			my $taxaId = $1;
			$goodTaxas->{$taxaId}++;
			$totalRefs++;
			if (exists $used->{$taxaId}) {
				next;
			}
			unless (exists $parents->{$taxaId}) {
				print STDERR "no parent for $taxaId\n";
				next;
			}
			my @chain;
			if (exists $wantedRanks->{$ranks->{$taxaId}}) {
				#push @chain, "$taxaId, $ranks->{$taxaId}";
				push @chain, "$taxaId";
			} else {
			#	push @chain, $taxaId;
			}
			my $node = $taxaId;
			my $escape = 0;
			while ($node != 1) {
				$escape++;
			if ($escape == 100) {
					last;
				}
				if (exists $parents->{$node}) {
					$node = $parents->{$node};
					if (exists $wantedRanks->{$ranks->{$node}}) {
						#push @chain, "$node, $ranks->{$node}";;
					push @chain, "$node";
					} else {
					#	push @chain, "$node";
					}
				} else {
					print STDERR "no parent for $node\n";
				}
			}
			$used->{$taxaId} = \@chain;
			#print "$ref\n";
			#browse(\@chain);
		
		}
	}
	my $votes = {};
	foreach my $taxaId (keys %$goodTaxas) {
		foreach my $node (@{$used->{$taxaId}}) {
			$votes->{$node} += $goodTaxas->{$taxaId};
		}
	}
	my $threshold = .9 * $totalRefs;
	foreach my $vote (keys %$votes) {
		if ($votes->{$vote} >= $threshold) {
			$taxaCounts->{$vote} += $size;
		}
	}
	#print "threshold = $threshold\n";
	#browse($votes);


}
__END__
#browse($used);
#browse($taxaCounts);

my $lines = {};
foreach my $taxaId (keys %$taxaCounts) {
	my $chain = {};
	if (exists $wantedRanks->{$ranks->{$taxaId}}) {
		$chain->{$ranks->{$taxaId}} = $wantedRanks->{$ranks->{$taxaId}} . $names->{$taxaId};
		#push @chain, "$taxaId";
	}
	my $node = $taxaId;
	my $selfRank = $ranks->{$taxaId};
	my $escape = 0;
	while ($node != 1) {
		$escape++; 
		if ($escape == 100) {
			last;  
		}
		if (exists $parents->{$node}) {
			$node = $parents->{$node};
			if (exists $wantedRanks->{$ranks->{$node}}) {
				$chain->{$ranks->{$node}} = $wantedRanks->{$ranks->{$node}}.$names->{$node};
			}
		} else {
			print STDERR "no parent for $node\n";
		}
	}	
	my @name;
#	browse($chain);
	foreach my $rank (@rankOrder) {
		if (exists $chain->{$rank}) {
			push @name, "$chain->{$rank}";
		} else {
			#push @name, $wantedRanks->{$rank} . "unclassified";
			push @name, "unclassified";
		}
		if ($rank eq $selfRank) {
			last;
		}
	}
	#my $name = join ";", @name;
	my $name = join "|", @name;
#	print "$name\t$taxaCounts->{$taxaId}\n";
	$lines->{$name} = $taxaCounts->{$taxaId};	






}

foreach my $taxaName (sort keys %$lines) {
	if ($taxaName =~ m/^Archaea/) {
		next;
	}
	print "$taxaName\t$lines->{$taxaName}\n";
}


sub makeTaxonomy {
	my $in = $_[0];
	open NODE, "$in";
	my $serealData;
	binmode(NODE);
	{
		local $/;
		undef $/;

		$serealData = <NODE>;
	}
	close NODE;
	my $data = $decoder->decode(decompress($serealData));

	my $parents = $data->{'parents'};
	my $ranks = $data->{'ranks'};
	my $names = $data->{'names'};
	my $children = $data->{'children'};
	return($parents, $names, $ranks, $children);

}
