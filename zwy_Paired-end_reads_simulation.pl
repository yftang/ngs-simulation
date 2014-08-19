#!/usr/bin/perl -w

=pod
This Perl script simulated paire-end short reads,
written by Yifei Tang,
from Systems Biology Research Center, Soochow University,
We used part of thought by Juliane Dohm and Claudio Lottaz.
=cut

$Usage = <<EOF;
Users can simply run this script followed by 6 arguments:
	1. Genome sequence file in Fasta format from which you will generate reads.
	2. Number of read pairs to be generated.
	3. Length of each read.
	4. Distance between two paired reads.
	5. Assumed sequencing error rate (%).
	6. Base name of output files. (<basename>_PE.fa, rev_<basename>_PE.fa)
Eg. perl simulation_PE.pl Test.fa 50000 36 500 0.8 simread
This command will generates 50000 36-length short read pairs from "Test.fa"
with insert size of 500 and sequencing error rate 0.8%, output to simread_PE.fa and rev_simread_PE.fa.
EOF

#Check if a leagal command line is given.
die "$Usage\n" if (@ARGV != 6);

#Define variables.
my $in = $ARGV[0];
my $totalreads = $ARGV[1];
my $readlen = $ARGV[2];
my $gap = $ARGV[3];
my $errorrate = $ARGV[4];
my $basename = $ARGV[5];
my (@seq, @rev_seq, @chro, @chro_len, @chro_reads);
my $chro_index = 0, $total_len = 0, $errors = 0, $read_index = 0;

#Open input and output files.
open (IN, $in) || die "Genome sequence file not found!\n";
open (OUT1, ">$basename\_PE.fa");
open (OUT2, ">rev_$basename\_PE.fa");

#Read the genome sequence file line by line.
print "Reading genome sequence...\n";
while (<IN>) {
	chomp;
	if ($_ =~ /\>/) {
		if ($chro_index != 0) {
			my $rev_seq = reverse ($seq[$chro_index]);
			$rev_seq =~ tr/ATCGU/TAGCA/;
			$rev_seq =~ tr/atcgu/tagca/;
			$rev_seq[$chro_index] = $rev_seq;
			$chro_len[$chro_index] = length ($seq[$chro_index]);
			$total_len += $chro_len[$chro_index];
			print "Length of chromsome $chro_index... \"$chro[$chro_index]\" is $chro_len[$chro_index]\n";
		}
		my $chro = substr ($_, 1, length ($_));
		$chro[++$chro_index] = $chro;
		print "Processing chromsome $chro_index... \"$chro\"\n";
	}elsif ($_ =~ /[ATCGUatcguNn]/) {
		$seq[$chro_index] .= $_;
	}
	if (eof) {
		my $rev_seq = reverse ($seq[$chro_index]);
		$rev_seq =~ tr/ATCGU/TAGCA/;
		$rev_seq =~ tr/atcgu/tagca/;
		$rev_seq[$chro_index] = $rev_seq;
		$chro_len[$chro_index] = length ($seq[$chro_index]);
		$total_len += $chro_len[$chro_index];
		print "Length of chromsome $chro_index... \"$chro[$chro_index]\" is $chro_len[$chro_index]\n";
	}
}
print "Complete with reading genome sequence with the total length of $total_len.\n";

#Determine reads to be generated of each chromsome.
print "Determining reads to be generated of each chromsome...\n";
for (my $n = 1; $n <= $chro_index; $n++) {
	$chro_reads[$n] = int (($chro_len[$n] / $total_len) * $totalreads);
	print "Chromsome $chro[$n] will generate $chro_reads[$n] read pairs.\n";
}
print "Complete with determing reads to be generated of each chromsome.\n";

#Main part simulating short reads.
print "Generating short read pairs...\n";
for (my $n = 1; $n <= $chro_index; $n++) {
	for (my $m = 0; $m < $chro_reads[$n]; $m++) {
		$read_index++;
		my $start = int (rand ($chro_len[$n] - ($readlen * 2 + $gap)));
		my $revstart = $chro_len[$n] - ($readlen * 2 + $gap + $start);
		my ($ifrev, $revifrev);
		my ($read, $revread, $alt_read, $alt_revread);
		if (rand (1) >= 0.5) {
			$read = substr ($seq[$n], $start, $readlen);
			$ifrev = "+";
			$revread = substr ($rev_seq[$n], $revstart, $readlen);
			$revifrev = "-";
		}else {
			$read = substr ($rev_seq[$n], $start, $readlen);
			$ifrev = "-";
			$revread = substr ($seq[$n], $revstart, $readlen);
			$revifrev = "+";
		}
		$alt_read = $read, $alt_revread = $revread;
		$read = &seq_error ($read), $revread = &seq_error ($revread);
		print OUT1 ">read_$read_index\_$chro[$n]_$ifrev\_$start\n$read\n";
		print OUT2 ">reverse-read_$read_index\_$chro[$n]_$revifrev\_$start\n$revread\n";
		++$errors if ($alt_read ne $read);
		++$errors if ($alt_revread ne $revread);
	}
	print "Complete with generating short read pairs from $chro[$n].\n";
}
my $allreads = $read_index * 2;
print "Complete with generating all short reads.\n";
print "You get $errors errors out of $allreads of all generated reads.\a\n";

#Close files opened.
close IN;
close OUT1;
close OUT2;

#Subroutine processing sequencing error.
sub seq_error {
	my $tmpread = $_[0];
	my @base = split (//, $tmpread);
	$tmpread = "";
	for (my $n = 0; $n < @base; $n++) {
		if (rand (100) > $errorrate) {
			$tmpread .= $base[$n];
		}else {
			my %fourbase = ("A", 0, "T", 0, "C", 0, "G", 0);
			delete $fourbase{$base[$n]};
			my @fourbase = keys (%fourbase);
			my $i = rand (1);
			if ($i < 0.33) {
				$tmpread .= $fourbase[0];
			}elsif ($i <0.66) {
				$tmpread .= $fourbase[1];
			}else {$tmpread .= $fourbase[2];}
		}
	}
	return $tmpread;
}
