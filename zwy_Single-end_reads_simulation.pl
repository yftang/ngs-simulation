#!/usr/bin/perl -w

=pod
This Perl script simulated single-end short reads,
written by Yifei Tang,
from Systems Biology Research Center, Soochow University,
We used part of thought by Juliane Dohm and Claudio Lottaz.
=cut

$Usage = <<EOF;
Users can simply run this script followed by 5 arguments:
	1. Genome sequence file in Fasta format from which you will generate reads.
	2. Number of reads to be generated.
	3. Length of each read.
	4. Assumed sequencing error rate (%).
	5. Base name of output files. (<basename>_SE.fa)
Eg. perl simulation_SE.pl Test.fa 50000 36 0.8 simread
This command will generates 50000 36-length short reads from "Test.fa"
with sequencing error rate 0.8%, output to simread_SE.fa.
EOF

#Check if a leagal command line is given.
die "$Usage\n" if (@ARGV != 5);

#Define variables.
my $in = $ARGV[0];
my $totalreads = $ARGV[1];
my $readlen = $ARGV[2];
my $errorrate = $ARGV[3];
my $basename = $ARGV[4];
my (@seq, @rev_seq, @chro, @chro_len, @chro_reads);
my $chro_index = 0, $total_len = 0, $errors = 0, $read_index = 0;

#Open input and output files.
open (IN, $in) || die "Genome sequence file not found!\n";
open (OUT, ">$basename\_SE.fa");

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
	print "Chromsome $chro[$n] will generate $chro_reads[$n] reads.\n";
}
print "Complete with determing reads to be generated of each chromsome.\n";

#Main part simulating short reads.
print "Generating short read pairs...\n";
for (my $n = 1; $n <= $chro_index; $n++) {
	for (my $m = 0; $m < $chro_reads[$n]; $m++) {
		$read_index++;
		my $start = int (rand ($chro_len[$n] - $readlen));
		my $ifrev;
		my ($read, $alt_read);
		if (rand (1) >= 0.5) {
			$read = substr ($seq[$n], $start, $readlen);
			$ifrev = "+";
		}else {
			$read = substr ($rev_seq[$n], $start, $readlen);
			$ifrev = "-";
		}
		$alt_read = $read;
		$read = &seq_error ($read);
		print OUT ">read_$read_index\_$chro[$n]_$ifrev\_$start\n$read\n";
		++$errors if ($alt_read ne $read);
	}
	print "Complete with generating short reads from $chro[$n].\n";
}
my $allreads = $read_index;
print "Complete with generating all short reads.\n";
print "You get $errors errors out of $allreads of all generated reads.\a\n";

#Close files opened.
close IN;
close OUT;

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
