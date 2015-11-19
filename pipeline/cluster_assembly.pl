use warnings;
use strict;

my $species = $ARGV[0];
my $dir = '/home/ssinghal/encelia/' . $species . '/';
my $seq = '/home/ssinghal/encelia/final_assemblies/' . $species . '.fa.final';
my $anno_db = '/home/ssinghal/encelia/annotation/mguttatus.fa';
my $blastx = '/home/ssinghal/bin/ncbi-blast-2.2.29/bin/blastx';
my $cap3 = '/home/ssinghal/bin/CAP3/cap3';
my $blast_out1 = run_blast($seq);
my $clustered_seq = cluster($seq, $blast_out1);

sub run_blast {
	my ($seq) = @_;
	my $out = $seq . '.blast.out';
	
	unless(-f $out) {
		my $call = system("$blastx -query $seq -db $anno_db -out $out -outfmt 6 -max_target_seqs 1 -num_threads 8 -evalue 1e-10");
		}
	
	return($out);
	}

sub cluster {
	my ($seq, $blast_out1) = @_;
	my $out = $seq . '.clustered.fa';
	
	#unless (-f $out) {
		my %seq;
		my $c;
		open(IN, "<$seq");
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$c = $1;
				}
			else {
				$seq{$c} .= $line;
				}
			}
		close(IN);
		
		my %matches;
		
		open(IN, "<$blast_out1");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\s+/, $line);
			unless ($matches{$d[0]}) {
				my $match = $1 if $d[1] =~ m/^([^\|]+)/;
				$matches{$d[0]} = $match;
				}
			}
		close(IN);
		
		my %rev_matches;
		foreach my $c (keys %matches) {
			my $m = $matches{$c};
			push(@{$rev_matches{$m}}, $c);
			}
		
		my $tracker = 1;
		
		open(OUT, ">$out");
		foreach my $match (keys %rev_matches) {
			my $tmp = $dir . $match . '.tmp.fa';
			open(TMP, ">$tmp");
			foreach my $c (@{$rev_matches{$match}}) {
				print TMP ">", $c, "\n", $seq{$c}, "\n";
				delete $seq{$c};
				}
			close(TMP);
			
			my $cap_call = system($cap3 . " " . $tmp . " -v 1 -u 1 -p 99 -o 16 -z 1 -e 11 -s 251 -h 99");
				
			my $assembled = $tmp . ".cap.contigs";
			my $singlets = $tmp . ".cap.singlets";
			
			open(SIN, "<$singlets");
			while(<SIN>) {
				chomp(my $line = $_);
				if ($line =~ m/>/) {
					print OUT ">contig", $tracker, "\n";
					$tracker++;
					}
				else {
					print OUT $line, "\n";
					}
				}
			close(SIN);
						
			open(CON, "<$assembled");
			while(<CON>) {
				chomp(my $line = $_);
				if ($line =~ m/>/) {
					print OUT ">contig", $tracker, "\n";
					$tracker++;
					}
				else {
					print OUT $line, "\n";
					}
				}
			close(CON);
			
			my $call = system("rm $tmp*");
			}
			
		foreach my $s (keys %seq) {
			print OUT ">contig", $tracker, "\n", $seq{$s}, "\n";
			$tracker++;
			}
		close(OUT);
	#	}

	return($out);
	}
