#!/usr/bin/env perl
use strict;
use warnings;

die usage() if @ARGV == 0;
my ($gff,$genome,$output) = @ARGV;

my %hash_LTR_sequence;
open NEW,"$gff" or die;
while(<NEW>){
	chomp;
	next if (/^#/);
	my @array = split /\s+/;
	if($array[2] eq "long_terminal_repeat"){
		my $start = $array[3] - 1;
		my $end = $array[4];

		system("samtools faidx $genome $array[0]:$start-$end > ./temp.fasta");

		my $seq_temp;
		open TEMP,"./temp.fasta" or die;
		while(<TEMP>){
			chomp;
			if(/^>/){
				next;
			}
			else{
				$seq_temp .= $_;
			}
		}
		close TEMP;

		system("rm ./temp.fasta");

		$array[-1] =~ /Parent=(.*);Method/;
		push @{$hash_LTR_sequence{$1}},$seq_temp;
	}
}
close NEW;

open OUT,">$output" or die;
print OUT "Element\tDistance\tInsertion_time\n";

foreach my $ltr(sort keys %hash_LTR_sequence){
	open NEW,">./seq.fasta" or die;
	print NEW ">seq1\n$hash_LTR_sequence{$ltr}[0]\n>seq2\n$hash_LTR_sequence{$ltr}[1]\n";
	close NEW;


	system("muscle -in ./seq.fasta -out ./seq.muscle -clwstrict");
	system("distmat -sequence ./seq.muscle -nucmethod 2 -outfile ./seq.dist");

	my $dist;
	open NEW1,"./seq.dist" or die;
	while(<NEW1>){
		chomp;
		next if (/Distance/ or /--/ or /Using/ or /^Gap/ or /^\s+$/ or /^$/);
		my @aa = split /\s+/;
		if($aa[-1] == 1){
			$dist = $aa[-3]/100;
			my $time = $dist/(2*1.3*10**-8);
			print OUT "$ltr\t$dist\t$time\n";
		}
	}
	close NEW1;

	system("rm ./seq.fasta ./seq.muscle ./seq.dist");

}

close OUT;



sub usage{
	my $die =<<DIE;
	usage : perl *.pl Me34V_genome_ctg.cor.pseudomolecule.fa.mod.EDTA.intact.gff Me34V_genome_ctg.cor.pseudomolecule.2bit output_file
DIE
}
