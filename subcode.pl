use strict;
use warnings;
use Getopt::Std;

# get options and input files
my %options=();
getopts("i:n:j:",\%options) or help();
my $file=$options{i}; 
my $n=$options{n};
help() unless $options{i};
help() unless $options{n};

#hashes, yeah
my %barcodes=();
my %secondbarcodes=();

#get barcodes from file of extra desired codes
#file format is barcode\tnumber
if($options{j}){
	my $twofile=$options{j};
	open FILE, "<$twofile" or die "no 2nd barcode file, your silly turnip!\n\n";
	while(<FILE>){
		chomp $_;
		my @line=split(/\t/,$_);
		$secondbarcodes{$line[0]}=$line[1];
	}
	close FILE;
}

#get barcodes from file of 96 possible barcodes
#file format is barcode\tnumber
open FILE, "<$file" or die "no barcode file, your silly turnip!\n\n";
while(<FILE>){
	chomp $_;
	my @line=split(/\t/,$_);
	$barcodes{$line[0]}=$line[1];
}
close FILE;

#gets barcode sequences and takes random sample of first $n
my @seqs=keys(%barcodes);
my @twoseqs=keys(%secondbarcodes);
my $req=@twoseqs;
my $diverse=0;
my @codes=();
my $it=1;

#loops through random sets of barcodes until finds ones that match 1) nucleotide diversity of all 4 nukes at each bp and 2) distance of ≥2 between all barcodes
until($diverse){
#	print STDERR "iteration $it\n" if $it % 1000==0 ;
	die "This has gone on too long. At this many iterations it is unlikely you will find a set of barcodes.\n\n" if $it > 1.10*fact(96)/(fact($n)*fact(96-$n));
	@codes=(@{shuffle($n-$req,\@seqs)},@twoseqs);
	$it++;		
	$diverse=1;

	#check to make sure 4 nukes represented for each bp of code
	my @nukes=@{count(\@codes)};
	foreach(@nukes){ $diverse=0 if $_<4; }
	
	#check hamming distance ≥2 for each
	foreach(@codes){
		$diverse=0 if hammy(\@codes)==0; 
	}	
}

print "@codes\n\n";


####SUBS

#shuffle barcodes. code lifted from:
#http://stackoverflow.com/questions/8963228/how-can-i-take-n-elements-at-random-from-a-perl-array
sub shuffle{
	my $picks_left = $_[0]; 
	my @rseqs=@{$_[1]};
	my $num_left = @rseqs;
	my @picks;
	my $idx = 0;
	# when we have all our picks, stop
	while($picks_left > 0 ) {  
		# random number from 0..$num_left-1
		my $rand = int(rand($num_left));
		# pick successful
		if( $rand < $picks_left ) {
			push @picks, $rseqs[$idx];
        		$picks_left--;
    		}
		$num_left--;
		$idx++;
	}
	return(\@picks);
}

#code for binomial coefficient from http://www.perlmonks.org/?node_id=777339
#sub bc($n,$k){[*]map ->$x, $y{$x/$y},(($n-$k+1..$n)Z(1..$k));}
#JRI's ludgy hack at it
sub fact{
	my $nn=shift();
	my $f = 1; $f *= $_ for 1..$nn;
	return($f);
}

#count hamming distance
sub hammy{
	my @hcodes=@{$_[0]};
	for my $start(0..$#hcodes-1){
		for my $comp ($start+1..$#hcodes){	
			my $diff=0;
			my $seq1=$hcodes[$start];
			my $seq2=$hcodes[$comp];
			for my $bp (0..4){
				$diff++ if substr($seq1,$bp,1) ne substr($seq2,$bp,1);				
			}
			return(0) unless $diff>=2;			
		}

	}
	return(1);
}

#counts of each barcode
sub count{
	my %agct=();
	my @nuke_clear=((0) x 5);
	my @count_codes=@{$_[0]};

	foreach my $c (@count_codes){
		#make uppercase
		$c=~tr/a-z/A-Z/;

		my @bp=split(//,$c);
		my $i=0;	
		foreach my $b (@bp){
			$i++;
			next if $agct{$i}->{$b};
			$nuke_clear[$i-1]++;
			$agct{$i}->{$b}=1;
		}
	}
	return(\@nuke_clear);
}

sub help{

	die "\nBarcode subthingy version 1.0

Usage: perl subcode.pl -i <barcode filename> -n <number of barcodes> -j <optional file of barcodes you want included>

Barcode file format should be <barcode>\t<number> on each line";

}
