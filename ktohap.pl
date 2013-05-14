use strict;
use warnings;

#get SNP information on positions etc.
open SNPDATA, "<snpdata.txt" or die "Need data on SNPs in snpdata.txt\n\n";
my %snpdata; my $disambiguity=0;
while(<SNPDATA>){
	chomp;
	#name\tchr\tpos\ttype -- last three are in "@stuff"
	my ($snp, @stuff)=split(/\t/,$_);
	$snpdata{$snp}=\@stuff;
}
close SNPDATA;


#get input and output file names
my ( $infile, $outfile );
my $format="hapmap"; my $institution="UCD"; my $panel="NA";
for my $arg (0..$#ARGV){
	if( $ARGV[$arg] eq "-i" ){ $infile=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "-o" ){ $outfile=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "-f" ){ $format=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "-d" ){ $disambiguity=1; }
	if( $ARGV[$arg] eq "--panel" ){ $panel=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "--inst" ){ $institution=$ARGV[$arg+1]; }
	if( $ARGV[$arg] eq "-h" || $ARGV[$arg] eq "--help" ){ die help(); }
}
die help() unless $infile;
die help() unless $outfile;
$outfile=$outfile . ".hmp.txt";


#reads in SNP data from Genome Studio file, changes hets to ambiguity codes, ditches 100% missing data
open IN, "<$infile";
my $bit=0; my @inds; my %genos; my %alleles;
while(<IN>){ 
	if( $_=~m/Data/ ){ $bit++; next; }
	next unless $bit;
	if( $bit == 1 ){ $bit++; chomp; @inds=split(/\t/,$_); shift(@inds); }
	else{
		chomp; $_=~s/\s+$//g;
		my $snp; my @data;
		($snp, @data)=split(/\t/, $_);
		$genos{$snp}=\@data;
		foreach(@data){ 
			if( $disambiguity ){
				$_=~s/--/NN/;
				my $first = substr($_,0,1);
				my $second = substr($_,1,1);
				next if $first=~/N/;
				$alleles{$snp}=$first unless $alleles{$snp}; 
				$alleles{$snp}.="/$first" unless $alleles{$snp}=~/$first/;
				$alleles{$snp}.="/$second" unless $alleles{$snp}=~/$second/;		
			}
			else{
				$_=ambiguity($_);
				next if $_=~/N/;
				$alleles{$snp}=$_ unless $alleles{$snp}; 
				$alleles{$snp}.="/$_" unless $alleles{$snp}=~/$_/;
			}
		}		
		$genos{$snp}=undef unless $alleles{$snp};
	}
}
close IN;


#print to hapmap format
if( $format eq "hapmap" ){
	open HAPOUT, ">$outfile";
	$"="\t";
	print HAPOUT "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\t@inds\n";
	#change array separator for printing;
	for my $c (0..10){ 
		my %temp;
		foreach(keys(%snpdata)){
			$temp{$_}=$snpdata{$_} if $snpdata{$_}[0]==$c;
		}
		my @sorted = sort { $temp{$a}[1] <=> $temp{$b}[1] } keys %temp; 
		foreach(@sorted){
			print HAPOUT "$_\t", $alleles{$_}, "\t", $snpdata{$_}[0], "\t", $snpdata{$_}[1], "\t+\t", $snpdata{$_}[2], "\t$institution\tNA\tNA\t$panel\t00\t@{$genos{$_}}\n" if $genos{$_}; 
		}
	}
	close HAPOUT;
	$"=" "; # change it back
}

#prints out help info
sub help
{
print "
 ########################################################################
#  Copyright (C) 2011 Jeffrey Ross-Ibarra <rossibarra\@ucdavis.edu>      #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License <http://www.gnu.org/licenses/>             #
#  for more details.                                                     #
#                                                                        #
 ########################################################################

	ktohap version 1.0 can be run as follows:
	
	ktohap.pl -i <INFILE> -o <OUTFILE> <OPTIONS>
	
	Where <INFILE> is a file in the default Genome Studio output format
	and <OUTFILE> is the base name of the hapmap-formatted output file.

	Optional flags include:
	
	-f	<FORMAT> Currently only accepts teh default \"hapmap\" but other formats may be available soon
	-d Makes ktohap report both alleles instead of treating heterozygotes as ambiguity codes.
	--inst <NAME> Defines institution column.  Default is UCD
	--panel <NAME> Defines panel column.  Default is \"NA\";
	\n\n";
}

#turns hets into ambiguity codes.  Woo!
sub ambiguity
{
	my $allele=$_[0];
	#ambiguize them damn codes
	$allele=~s/GT|TG/K/; $allele=~s/--/N/; $allele=~s/AC|CA/M/; $allele=~s/AG|GA/R/; $allele=~s/CT|TC/Y/; $allele=~s/CG|GC/Y/; $allele=~s/AT|TA/W/; 
	return(substr($allele,0,1));
}