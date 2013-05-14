###Jeff Ross-Ibarra
###Nov. 2012
###Code samples sequences randomly from a fastq file.  
###Run as "samperl.pl <FASTQ FILE> <NUMBER OF READS>"

#!/usr/bin/perl
use POSIX;

#get file name, sample size, die and show usage
die "usage: samperl.pl <FASTQ FILE> <NUMBER OF READS>\n\n" unless $#ARGV==1;
my $file=$ARGV[0];
my $samplesize=$ARGV[1];

#length of file
my $long=`wc -l $file | awk '{print \$1 }'`;
#number reads
$long=$long/4;

#get random sample of size $samplesize
my @lines=();
push(@lines,floor(rand()*$long)*4) for 0..$samplesize-1;
@lines = sort {$a <=>$b} @lines;

#print reads corresponding to random sample
open FILE, "<$file";

my $counter=0;
my $next=0;
while(<FILE>){
	$counter++;
	if($counter > $lines[$next]){
		print $_;
		$next++ if $counter==$lines[$next]+4;	
	}
	last if $next>$#lines;
}
