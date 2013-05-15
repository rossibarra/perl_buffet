# geneious_to_fasta.pl:

A simple perl script that converts files in the .genious format to fasta format. 

Runs with the following command:  

genious_to_fasta.pl -i <filename>


# SFS.pl
Usage:	SFS.pl -F <number of files> <filenames>

the <filenames> argument accepts the wildcard * character for multiple files

Optional command line arguments:

-q removes header (only if -f is not also chosen)
-n new sample size
-f prints out frequency counts of derived mutations at every site
--condition conditions on observing a site as polymorphic in the new sample size
--unrounded returns unrounded frequencies
--folded calculates folded SFS
--version prints version and license information


fastU.pl

Required command line parameters:

               -I 'filename' : input name of fasta file(s) (accepts * wildcard).

Optional command line parameters:

               -V : for verbose output
               -N : turn end gaps into N's
               -M : turn all characters except A,C,G,T,_ into N's
               -O : outputs to file 'filename'.out instead of stdout
               -C : changelog since version 1.0
               -S : splits alignment into separate files for groups of sequence ID's
               -G : replaces Genbank ID's with species name (requires Bioperl installation)


ktohap.pl
The default format looks something like:

[Header]
GSGT Version	1.1.9
Processing Date	5/3/2011 10:10 AM
Content		MaizeSNP50_A.bpm
Num SNPs	55126
Total SNPs	56110
Num Samples	94
Total Samples	94
[Data]
	RIMMA0662.1	RIMMA0665.1	RIMMA0394.2	RIMMA0404.1
abph1.15	GG	GG	AA	AA	
abph1.22	AA	AA	AA	AA	

ktohap version 1.0 can be run as follows:
	
	ktohap.pl -i <INFILE> -o <OUTFILE> <OPTIONS>
	
	Where <INFILE> is a file in the default Genome Studio output format
	and <OUTFILE> is the base name of the hapmap-formatted output file.

	Also requires a file called “snpdata.txt” that has information on SNP locations and types, in the following format:

	PZE0000200501   0       4894262 hapmap
	PZE0000344524   0       4921654 hapmap
	PZE0008211605   0       2572930 hapmap

	Where first column is SNP name, then chromosome, then position on chromosome, then SNP type.

	Optional flags include:
	
	-f	<FORMAT> Currently only accepts teh default \"hapmap\" but other formats may be available soon
	-d Makes ktohap report both alleles instead of treating heterozygotes as ambiguity codes.
	--inst <NAME> Defines institution column.  Default is UCD
	--panel <NAME> Defines panel column.  Default is \"NA\";


HapToTab.pl

	HapToTab version 2.1 can be run as follows:
	
	HapToTab.pl <INFILE> <--snp or --bp> <INT> <OPTIONS>
	
	Where <INFILE> is a file in hapmap format.  
	HapToTab needs to know the window size either in SNPs or bp (e.g. --snp 100 or --bp 10000).  
	A negative number for SNPs will run HapToTab using only one window for the whole file.  
	The following options are available:
	
	--g <NAME> name the group you are running haptotab on.  default is \"all\"
	--maxdip format the output files for Hudson's maxdip program
	--maxhap format the output files for Hudson's maxhap program
	--nohets remove all heterozygotes (code as missing data)
	--out file has an outgroup (assumes the last column in the hapmap file)
	--freq <INT> filters out all SNPs with a frequency less than INT in the data\n\n";


samperl.pl
Code samples sequences randomly from a fastq file.  
Run as "samperl.pl <FASTQ FILE> <NUMBER OF READS>"

