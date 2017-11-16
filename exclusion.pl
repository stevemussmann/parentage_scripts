#! /usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'ho:s:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $str, $out ) = &parsecom( \%opts );

# declare variables
my @strlines; # array to hold lines from structure file
my %hohoa; #hash of hashes of arrays to hold structure file contents
my @inds; #array to hold list of individuals
my @combs;
my @pair1;
my @pair2;

# put files into array
&filetoarray( $str, \@strlines );

# make structure file into single line
for(my $i = 0; $i < @strlines; $i++){
	my @alleles = split(/\s+/, $strlines[$i]);
	my $ind = shift(@alleles);
	if($i%2==0){ #if even, put into allele 1
		foreach my $allele( @alleles ){
			push(@{$hohoa{$ind}{"allele1"}}, $allele);
		}
	}else{ #if odd, put into allele 2
		foreach my $allele( @alleles ){
			push(@{$hohoa{$ind}{"allele2"}}, $allele);
		}
	}
}

#get list of all individuals in hash
foreach my $ind(sort keys %hohoa){
	push(@inds, $ind);

}

# find all unique combinations of two individuals in the list
push(@combs, @$_) for &combinations( \@inds, 2 );

# split the output of the combinations subroutine into two arrays
for(my $i=0; $i<@combs; $i+=2){
	push(@pair1, $combs[$i]);
	push(@pair2, $combs[$i+1]);
	#print $combs[$i], ",", $combs[$i+1], "\n";

}

open(OUT, '>', $out) or die "Can't open $out: $!\n\n";

print OUT "Sample_1", "\t", "Sample_2", "\t", "Misssing_loci", "\t", "Incompatible_loci", "\t", "\%Missing", "\t", "\%Incompatible", "\n";

for(my $i=0; $i<@pair1; $i++){
	print "Working on pair $i of ", scalar(@pair1), ".\n";
	#print $pair1[$i], "\n";
	my $missing=0;
	my $incompatible=0;
	for(my $j=0; $j<@{$hohoa{$pair1[$i]}{"allele1"}}; $j++){
		#print $j, "\n";
		#print $hohoa{$pair1[$i]}{"allele1"}[$j], "\n";
		#check if either individual has missing data
		if( $hohoa{$pair1[$i]}{"allele1"}[$j]==-9 || $hohoa{$pair2[$i]}{"allele1"}[$j]==-9 || $hohoa{$pair1[$i]}{"allele2"}[$j]==-9 || $hohoa{$pair2[$i]}{"allele2"}[$j]==-9){
			$missing++;
			next;
		}else{
			if($hohoa{$pair1[$i]}{"allele1"}[$j]==$hohoa{$pair2[$i]}{"allele1"}[$j] || $hohoa{$pair1[$i]}{"allele1"}[$j]==$hohoa{$pair2[$i]}{"allele2"}[$j]){
				next;
			}elsif($hohoa{$pair1[$i]}{"allele2"}[$j]==$hohoa{$pair2[$i]}{"allele1"}[$j] || $hohoa{$pair1[$i]}{"allele2"}[$j]==$hohoa{$pair2[$i]}{"allele2"}[$j]){
				next;
			}else{
				$incompatible++;
			}
		}
	}
	
	my $loci=scalar(@{$hohoa{$pair1[$i]}{"allele2"}});
	
	my $percentmissing=($missing/$loci)*100;
	
	my $percentincompatible=($incompatible/($loci-$missing))*100;
	
	print OUT $pair1[$i], "\t", $pair2[$i], "\t", $missing, "\t", $incompatible, , "\t";
	printf OUT "%.4f", $percentmissing;
	print OUT "\t";
	printf OUT "%.4f", $percentincompatible;
	print OUT "\n";

}

close OUT;

#print Dumper(\%hohoa);

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nexclusion.pl is a perl script developed by Steven Michael Mussmann\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Program Options:\n";
  print "\t\t[ -h | -o | -s ]\n\n";
  print "\t-h:\tUse this flag to display this help message.\n";
  print "\t\tThe program will die after the help message is displayed.\n\n";
  print "\t-o:\tUse this flag to specify the output file name.\n";
  print "\t\tIf no name is provided, the file extension \".pairs\" will be appended to the input file name.\n\n";
  print "\t-s:\tUse this flag to specify the name of the structure file produced by pyRAD.\n\n";
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
  my( $params ) =  @_;
  my %opts = %$params;
  
  # set default values for command line arguments
  my $str = $opts{s} || die "No input file specified.\n\n"; #used to specify input file name.  This is the input snps file produced by pyRAD
  my $out = $opts{o} || "$str.pairs"  ; #used to specify output file name.  If no name is provided, the file extension ".out" will be appended to the input file name.


  return( $str, $out );

}

#####################################################################################################
# subroutine to put file into an array

sub filetoarray{

  my( $infile, $array ) = @_;

  
  # open the input file
  open( FILE, $infile ) or die "Can't open $infile: $!\n\n";

  # loop through input file, pushing lines onto array
  while( my $line = <FILE> ){
    chomp( $line );
    next if($line =~ /^\s*$/);
    #print $line, "\n";
    push( @$array, $line );
  }

  close FILE;

}

#####################################################################################################
# subroutine to get all unique pairs of individuals from an array

sub combinations{

	my( $array, $n ) = @_;
	
	if($n > scalar(@$array)){
		die "Number of individuals in array less than $n";
	}
	
	if($n <= 1){
		return map [$_], @$array;
	}
	
	my @comb;
	
	for(my $i=0; $i+$n <= @$array; $i++){
		my $val = $$array[$i];
		my @rest = @$array[$i+1..$#$array];
		push @comb, [$val, @$_] for combinations(\@rest, $n-1);
	}
	
	return @comb;

}
#####################################################################################################
