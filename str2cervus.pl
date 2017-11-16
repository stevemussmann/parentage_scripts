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
getopts( 'c:e:hH:o:m:M:s:x:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $str, $out, $map, $missing, $error, $columns, $maxMissing, $minHet ) = &parsecom( \%opts );

# declare variables
my @strlines; # array to hold lines from structure file
my @maplines; # array to hold lines from species map file
my %hohoa; #hash of hashes of arrays to hold structure file contents
my %maphash; #key=individual, value=sex or status as offspring
my @inds; #array to hold list of individuals

# get number of loci in file
my $loci = `head -n1 $str | sed 's/\\t/\\n/g' | wc -l`;
chomp($loci);
$loci = $loci-$columns;

# put files into array
&filetoarray( $str, \@strlines );
&filetoarray( $map, \@maplines );

# put structure file into hash
&strToHash( \@strlines, \%hohoa );

&mapToHash( \@maplines, \%maphash );

my $badHash = &arrayCheck( \%hohoa, $maxMissing, $minHet );

#get list of all individuals in hash
foreach my $ind(sort keys %hohoa){
	push(@inds, $ind);
	#print $ind, "\n";
}


&printFile( $out, $loci, $missing, $error, \@inds, \%hohoa, \%maphash, $badHash, $str );

#print Dumper(\%maphash);
#print Dumper(\%hohoa);

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
	print "str2cervus.pl is a perl script developed by Steven Michael Mussmann\n\n";
	print "To report bugs send an email to mussmann\@email.uark.edu\n";
	print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
	print "Program Options:\n";
	print "\t\t[ -c | -e | -h | -H | -m | -M | -o | -s | -x ]\n\n";
	print "\t-c:\tSpecify number of columns preceding the first locus.  Default=6\n\n";
	print "\t-e:\tSpecify error rate. Default=0.005\n\n";
	print "\t-h:\tUse this flag to display this help message.\n";
	print "\t\tThe program will die after the help message is displayed.\n\n";
	print "\t-H:\tSpecify the minimum level of acceptable heterozygosity for a locus to be kept. Default=0.20\n\n";
	print "\t-m:\tUse this flag to specify the name of the map file.\n";
	print "\t\tIf no name is provided, the program will fail to execute.\n\n";
	print "\t-M:\tSpecify the missing data value that is used in the structure file. Default=-9\n\n";
	print "\t-o:\tUse this flag to specify the output file name.\n";
	print "\t\tIf no name is provided, the file extension \".csv\" will be appended to the input file name.\n\n";
	print "\t-s:\tUse this flag to specify the name of the structure file produced by pyRAD.\n\n";
	print "\t-x:\tSpecify the maximum number of individuals allowed to have missing data for a column to be retained. Default=10\n\n";
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
	my( $params ) =  @_;
	my %opts = %$params;
  
	# set default values for command line arguments
	my $str = $opts{s} || die "No input file specified.\n\n"; #used to specify input file name.  This is the input snps file produced by pyRAD
	my $out = $opts{o} || "$str.csv"  ; #used to specify output file name.  If no name is provided, the file extension ".csv" will be appended to the input file name.
	my $map = $opts{m} || die "No map file specified.\n\n"; #used to specify map file that will separate adult males, adult females, and juveniles
	my $missing = $opts{M} || "-9";
	my $error = $opts{e} || "0.005";
	my $columns = $opts{c} || "6";
	my $maxMissing = $opts{x} || "10"; #max number of individuals allowed to have missing data for a column.
	my $minHet = $opts{H} || "0.20"; #minimum level of acceptable heterozygosity for a locus to be kept

	return( $str, $out, $map, $missing, $error, $columns, $maxMissing, $minHet );

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
# put structure file into hash

sub strToHash{

	my( $array, $hohoa ) = @_;

	my $counter = 0; #counter that will be returned with the number of loci in the file
	
	for(my $i = 0; $i < @$array; $i++){
		my @alleles = split(/\s+/, $$array[$i]);
		my $ind = shift(@alleles);
		if($i%2==0){ #if even, put into allele 1
			foreach my $allele( @alleles ){
				if($allele != -9){
					$allele+=100;
				}else{
					$allele=0;
				}
				push(@${$hohoa{$ind}{"allele1"}}, $allele);
			}
		}else{ #if odd, put into allele 2
			foreach my $allele( @alleles ){
				if($allele != -9){
					$allele+=100;
				}else{
					$allele=0;
				}
				push(@${$hohoa{$ind}{"allele2"}}, $allele);
			}
		}
	}
}

#####################################################################################################
# put map file into hash

sub mapToHash{

	my( $array, $hash ) = @_;
	
	foreach my $item( @$array ){
		my @pair = split(/\s+/, $item);
		my $ind = $pair[0];
		my $status = $pair[1];
		$$hash{$ind} = $status;
	}
}

#####################################################################################################
# checks for loci with more than one allele

sub arrayCheck{

	my( $hohoa, $miss, $minHet ) = @_;
	
	my %allelecount;
	my %missingData;
	my %badHash;
	my %hetHash;
	
	foreach my $ind( sort keys %$hohoa ){
		for( my $i=0; $i<@${$hohoa{$ind}{"allele1"}}; $i++){
			$hetHash{$i}{"het"}=0;
			$hetHash{$i}{"hom"}=0;
		}
	}
	
	foreach my $ind( sort keys %$hohoa ){
		for( my $i=0; $i<@${$hohoa{$ind}{"allele1"}}; $i++){
			if($${$hohoa{$ind}{"allele1"}}[$i] != 0){
				$allelecount{$i}{$${$hohoa{$ind}{"allele1"}}[$i]}++;
			}else{
				$missingData{$i}{$${$hohoa{$ind}{"allele1"}}[$i]}++;
			}
			if($${$hohoa{$ind}{"allele2"}}[$i] != 0){
				$allelecount{$i}{$${$hohoa{$ind}{"allele2"}}[$i]}++;
			}
			if( $${$hohoa{$ind}{"allele1"}}[$i] != 0 && $${$hohoa{$ind}{"allele2"}}[$i] != 0 && $${$hohoa{$ind}{"allele1"}}[$i] != $${$hohoa{$ind}{"allele2"}}[$i] ){
				$hetHash{$i}{"het"}++;
			}elsif( $${$hohoa{$ind}{"allele1"}}[$i] != 0 && $${$hohoa{$ind}{"allele2"}}[$i] != 0 && $${$hohoa{$ind}{"allele1"}}[$i] == $${$hohoa{$ind}{"allele2"}}[$i] ){
				$hetHash{$i}{"hom"}++;
			}
		}
	}

# 	foreach my $locus( sort keys %allelecount ){
# 		my $n = keys %{$allelecount{$locus}};
# 		if($n != 2){
# 			print $n, "\n";
# 			$badHash{$locus}++;
# 		}
# 		#print $n, "\n";
# 	}
	
	foreach my $locus( sort keys %missingData ){
		if($missingData{$locus}{0} > $miss){
			print $missingData{$locus}{0}, "\n";
			$badHash{$locus}++;
		}
	}
	
	foreach my $locus( sort keys %hetHash ){
		my $sum = $hetHash{$locus}{"het"} + $hetHash{$locus}{"hom"};
		if( $hetHash{$locus}{"het"} == 0){
			$hetHash{$locus}{"hz"}=0;
		}else{
			my $hz = $hetHash{$locus}{"het"} / $sum;
			$hetHash{$locus}{"hz"} = $hz;
		}
	}
	
	foreach my $locus( sort keys %hetHash ){
		if( sprintf( "%.2f", $minHet ) > sprintf( "%.2f", $hetHash{$locus}{"hz"} ) ){
			$badHash{$locus}++;
		}
	}
	
	
	#print Dumper(\%allelecount);
	#print Dumper( \%badHash );
	#print Dumper( \%hetHash );
	
	return( \%badHash );
}



#####################################################################################################
# print output file

sub printFile{

	my( $file, $loci, $missing, $error, $indsArray, $hohoa, $maphash, $badHash, $str ) = @_;

	my $extra = keys %$badHash;
	print $extra, "\n";
	
	$loci = $loci-$extra;
	
	open( OUT, '>', $file ) or die "Can't open $file: $!\n\n";
	open( MALELIST, '>', "$str.list.male.csv" ) or die "Can't open output file: $!\n\n";
	open( FEMALELIST, '>', "$str.list.female.csv" ) or die "Can't open output file: $!\n\n";
	open( OFFLIST, '>', "$str.list.offspring.csv" ) or die "Can't open output file: $!\n\n";
	open( MALEDATA, '>', "$str.data.male.csv" ) or die "Can't open output file: $!\n\n";
	open( FEMALEDATA, '>', "$str.data.female.csv" ) or die "Can't open output file: $!\n\n";
	open( OFFDATA, '>', "$str.data.offspring.csv" ) or die "Can't open output file: $!\n\n";
	
	
	print OUT "Individual,Sex";
	print MALEDATA "Individual,Sex";
	print FEMALEDATA "Individual,Sex";
	print OFFDATA "Individual,Sex";
	for(my $i=0; $i<$loci; $i++){
		print OUT ",Locus", $i+1, "a", ",Locus", $i+1, "b";
		print MALEDATA ",Locus", $i+1, "a", ",Locus", $i+1, "b";
		print FEMALEDATA ",Locus", $i+1, "a", ",Locus", $i+1, "b";
		print OFFDATA ",Locus", $i+1, "a", ",Locus", $i+1, "b";
	}
	print OUT "\n";
	print MALEDATA "\n";
	print FEMALEDATA "\n";
	print OFFDATA "\n";

	foreach my $ind( @$indsArray ){
		if( $$maphash{$ind} eq "male" ){
			print OUT $ind, ",M";
			print MALEDATA $ind, ",M";
			print MALELIST $ind, "\n";
			for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
				if(!(exists $$badHash{$i})){
					print OUT ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
					print MALEDATA ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
				}
			}
			print OUT "\n";
			print MALEDATA "\n";
		}elsif( $$maphash{$ind} eq "female" ){
			print OUT $ind, ",F";
			print FEMALEDATA $ind, ",F";
			print FEMALELIST $ind, "\n";
			for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
				if(!(exists $$badHash{$i})){
					print OUT ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
					print FEMALEDATA ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
				}
			}
			print OUT "\n";
			print FEMALEDATA "\n";
		}
	}
	foreach my $ind( @$indsArray ){
		if( $$maphash{$ind} eq "offspring" ){
			print OUT $ind, ",O";
			print OFFDATA $ind, ",O";
			print OFFLIST $ind, "\n";
			for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
				if(!(exists $$badHash{$i})){
					print OUT ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
					print OFFDATA ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
				}
			}
			print OUT "\n";
			print OFFDATA "\n";
		}
	}
	
	close OUT;
	close MALEDATA;
	close MALELIST;
	close FEMALEDATA;
	close FEMALELIST;
	close OFFDATA;
	close OFFLIST;
}

#####################################################################################################

# Individual ID,Sex,FCB193a,FCB193b,FCB304a,FCB304b,JP15a,JP15b,...
# 6MEL,2,101,111,131,145,113,123,...
# 7AVO,2,0,0,133,133,113,113,...
# 7CFU,2,99,99,129,133,111,119,...
# 7ANK,2,111,119,133,145,113,113,...
# 7MUS,2,101,111,131,131,0,0,...
# 8EGR,2,99,111,131,133,123,123,...This comma-separated (CSV) file contains two initial columns Individual ID and Sex, followed by six columns giving the genotypes at three microsatellite loci, FCB193, FCB304 and JP15. The demo file in fact contains data for a total of nine microsatellite loci.

# The Sex column indicates if the individual is male or female (in this case they are all males). The Sex column is optional. It is not used by Cervus during parentage analysis, although it is useful for identity analysis. You can include as many additional columns as you like.

# In this example alleles are coded numerically by microsatellite allele length and missing data is coded as 0 (zero); elsewhere in the actual file some other missing data is coded as blank.
