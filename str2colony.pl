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
getopts( 'A:c:d:e:hH:o:m:M:s:x:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $str, $out, $map, $missing, $error, $columns, $maxMissing, $minHet, $dominance, $dropout, $maf ) = &parsecom( \%opts );

# declare variables
my @strlines; # array to hold lines from structure file
my @maplines; # array to hold lines from species map file
my %hohoa; #hash of hashes of arrays to hold structure file contents
my %maphash; #key=individual, value=sex or status as offspring
my @inds; #array to hold list of individuals

# get number of loci in file
my $loci = `head -3 $str | tail -1 | sed 's/\\t/\\n/g' | wc -l`;
chomp($loci);
$loci = $loci-$columns;

# put files into array
&filetoarray( $str, \@strlines );
&filetoarray( $map, \@maplines );

# remove header lines
shift( @strlines );
shift( @strlines );

# put structure file into hash
&strToHash( \@strlines, \%hohoa, $missing );

&mapToHash( \@maplines, \%maphash );

my $badHash = &arrayCheck( \%hohoa, $maxMissing, $minHet, $maf );

#get list of all individuals in hash
foreach my $ind(sort keys %hohoa){
	push(@inds, $ind);
	#print $ind, "\n";
}

#&printInfo( $badHash, $loci, $dominance, $dropout, $error );
#&printFile( $out, $loci, $missing, $error, \@inds, \%hohoa, \%maphash, $badHash, $str );
my $remove = &blackList( \@inds, \%hohoa, $badHash, $loci, $str );

&printDat( $badHash, $loci, $dominance, $dropout, $error, \@inds, \%hohoa, \%maphash, $str, $remove );

#print Dumper(\%maphash);
#print Dumper(\%hohoa);

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
	print "str2colony.pl is a perl script developed by Steven Michael Mussmann\n\n";
	print "To report bugs send an email to mussmann\@email.uark.edu\n";
	print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
	print "Program Options:\n";
	print "\t\t[ -A | -c | -d | -e | -h | -H | -m | -M | -o | -s | -x ]\n\n";
	print "\t-A:\tUse this to specify the minimum allowed minor allele frequency.\n";
	print "\t\tDefault value is 0.05.\n\n";
	print "\t-c:\tUse this to specify the number of columns preceeding the first allele in the input structure file.\n";
	print "\t\tDefault value is 6 columns.\n\n";
	print "\t-d:\tUse this to specify whether loci are codominant (0) or dominant (1).\n";
	print "\t\tDefault value is 0.\n\n";
	print "\t-e\tUse this to specify the genotyping error rate for all loci in the file.\n";
	print "\t\tDefault value is 0.01.\n\n";
	print "\t-h:\tUse this flag to display this help message.\n";
	print "\t\tThe program will die after the help message is displayed.\n\n";
	print "\t-H:\tUse this to specify the minimum acceptable heterozygosity for a locus.\n";
	print "\t\tDefault value is 0.10.\n\n";
	print "\t-m:\tUse this flag to specify the name of the map file.\n";
	print "\t\tIf no name is provided, the program will fail to execute.\n\n";
	print "\t-M:\tUse this to specify the value of missing data in the input file.\n";
	print "\t\tDefault value is 0.\n\n";
	print "\t-o:\tUse this flag to specify the output file name.\n";
	print "\t\tIf no name is provided, the file extension \".csv\" will be appended to the input file name.\n\n";
	print "\t-s:\tUse this flag to specify the name of the structure file produced by Stacks.\n";
	print "\t\tIf no name is provided, the program will fail to execute.\n\n";
	print "\t-x:\tUse this to specify the maximum number of individuals allowed to have missing data at a locus.\n";
	print "\t\tDefault value is 60.\n\n";
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
	my( $params ) =  @_;
	my %opts = %$params;
  
	# set default values for command line arguments
	my $str = $opts{s} || die "No input file specified.\n\n"; #used to specify input file name.  This is the input snps file produced by pyRAD
	my $out = $opts{o} || "$str.csv"  ; #used to specify output file name.  If no name is provided, the file extension ".out" will be appended to the input file name.
	my $map = $opts{m} || die "No map file specified.\n\n"; #used to specify map file that will separate adult males, adult females, and juveniles
	my $missing = $opts{M} || "0";
	my $error = $opts{e} || "0.01";
	my $columns = $opts{c} || "6";
	my $maxMissing = $opts{x} || "60"; #max number of individuals allowed to have missing data for a column.
	my $minHet = $opts{H} || "0.10"; #minimum level of acceptable heterozygosity for a locus to be kept
	my $dominance = $opts{d} || "0"; #dominant (0) or codominant (1)
	my $dropout = $opts{D} || "0.000"; #dropout rate for locus
	my $maf = $opts{A} || "0.05"; #minimum cutoff for minor allele frequency

	return( $str, $out, $map, $missing, $error, $columns, $maxMissing, $minHet, $dominance, $dropout, $maf );

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

	my( $array, $hohoa, $missing ) = @_;

	my $counter = 0; #counter that will be returned with the number of loci in the file
	
	for(my $i = 0; $i < @$array; $i++){
		my @alleles = split(/\s+/, $$array[$i]);
		my $ind = shift(@alleles);
		my $pop = shift(@alleles);
		if($i%2==0){ #if even, put into allele 1
			foreach my $allele( @alleles ){
				if($allele != $missing){
					$allele+=100;
				}else{
					$allele=0;
				}
				push(@${$hohoa{$ind}{"allele1"}}, $allele);
			}
		}else{ #if odd, put into allele 2
			foreach my $allele( @alleles ){
				if($allele != $missing){
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

	my( $hohoa, $miss, $minHet, $maf ) = @_;
	
	my %allelecount;
	my %missingData;
	my %badHash;
	my %hetHash;
	my %mafHash;
	
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
	
	# missing data filter
	foreach my $locus( sort keys %missingData ){
		if($missingData{$locus}{0} > $miss){
			print $missingData{$locus}{0}, "\n";
			$badHash{$locus}++;
		}
	}
	
	# calculate heterozygosity for each locus
	foreach my $locus( sort keys %hetHash ){
		my $sum = $hetHash{$locus}{"het"} + $hetHash{$locus}{"hom"};
		if( $hetHash{$locus}{"het"} == 0){
			$hetHash{$locus}{"hz"}=0;
		}else{
			my $hz = $hetHash{$locus}{"het"} / $sum;
			$hetHash{$locus}{"hz"} = $hz;
		}
	}
	
	# Heterozygosity filter
	foreach my $locus( sort keys %hetHash ){
		if( sprintf( "%.2f", $minHet ) > sprintf( "%.2f", $hetHash{$locus}{"hz"} ) ){
			$badHash{$locus}++;
		}
	}
	
	my $totHet = 0.0;
	# calc heterozygosity
	foreach my $locus( sort keys %hetHash ){
		$totHet += $hetHash{$locus}{"hz"};
	}
	my $count = keys%hetHash;
	my $het = $totHet/$count;
	print "Heterozygosity for this file is ", $het, "\n";
	
	# singleton and MAF filter
	foreach my $locus( sort keys %allelecount ){
		my @tempArr;
		foreach my $allele( sort keys %{$allelecount{$locus}}){
			push(@tempArr, $allelecount{$locus}{$allele});
		}
		@tempArr = sort{ $a <=> $b } @tempArr;

		#filter out any locus that is a singleton
		if( $tempArr[0] == 1 ){
			$badHash{$locus}++;
		}
		
		#filter out any locus that has maf below minimum value
		my $total = 0;
		foreach my $num( @tempArr ){
			$total += $num;
		}
		my $freq = $tempArr[0] / $total;
		if( sprintf( "%.2f", $maf ) > sprintf( "%.2f", $freq ) ){
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
	#print $extra, "\n";
	
	$loci = $loci-$extra;
	
	open( OUT, '>', $file ) or die "Can't open $file: $!\n\n";
	open( MALEDATA, '>', "$str.data.male.csv" ) or die "Can't open output file: $!\n\n";
	open( FEMALEDATA, '>', "$str.data.female.csv" ) or die "Can't open output file: $!\n\n";
	open( OFFDATA, '>', "$str.data.offspring.csv" ) or die "Can't open output file: $!\n\n";

	foreach my $ind( @$indsArray ){
		if( $$maphash{$ind} eq "male" ){
			print OUT $ind;
			print MALEDATA $ind;
			for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
				if(!(exists $$badHash{$i})){
					print OUT ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
					print MALEDATA ",", $${$hohoa{$ind}{"allele1"}}[$i], ",", $${$hohoa{$ind}{"allele2"}}[$i];
				}
			}
			print OUT "\n";
			print MALEDATA "\n";
		}elsif( $$maphash{$ind} eq "female" ){
			print OUT $ind;
			print FEMALEDATA $ind;
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
			print OUT $ind;
			print OFFDATA $ind;
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
	close FEMALEDATA;
	close OFFDATA;
}

#####################################################################################################
# print marker info file

sub printInfo{

	my( $badhash, $loci, $dominance, $dropout, $err ) = @_;

	my $extra = keys %$badHash;
	#print $extra, "\n";
	
	$loci = $loci-$extra;
	
	open( OUT, '>', "MarkerInfo.txt") or die "Can't open MarkerInfo.txt for writing: $!\n\n";
	
	for( my $i=0; $i<$loci; $i++){
		print OUT "mk", $i+1, "\t";
	}
	print OUT "\n";
	
	for( my $i=0; $i<$loci; $i++){
		print OUT $dominance, "\t";
	}
	print OUT "\n";
	
	for( my $i=0; $i<$loci; $i++){
		print OUT $dropout, "\t";
	}
	print OUT "\n";
	
	for(my $i=0; $i<$loci; $i++){
		print OUT $err, "\t";
	}
	print "\n";
	
	close OUT;
	
	
	
}
#####################################################################################################
# subroutine to remove individuals with no data at kept loci

sub blackList{

	my( $indsArray, $hohoa, $badHash, $loci, $str ) = @_;
	
	my @fileArr = split(/\./, $str);
	my $file = $fileArr[0] . ".removed";
	
	my $extra = keys %$badHash;
	#print $extra, "\n";
	
	$loci = $loci-$extra;
	
	my %blacklist;
	my %remove;

	foreach my $ind( @$indsArray ){
		for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
			if(!(exists $$badHash{$i})){
				if($${$hohoa{$ind}{"allele1"}}[$i] == 0 && $${$hohoa{$ind}{"allele2"}}[$i] == 0){
					$blacklist{$ind}++;
				}
			}
		}
	}
	
	open( OUT, '>', $file ) or die "Can't open $file: $!\n\n";
	
	foreach my $ind( sort keys %blacklist ){
		if( $blacklist{$ind} == $loci ){
			print OUT "No data for $ind\n";
			$remove{$ind}++;
		}
	
	}
	
	close OUT;
	
	return( \%remove ); 
	
	
}

#####################################################################################################
# subroutine to print .dat file for colony2
sub printDat{

	my( $badHash, $loci, $dominance, $dropout, $error, $indsArray, $hohoa, $maphash, $str, $remove ) = @_;
	
	my @fileArr = split(/\./, $str);
	
	my $file = $fileArr[0] . ".dat";
	my $extra = keys %$badHash;
	#print $extra, "\n";
	
	$loci = $loci-$extra;
	
	my $males=0;
	my $females=0;
	my $offspring=0;
	
	foreach my $ind( @$indsArray ){
		if(!(exists $$remove{$ind})){
			if( $$maphash{$ind} eq "offspring" ){
				$offspring++;
			}elsif( $$maphash{$ind} eq "male" ){
				$males++;
			}elsif( $$maphash{$ind} eq "female" ){
				$females++;
			}
		}
	}
	
	open( OUT, '>', $file) or die "Can't open $file d-bag: $!\n\n";
	
	print OUT $fileArr[0], "\n";
	print OUT $fileArr[0], "\n";
	print OUT "$offspring        ! Number of offspring in the sample\n";
	print OUT "$loci       ! Number of loci\n";
	print OUT "1234      ! Seed for random number generator\n";
	print OUT "0         ! 0/1=Not updating/updating allele frequency\n";
	print OUT "2         ! 2/1=Dioecious/Monoecious species\n";
	print OUT "0         ! 0/1=Inbreeding absent/present\n";
	print OUT "0         ! 0/1=Diploid species/HaploDiploid species\n";
	print OUT "0  0      ! 0/1=Polygamy/Monogamy for males & females\n";
	print OUT "0         ! 0/1 = Clone inference = No/Yes\n";
	print OUT "1         ! 0/1=Scale full sibship=No/Yes\n";
	print OUT "1 1 1     ! 0/1/2/3/4=No/Weak/Medium/Strong sibship prior; 4=Optimal sibship prior for Ne\n";
	print OUT "0         ! 0/1=Unknown/Known population allele frequency\n";
	print OUT "1         ! Number of runs\n";
	print OUT "3         ! 1/2/3/4 = Short/Medium/Long/VeryLong run\n";
	print OUT "0         ! 0/1=Monitor method by Iterate#/Time in second\n";
	print OUT "10000     ! Monitor interval in Iterate# / in seconds\n";
	print OUT "0         ! 0/1=DOS/Windows version\n";
	print OUT "1         ! 0/1/2=Pair-Likelihood-Score(PLS)/Full-Likelihood(FL)/FL-PLS-combined(FPLS) method\n";
	print OUT "2         ! 0/1/2/3=Low/Medium/High/VeryHigh precision\n";
	print OUT "\n";
	
	#print marker information
	for( my $i=0; $i<$loci; $i++){
		print OUT "mk", $i+1, "\t";
	}
	print OUT "\n";
	
	for( my $i=0; $i<$loci; $i++){
		print OUT $dominance, "\t";
	}
	print OUT "\n";
	
	for( my $i=0; $i<$loci; $i++){
		print OUT $dropout, "\t";
	}
	print OUT "\n";
	
	for(my $i=0; $i<$loci; $i++){
		print OUT $error, "\t";
	}
	print OUT "\n";
	print OUT "\n";
	
	#print offspring data
	foreach my $ind( @$indsArray ){
		if(!(exists $$remove{$ind})){
			if( $$maphash{$ind} eq "offspring" ){
				print OUT $ind;
				for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
					if(!(exists $$badHash{$i})){
						print OUT " ", $${$hohoa{$ind}{"allele1"}}[$i], " ", $${$hohoa{$ind}{"allele2"}}[$i];
					}
				}
				print OUT "\n";
			}
		}
	}
	print OUT "\n";
	
	print OUT ".7 .7\n";
	print OUT "$males $females\n";
	print OUT "\n";
	
	foreach my $ind( @$indsArray ){
		if(!(exists $$remove{$ind})){
			if( $$maphash{$ind} eq "male" ){
				print OUT $ind;
				for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
					if(!(exists $$badHash{$i})){
						print OUT " ", $${$hohoa{$ind}{"allele1"}}[$i], " ", $${$hohoa{$ind}{"allele2"}}[$i];
					}
				}
				print OUT "\n";
			}
		}
	}
	print OUT "\n";
	
	foreach my $ind( @$indsArray ){
		if(!(exists $$remove{$ind})){
			if( $$maphash{$ind} eq "female" ){
				print OUT $ind;
				for( my $i=0; $i< @${$hohoa{$ind}{"allele1"}}; $i++ ){
					if(!(exists $$badHash{$i})){
						print OUT " ", $${$hohoa{$ind}{"allele1"}}[$i], " ", $${$hohoa{$ind}{"allele2"}}[$i];
					}
				}
				print OUT "\n";
			}
		}
	}
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	print OUT "0\n";
	print OUT "\n";
	
	
	close OUT;
}
#####################################################################################################
