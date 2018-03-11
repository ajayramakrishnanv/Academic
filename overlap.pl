#!/usr/bin/perl -w
use strict;
use warnings;
##Defining functions
sub min {
 	my $a=$_[0];
 	my $b=$_[1];
 	if ($a > $b) {
 		return $b;
 	}
 	else{
 		return $a;
 	}
}	

sub max {
	my $c=$_[0];
	my $d=$_[1];
	if ($c > $d) {
		return $c;
	}
	else{
		return $d;
	}
}	 

##Defining initial variables
my $file1;
my $file2;
my $perc;
my $out;
my $join;
my $fh3;
my $cut;

#Getopts emulation
for (my $var = 0; $var < scalar @ARGV; $var++) {
	if ($ARGV[$var] =~ /i1/) {
		$file1 = $ARGV[$var+1];
	}
	elsif ($ARGV[$var] =~ /i2/){
		$file2=$ARGV[$var+1];
	}
	elsif($ARGV[$var] =~ /\-m/){
		$perc=$ARGV[$var+1]/100;
	}
	elsif($ARGV[$var] =~ /\-j/){
		$join = 0;
	}
	elsif($ARGV[$var] =~ /\-o/){
		$out=$ARGV[$var+1];
	}
}

#Checks if files exist
if (! -e $file1) {
	print $file1," does not exist","\n";
}

if (! -e $file2) {
	print $file2," does not exist","\n";
}

#If an output file is defined, write to that, print to STDOUT by default
if (defined $out) {
	open($fh3,'>>',$out);
}
open(my $fh1,'<',$file1);
open(my $fh2,'<',$file2);
#check if file is sorted by coordinates or if there are interspersed lines 
my %seen;
my %chroms;
my %seen1;
my %chroms1;
my @lines=<$fh1>;
my @lines1=<$fh2>;
my $prev=(split /\t/, $lines[0])[0];
my $prev1=(split /\t/,$lines1[0])[0];
print "Checking if entries are sorted.....","\n";
foreach my $x (@lines) {
	chomp $x;
	my $chrom=(split /\t/, $x)[0];
	my $pos=(split /\t/, $x)[1];
	if (! $chroms{$chrom}){
		$chroms{$chrom}=$pos;
	}
	else{
		if ($pos < $chroms{$chrom}) {
			print $file1," not sorted by co-ordinates","\n";
		}
		else{
			$chroms{$chrom}=$pos;
		}
	}
	if ($chrom !~ /^$prev$/) {
		if (! $seen{$prev}) {
			$seen{$prev}=1;
			if ($seen{$chrom}) {
				print "Chromosomes appear to be interspersed. Does ",$file1," need to be sorted?","\n";
			}
		$prev=$chrom;
	}
}
}
#checking if entries are interspersed
foreach my $x1 (@lines1) {
	chomp $x1;
	my $chrom1=(split /\t/, $x1)[0];
	my $pos1=(split /\t/, $x1)[1];
	if (! $chroms1{$chrom1}){
		$chroms1{$chrom1}=$pos1;
	}
	else{
		if ($pos1 < $chroms1{$chrom1}) {
			print $file2," not sorted by co-ordinates","\n";
		}
		else{
			$chroms1{$chrom1}=$pos1;
		}
	}
	if ($chrom1 !~ /^$prev1$/) {
		if (! $seen1{$prev1}) {
			$seen1{$prev1}=1;
			if ($seen1{$chrom1}) {
				print "Chromosomes appear to be interspersed. Does the ",$file2," need to be sorted?","\n";
			}
		$prev1=$chrom1;
	}
}
}
#Get unique chromosome values and for each chromosome and compute the members that overlap
foreach my $test(sort keys %chroms){
print "Parsing ",$test,"\n";
	seek $fh1,0,0;
	seek $fh2,0,0;
	my @x=[];
	my @y=[];
	my $co=0;

	while (my $row1= <$fh1>) {
		chomp $row1;
		my @arr=split /\t/,$row1;
		if ($arr[0] =~ /^$test$/) {
			$x[$co][0] = $arr[1];
			$x[$co][1] = $arr[2];
			$co=$co+1;
			} 								
	}
	my $co1=0;
	while (my $row2= <$fh2>) {
			chomp $row2;
		   my @arrb=split /\t/,$row2;
		   if ($arrb[0] =~ /^$test$/) {
		   	$y[$co1][0] = $arrb[1];
		   	$y[$co1][1] = $arrb[2];
		   	$co1=$co1+1;
		    } 								
	}

my $len1=scalar @x;
my $len2=scalar @y;
my $x1=0;
my $y1=0;
#Run loop for each chromosome while the members exist
while ($x1 < $len1 and $y1 < $len2){
		#print "Comparing ",$x[$x1][0]," ",$x[$x1][1]," and ",$y[$y1][0]," ",$y[$y1][1],"\n";
		if ($y[$y1][0] > $x[$x1][1]) {
			$x1=$x1+1;
		}
		elsif($x[$x1][0] > $y[$y1][1]){
			$y1=$y1+1;
		}
			else{
		my $overlap;
		my $c=$y1;
		if (! defined $join) {
 			$cut=0;
 			}
#If an overlap is found, proceed downward in the file until the start co-ordinate
#of the second file is greater than the stop co-ordinate of the first. Calculate 
#all overlaps in this window and then increment the TE file by one while keeping
#the Intron file same
		while($y[$c][0] and ($x[$x1][1] >= $y[$c][0])){	
		$overlap = min($x[$x1][1],$y[$c][1])-max($x[$x1][0],$y[$c][0]);		
		
		if ($overlap > 0){
			my $temp = 	$x[$x1][1] - $x[$x1][0];
				if ($overlap/$temp >= $perc and defined $join){
					#print "Comparing ",$x[$x1][0]," ",$x[$x1][1]," and ",$y[$c][0]," ",$y[$c][1],"\n";
 					if (defined $out) {
 						print $fh3 $test,"\t",$x[$x1][0],"\t",$x[$x1][0]+$overlap,"\t",$test,"\t",$y[$c][0],"\t",$y[$c][1],"\n";
 						}	
 					print $test,"\t",$x[$x1][0],"\t",$x[$x1][0]+$overlap,"\t",$test,"\t",$y[$c][0],"\t",$y[$c][1],"\n";
 					}
 				elsif($overlap/$temp >= $perc and ! defined $join){
 					#print "Comparing ",$x[$x1][0]," ",$x[$x1][1]," and ",$y[$c][0]," ",$y[$c][1],"\n";
 					if (defined $out) {
 					print $fh3 $test,"\t",$x[$x1][0],"\t",$x[$x1][0]+$overlap,"\n";	
 					}
 					print $test,"\t",$x[$x1][0],"\t",$x[$x1][0]+$overlap,"\n";
 					$cut=$cut+1;	
 					}	
			}
		if(defined $cut and $cut !=0){			
 			last;
 			}
 		else{
 			$c=$c+1;
 			}			
		}
 		$x1=$x1+1;		
		}		
	}
}