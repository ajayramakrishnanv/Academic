#!/usr/bin/perl -w
use strict;
my $file1=$ARGV[0];
my $file2=$ARGV[1];
#Open files to read sequences
open(my $fh1,'<',$file1);
open(my $fh2,'<',$file2);
my $a='';
my $b='';
while (my $row1 = <$fh1>) {
	chomp $row1;
	if ($row1 =~ /^[^>]/) {
		$a=$a.''.$row1;
	}
}
while (my $row2 = <$fh2>) {
	chomp $row2;
	if ($row2 =~ /^[^>]/) {
		$b=$b.''.$row2;
	}
}
#Initialize matrix based on sequence length
my @matrix;
$matrix[0][0]=0;
for (my $x = 1; $x < length($a)+1; $x++) {
	$matrix[0][$x]=$matrix[0][$x-1]-1;
}
for (my $y = 1; $y < length($b)+1; $y++) {
	$matrix[$y][0]=$matrix[$y-1][0]-1;
	#print $matrix[$y][0];
}
#Fill matrix with greatest scores for each base position
my $startx=1;my $starty=1;my $score;
foreach my $one(split //,$b){
	foreach my $two(split //, $a){
		if ($one =~ /$two/) {
			$score=1;
		 } 
		else{
			$score=-1;
		}
		my $score1=$score+$matrix[$startx-1][$starty-1];
		my $score2=$score+$matrix[$startx][$starty-1];
		my $score3=$score+$matrix[$startx-1][$starty];	

		my $greatest=$score1;
		if ($score2 > $greatest) {
			$greatest=$score2;
		}
		elsif($score3>$greatest){
			$greatest=$score3;	
			}
		$matrix[$startx][$starty]=$greatest;	
		#print $one," ",$two," ",$score1," ",$score2," ",$score3," ",$greatest,"\n";
		$starty=$starty+1;
		}
$starty=1;	
#print "\n";	
$startx=$startx+1;	
}
#Traceback, store path contents in an array (Partly because of my stupidity)
my $xb=length($b);my $yb=length($a);
my @seqa;
my @seqb;

while ($xb != 0 && $yb !=0 ) {
#print "Point on ",$xb," ",$yb,"\n";
my $score1b=$matrix[$xb-1][$yb-1];
my $score2b=$matrix[$xb-1][$yb];
my $score3b=$matrix[$xb][$yb-1];

my $greatest=$score1b;
if ($score2b > $greatest) {
	$greatest=$score2b;# body...
	}
elsif($score3b>$greatest){
	$greatest=$score3b;	
	}

if ($score1b==$greatest) {
		unshift(@seqa,$yb);
		unshift(@seqb,$xb);
		$xb=$xb-1;
		$yb=$yb-1;
}
elsif($score2b==$greatest){
		unshift(@seqa,"-");
		unshift(@seqb,$xb);
		$xb=$xb-1;
		$yb=$yb;
}
elsif($score3b==$greatest){
		unshift(@seqa,$yb);
		unshift(@seqb,"-");
		$xb=$xb;
		$yb=$yb-1;
}
#print $greatest," ",$xb," ",$yb,"\n";
}
#Fill gaps if sequences are vastly dissimilar
while($yb!=0){
	unshift(@seqb,"-");
	unshift(@seqa,$yb);
	$yb=$yb-1;
}

while($xb!=0){
	unshift(@seqa,"-");
	unshift(@seqb,$xb);
	$xb=$xb-1;
}
#Initialize arrays to hold print strings
my @mid;
my @pra;my @prb;

foreach my $x(@seqa){
	if($x =~ /\-/){
		push(@pra,$x);
	}
	else{
	my $sub1 =(substr($a,$x-1,1));	
	push(@pra,$sub1);
}
}
foreach my $y(@seqb){
	if($y =~ /\-/){
		push(@prb,$y);
	}
	else{
	my $sub2=(substr($b,$y-1,1));	
	push(@prb,$sub2);
}
}
#Calculate alignment score and print
my $aln=0;
for (my $x = 0; $x < scalar @pra; $x++) {
	if ($pra[$x] =~ $prb[$x]) {
		push(@mid,"|");
		$aln=$aln+1;
	}
	else{
		push(@mid," ");
		$aln=$aln-1;
	}
}
foreach my $a(@pra){print $a;} print "\n"; foreach my $b(@mid){print $b;} print "\n";foreach my $c(@prb){print $c;} print "\n";
print "Alignment score is : ",$aln,"\n";