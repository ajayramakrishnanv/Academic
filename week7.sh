#!/bin/bash 

#GETOPTS BLOCK####################################################################################################

usage="Usage : week7.sh -a -b -r -o -e -z -i -v -h \n
		\t-a [Read pair one] \n
		\t-b [Read pair	 two] \n
		\t-r [reference] \n
		\t-o [output] \n 
		\t-e [realign] \n
		\t-z [zipped output] \n
		\t-i [index bam] \n
		\t-v [verbose] \n
		\t-h [help]"

while getopts ":a:b:r:o:eczivh" opt
do
  case $opt in
    a)	
	f1=$OPTARG
	;;		
    b)
        f2=$OPTARG 
	;;
    r)
	ref=$OPTARG
	;;
    o)
	out=$OPTARG	
	;;
    e)
	realign=1
	;;
    i)
	index=1
	;;
    z)
	zipped=1
	;;
    v)
	ver=1 
	;;
    ?)
	echo -e $usage
	exit 
	;;			 		
    h) 
	echo -e $usage
	exit
	;;			
  esac
done
shift $[ $OPTIND-1 ]

#Check Block###################################################################################################
#Checking for file existances and prompting if absent. I have accounted for file names to be supplied as paths as well.
if [ -z $f1 ] || [ -z $f2 ] || [ -z $ref ] || [ -z $out ]
then
	echo -e $usage
exit
fi

if  [[ ! ( -f ./data/$f1 || -f $f1 ) ]] 
then	
	echo "$f1 does not exist or is not specified"
exit 2
fi


if  [[ ! ( -f $f2  || -f ./data/$f2 ) ]] 
then	
	echo "$f2 does not exist"
exit 2
fi

if  [[ ! ( -f $ref  || -f ./data/$ref ) ]] 
then	
	echo "$ref does not exist"
exit 2
fi

if  [[ ( -f $out.vcf  || -f ./output/$out.vcf || -f $out || -f ./output/$out ) ]] 
then	
	read -p "The file already exists, do you want to overwrite?(Y to ovrwrite/N to exit): " name
	if [[ $name == "Y" ]]
	then
	     :
	else
	     exit
	fi
fi

if [ -f ./data/$f1 ]
then
	f1=./data/$f1 
fi
if [ -f ./data/$f2 ] 
then
	f2=./data/$f2
fi
if [ -f ./data/$ref ]
then
	ref=./data/$ref	
fi

if [[ $out =~ .*\.vcf.gz$ ]] 
then	
	out=`echo $out | sed -n 's/\(.*\.vcf\\)\.gz/\1/p'`
fi

if ! [[ $out =~ .*\.vcf$ ]] 
then	
	out=$out.vcf
fi	


if [[ $out =~ .*/.* ]]
then	
	:
else
	out=./output/$out
fi

##Check Block###################################################################################################
#Initializing a to extract part of reference file to create dictionary (E.g : 'chr17' in 'chr17.fa' )
a=`echo $ref | sed -n 's/.*\(\/.*\)\(\.fa$\)/\1/p'| cut -b 2-`
################################################################################################################

##Alignment

SM="Test"

	if [ ! -z  $ver ]
	then
	 	echo -e "\nCreating Alignment\n" 
	fi
	
	bwa index $ref
	#./lib/bwa index $ref ##to be used if for some reason bwa is not in the path of the running system. Likewise for 	the rest

	if [ ! -z $ver ]
	then
	 	echo -e "\nGenerating SAM file\n" 
	fi	

	bwa mem -R '@RG\tID:Read\tSM:Samp\tLB:Lib' $ref $f1 $f2 > ./tmp/$SM.sam 
		

	if [ ! -z $ver ]
	then
	 	echo -e "\nFixing the SAM file\n" 
	fi

	
	samtools fixmate -O bam ./tmp/$SM.sam ./tmp/$SM.fixed.bam
	 
	if [ ! -z $ver ]
	then
	 	echo -e "\nSorting the BAM file\n" 
	fi
	
	
	samtools sort -O bam -o ./tmp/$SM.sorted.bam -T ./tmp/$SMtmp ./tmp/$SM.fixed.bam

if [ ! -z $ver ]
then
	echo -e "\nIndexing and alignment completed.........................\n"
fi

if [[ $index == 1 && $realign != 1  ]]
then
	
	samtools index ./tmp/$SM.sorted.bam
	mv ./tmp/$SM.sorted.bam.bai ./output/output.bam.bai
fi

if [ ! -z $ver ]
then
	echo -e "\nYour file has been indexed.........................\n"
fi

###############################################################################################################

##Improvement
##includes duplicate marking
if [[ $realign == 1 ]]
then
	a=`echo $ref | sed -n 's/.*\(\/.*\)\(\.fa$\)/\1/p'| cut -b 2-`
	
	if [! -z $ver ]
	then
	 	echo -e "\nCreating reference FASTA index\n" 
	fi

	
	samtools faidx $ref
 	
	if [ ! -z $ver ]
	then
	 	echo -e "\nBAM index for re alignment index\n" 
	fi

	
	samtools index ./tmp/$SM.sorted.bam
	
	if [ ! -z $ver ]
	then
	 	echo -e "\nCreating reference Sequence Dicitonary\n" 
	fi
	
	java -jar ./lib/picard.jar CreateSequenceDictionary REFERENCE=$ref OUTPUT=./data/$a.dict
	
	if [ ! -z $ver ]
	then
	 	echo -e "\nInitializing re-alignment\n" 
	fi

	java -Xmx2g -jar ./lib/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I ./tmp/$SM.sorted.bam -o ./tmp/$SM.intervals -known ./data/Mills1000G.b38.vcf
#### The file I used as a reference was Mills_and_1000G_gold_standard.indels.hg38.vcf available in my lib
	java -Xmx4g -jar ./lib/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I ./tmp/$SM.sorted.bam -targetIntervals ./tmp/$SM.intervals -known ./data/Mills1000G.b38.vcf -o ./tmp/$SM.realigned.bam
	
	if [ ! -z $ver ]
	then
	 	echo -e "\nMarking Duplicates\n" 
	fi	

	java -jar ./lib/picard.jar MarkDuplicates I=./tmp/$SM.realigned.bam O=./tmp/$SM.realigned1.bam M=./tmp/Dups.txt
	
	
	samtools index ./tmp/$SM.realigned1.bam

	 if [ ! -z $ver ]
	then
	 	echo -e "\nIndexing File\n" 
	fi
	
	if [[ $index == 1 ]]
	then
	mv ./tmp/$SM.realigned1.bam.bai ./output/output.bam.bai
	fi
fi


###########################################################################################IMPROVEMENT SECTION

##########################VARIANT CALLING#####################################################################
if [ ! -z $ver ]
then
	echo -e "\nGenerating bcf file for calling variants.....\nGenerating output VCF file....\n"
fi

if [[ $realign == 1 ]]
then
	
	samtools mpileup -go ./tmp/output.bcf -f $ref ./tmp/$SM.realigned1.bam
	bcftools call --skip-variants indels --multiallelic-caller --variants-only -Ov -o $out ./tmp/output.bcf
	mv ./tmp/$SM.realigned1.bam ./output/output.bam
	
else
	samtools mpileup -go ./tmp/output.bcf -f $ref ./tmp/$SM.sorted.bam
	bcftools call --skip-variants indels --multiallelic-caller --variants-only -Ov -o $out ./tmp/output.bcf
	mv ./tmp/$SM.sorted.bam ./output/output.bam
fi

#############################COMPRESS OUTPUT IF SPECIFIED#####################################################
if [ ! -z $ver ]
then
	echo -e "\nCompressing output..Reorganizing files\n" 
fi


if [[ $zipped == 1 ]]
then
	gzip $out
fi

##Cleaning everything out from ./data
ls ./data/*.fa.* | xargs -I {} mv {} ./tmp/

if [ -f ./data/$a.dict ]
then
mv ./data/$a.dict ./tmp/
fi 
