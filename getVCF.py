import os
import subprocess as sp
import argparse as ap
import random
import string
import pysam
import glob

#Read in command line arguments
parser = ap.ArgumentParser()
parser.add_argument("-r","--reference",required=True,help ="Add path to reference ",type = str )
parser.add_argument("-p1","--parent1",required=True,help="Add path to parent Read1.fq",type= str)
parser.add_argument("-p2","--parent2",required=True,help="Add path to parent Read2.fq",type= str)
parser.add_argument("-d1","--daughter1",required=True,help="Add path to daughter Read1.fq",type= str)
parser.add_argument("-d2","--daughter2",required=True,help="Add path to daughter Read2.fq",type= str)
parser.add_argument("-g","--gff",required=True,help="Add path to reference GFF",type=str)
parser.add_argument("-retain",action='store_true',help="Retain intermediate files")
args=parser.parse_args()

#Create temp dir to hold processing files
dirname = 'run'+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(4))
sp.call(['mkdir',dirname])
tpath=dirname+'/'

#assumes all binaries are already in path

def getVariants(ref,seq1,seq2,out_name):
     if not(os.path.isfile(ref) or os.path.isfile(seq1) or os.path.isfile(seq2)):
         print('Missing file, please recheck entries')
#     #index the genome
     sp.call(['bwa','index',ref])

     #align the reads
     align = sp.call(['bwa','mem',ref,seq1,seq2,'-o',tpath+out_name+'.sam'])
     print('Done aligning')
     #sort and convert to BAM
     pysam.sort("-O","bam","-o",tpath+out_name+'.bam',tpath+out_name+'.sam')
     print("done sorting")

     #call mpileup
     mpileup = sp.Popen(['bcftools','mpileup','-Ou','-f',ref,tpath+out_name+'.bam'],stdout=sp.PIPE)
     vcfc = sp.Popen(['bcftools','call','-vmO','z','-o',tpath+out_name+'.vcf.gz'],stdin=mpileup.stdout)
     mpileup.wait()
     print("done calling variants")
     return tpath+out_name+'.vcf.gz'
#
def getConsensus(ref,varpath,consensus):
     if not(os.path.isfile(varpath) or os.path.isfile(ref)):
         print("missing file, please recheck entries")
#
#     #index vcf

     sp.call(['tabix','-p','vcf',varpath])
     print("done indexing")
#
#     #get consensus sequence
     sp.call(['bcftools','consensus','-f',ref,varpath,'-o',tpath+consensus])
     print("done getting consensus")
     return tpath+consensus

def annotate(gff,vcf_to_be_annotated):
    sep='\t'
    #extract relevant fields from gff sort GFF
    with open(tpath+'temp.sorted.gff','w') as f:
        sp.call(['bedtools','sort','-i',gff],stdout=f)
    f.close

    if not(os.path.isfile(gff)):
        print("missing file, please recheck entries")
    else:
        with open(tpath+'ann.bed','w') as bed:
            with open(tpath+'temp.sorted.gff') as gff:
                for each in gff:
                    if not((each[0] == '#') or (each.split()[2] == 'region')):
                        out=bed.write(each.split()[0]+sep+each.split()[3]+sep+each.split()[4]+sep+each.split()[6]+sep+each.split()[8]+'\n')
    bed.close()
    gff.close()

    with open(tpath+'ann.hdr','w') as f:
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="GENE">'+'\n')
        f.write('##INFO=<ID=STR,Number=1,Type=String,Description="STR">'+'\n')
    f.close()
    #bgzip and index bed for annotation
    print('indexing for annotation')
    sp.call(['bgzip',tpath+'ann.bed'])
    sp.call(['tabix','-p','bed',tpath+'ann.bed.gz'])
    print("preparing to annotate file")
    #prepare header to add to parent_vcf


    ##annotate vcf
    sp.call(['bcftools','annotate','-a',tpath+'ann.bed.gz','-h',tpath+'ann.hdr','-c','CHROM,FROM,TO,STR,GENE','-o','ann.vcf',vcf_to_be_annotated])

    ##cleanup
def clean():
    if(args.retain):
        pass
    else:
        indpath=args.reference+'.*'
        sp.call(['rm','-r',tpath])
        for f in glob.glob(indpath):
            os.remove(f)

parent_vcf = getVariants(args.reference,args.parent1,args.parent2,'parent')
consensus_seq = getConsensus(args.reference,parent_vcf,'consensus.fa')
daughter_vcf = getVariants(consensus_seq,args.daughter1,args.daughter2,'daughter')
annotate(args.gff,daughter_vcf)
clean()
