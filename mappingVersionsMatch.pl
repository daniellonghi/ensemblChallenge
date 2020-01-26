#!/usr/bin/perl
use Bio::DB::Fasta;

#currentGFF - This file must be a GFF3, its file will be loaded and 
#the new coordinates will starting to be search in a older assembly.
my $currentGFF = $ARGV[0];
my $curretAssembly = $ARGV[1];
my $oldAssembly = $ARGV[2];

#REMOVE ALL LINES WITH STARTING WITH #
#system("sed -i '' '/^#/ d' " . $currentGFF);
#FIST LINE IS CHR SIZE
#system("sed -i '' '1d' " . $currentGFF);

my $FileGFF;
my $currentGenome;
my $olderGenome;

open($FileGFF, "<$currentGFF") or die "Please, check your file.";
open($currentGenome, "<$curretAssembly") or die "Please, check your file.";
open($olderGenome, "<$oldAssembly") or die "Please, check your file.";

#PUT IN MEMORY DNA SEQUENCE
my $sFile;
my @idSequence;
my $id;
$sFile = Bio::DB::Fasta->new("./" . $oldAssembly);
@idSequence = $sFile->get_all_primary_ids;

my $fullHeader;
my $positiveDNA;
my $negativeDNA;

print("#######\n");
print("READING DNA\n");
foreach $id (@idSequence){
	$fullHeader = $sFile->header($id);
	$positiveDNA = $sFile->get_Seq_by_id($id);
	$positiveDNA = uc $positiveDNA->seq;
}

# $negativeDNA = reverse $positiveDNA;
# #MAKE SURE TO DO NOT REPLACE T for A and C for G.
# $negativeDNA =~ tr/A/X/;
# $negativeDNA =~ tr/T/A/;
# $negativeDNA =~ tr/X/T/;
# $negativeDNA =~ tr/G/X/;
# $negativeDNA =~ tr/C/G/;
# $negativeDNA =~ tr/X/C/;

print("#######\n\n");

#GLOBAL VARIABLE
my $lineRecord;
my $lineRecordFormatted;
my @lineSplitted;
my $fullSeqToSplit;
my @dnaSplitted;
my $outPutAnnotationGene = "";
#GLOBAL VARIABLE

print("#######\n");
print("STARTING PROCESS REMAPPING\n");
my $i = 1;
foreach $lineRecord (<$FileGFF>){
	print("Formatting Sequence: ". $i . "\n");
	#FORCE REBUILD DB::FASTA BECASE OF TEMP FILES
	system("rm *.index");

	#FORMATTING LINE AS RECORD
	chomp $lineRecord;
	@lineSplitted = split(/\t/, $lineRecord);
	$lineRecordFormatted = @lineSplitted[0] . "\t" . @lineSplitted[1] . "\t" . @lineSplitted[2] . "\t" . @lineSplitted[3] . "\t" . @lineSplitted[4] . "\t";
	$lineRecordFormatted .= @lineSplitted[5] . "\t" . @lineSplitted[6] . "\t" . @lineSplitted[7] . "\t" . @lineSplitted[8]; 

	#SAVE SEQUENCE
	system("echo '" . $lineRecordFormatted . "' > tempSeq.gff3");

	#GET FASTA OF THE SEQUENCE
	system("bedtools getfasta -fi ./" . $curretAssembly . " -bed ./tempSeq.gff3 > tempSeq.fa");

	#READ THE FASTA SEQUENCE TO SPLIT
	$sFile = Bio::DB::Fasta->new("./tempSeq.fa");
	@idSequence = $sFile->get_all_primary_ids;

	if (scalar @idSequence >= 1){
		foreach $id (@idSequence){
			$fullHeader = $sFile->header($id);
			$fullSeqToSplit = $sFile->get_Seq_by_id($id);
			$fullSeqToSplit = uc $fullSeqToSplit->seq;
		}
	}

	#STARTING SPLIT TO GET NEW COORDINATES
	# if(@lineSplitted[6] == "-" || @lineSplitted[6] == "."){
	# 	@dnaSplitted = split(/$fullSeqToSplit/, $negativeDNA);
	# }elsif (@lineSplitted == "+"){
	# 	@dnaSplitted = split(/$fullSeqToSplit/, $positiveDNA);
	# }

	@dnaSplitted = split(/$fullSeqToSplit/, $positiveDNA);

	if (scalar @dnaSplitted >= 1){
		@lineSplitted[3] = length(@dnaSplitted[0]);
		@lineSplitted[4] = length(@dnaSplitted[0]) + length($fullSeqToSplit);

		$outPutAnnotationGene .= @lineSplitted[0] . "\t" . @lineSplitted[1] . "\t" . @lineSplitted[2] . "\t" . (@lineSplitted[3] + 1) . "\t" . @lineSplitted[4] . "\t";
		$outPutAnnotationGene .= @lineSplitted[5] . "\t" . @lineSplitted[6] . "\t" . @lineSplitted[7] . "\t" . @lineSplitted[8] . "\n";
	}else{
		@lineSplitted[3] = "____";
		@lineSplitted[4] = "____";

		$outPutAnnotationGene .= @lineSplitted[0] . "\t" . @lineSplitted[1] . "\t" . @lineSplitted[2] . "\t" . (@lineSplitted[3] + 1) . "\t" . @lineSplitted[4] . "\t";
		$outPutAnnotationGene .= @lineSplitted[5] . "\t" . @lineSplitted[6] . "\t" . @lineSplitted[7] . "\t" . @lineSplitted[8] . "\n";
	}
	$i = $i+1;
}

system("echo '" . $outPutAnnotationGene . "' > newCoordinatesMatch.gff3");

print("#######\n");
system("rm temp*");

#print $GFFFamiliaClassifica $GFFOutPut;
#print $FASTAFamiliaClassifica $FASTAOutPut;

#close($GFFFamiliaClassifica);
#close($FASTAFamiliaClassifica);

exit 0;