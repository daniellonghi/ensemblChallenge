#!/usr/bin/perl
use Bio::DB::Fasta;
use Path::Tiny;
use strict;

#currentGFF - This file must be a GFF3, its file will be loaded and 
#the new coordinates will starting to be search in a older assembly
my $currentGFF = $ARGV[0];
my $currentAssembly = $ARGV[1];
my $olderAssembly = $ARGV[2];

#REMOVE ALL LINES WITH STARTING WITH #
#system("sed -i '' '/^#/ d' " . $currentGFF);
#FIST LINE IS CHR SIZE
#system("sed -i '' '1d' " . $currentGFF);

open(my $file1, "<$currentGFF") or die "Please, check your file.";

system("rm ./newCoordinatesBlast.gff3");
open(my $newGeneCoordinates, ">:encoding(UTF-8)", "./newCoordinatesBlast.gff3");

print("#######\n");
print("STARTING PROCESS\n");

#GLOBAL VARIABLE
my $lineRecords = "";
my $returnBLAST = "";
my $returnFASTA = "";
my $returnWrite = "";
my $newOutputAnnotation = "";
my @oldCoordinates;
my @newCoordinates;
#GLOBAL VARIABLE
my $i = 0;
print("CONVERTING\n");
while ($lineRecords = <$file1>){
	print("Sequence: " . ($i+1) . "\n");
	chomp $lineRecords;
	#SAVE COORDINATES INTO TEMP FILE
	$returnWrite = system('echo "' . $lineRecords . '" > tempSubquery.gff3');
	@oldCoordinates = split(/\t/, $lineRecords);

	#CALL FOR EXTRACT FASTA AND SAVE TEMP FILE
	$returnFASTA = system('bedtools getfasta -fi ./' . $currentAssembly . ' -bed ./tempSubquery.gff3 -name+ > tempSubquery.fa');
	
	#CALL FOR BLAST
	$returnBLAST = system('blastn -query ./tempSubquery.fa -subject ./' . $olderAssembly . ' -outfmt "6 sseqid sstart send pident" -max_target_seqs 1 -max_hsps 1 -out tempBlast.txt');
	
	#GET FIRST LINE OF THE NEW COORDINATE
	@newCoordinates = path('./tempBlast.txt')->lines( { chomp => 1, count => 1 } );
	
	#CHECK IF A HIT WAS FOUND
	if (scalar @newCoordinates >= 1){
		@newCoordinates = split(/\t/, @newCoordinates[0]);
		#REPLACE START-END COORDINATE
		@oldCoordinates[3] = @newCoordinates[1];
		@oldCoordinates[4] = @newCoordinates[2];
		@oldCoordinates[8] = "NCI=" . @newCoordinates[3] . ";" . @oldCoordinates[8];
	}else{
		#REPLACE START-END COORDINATE
		@oldCoordinates[3] = "___";
		@oldCoordinates[4] = "___";
		@oldCoordinates[8] = "NCI=0;" . @oldCoordinates[8];
	}

	$newOutputAnnotation .= @oldCoordinates[0] . "\t" . @oldCoordinates[1] . "\t" . @oldCoordinates[2] . "\t";
	$newOutputAnnotation .= @oldCoordinates[3] . "\t" . @oldCoordinates[4] . "\t" . @oldCoordinates[5] . "\t";
	$newOutputAnnotation .= @oldCoordinates[6] . "\t" . @oldCoordinates[7] . "\t" . @oldCoordinates[8];
	$newOutputAnnotation .= "\n";
	$i = $i + 1;
}
print("CONVERTING DONE\n");
print("FINISHING PROCESS\n");
print("#######\n");

print $newGeneCoordinates $newOutputAnnotation;
close($newGeneCoordinates);
system("rm temp*");

exit 0;
