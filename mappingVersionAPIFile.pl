#!/usr/bin/perl
use Bio::EnsEMBL::Registry;
use Data::Dumper;

#currentGFF - This file must be a GFF3, its file will be loaded and 
#the new coordinates will starting to be search in a older assembly.
my $currentGFF = $ARGV[0];

#REMOVE ALL LINES WITH STARTING WITH #
#system("sed -i '' '/^#/ d' " . $currentGFF);
#FIST LINE IS CHR SIZE
#system("sed -i '' '1d' " . $currentGFF);

my $registry = 'Bio::EnsEMBL::Registry';

print("\n#######");
print("\nSTARTING");
print("\n#######\n");

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 3337
);

my $FileGFF;
open($FileGFF, "<$currentGFF") or die "Please, check your file.";

my $columnEightType;
my @auxSplit;
my $subTypeGene;
my $subTypeTranscript;
my $outPutAnnotationGene;
foreach $lineRecord (<$FileGFF>){
	#FORMATTING LINE AS RECORD
	chomp $lineRecord;
	@lineSplitted = split(/\t/, $lineRecord);
	$lineRecordFormatted = @lineSplitted[0] . "\t" . @lineSplitted[1] . "\t" . @lineSplitted[2] . "\t" . @lineSplitted[3] . "\t" . @lineSplitted[4] . "\t";
	$lineRecordFormatted .= @lineSplitted[5] . "\t" . @lineSplitted[6] . "\t" . @lineSplitted[7] . "\t" . @lineSplitted[8]; 

	#GET GENE OR TRANSCRIPT SEQUENCES FROM COLUMN 8th
	$columnEightType = @lineSplitted[8];
		
	$subTypeGene = "ID=gene";
	if (index($columnEightType, $subTypeGene) != -1) {
		
    	@auxSplit = split(";", $columnEightType);
    	@auxSplit = split(":", @auxSplit[0]);
    	$subTypeGene = @auxSplit[1];

		#GET GENE INFORMATION FROM ENSEMBL
		my $gene_adaptor  = $registry->get_adaptor('Human', 'Core', 'Gene');
		my $gene = $gene_adaptor->fetch_by_stable_id($subTypeGene);

		#GET TRANSCRIPT INFORMATION FROM ENSEMBL
		# my $transcript_adaptor  = $registry->get_adaptor('Human', 'Core', 'Transcript');
		# my $transcript = $transcript_adaptor->fetch_by_stable_id($subTypeTranscript);

		if ($gene){
			print("\nFormatting new sequence.");
			print(" -> " . $subTypeGene);

			@lineSplitted[3] = $gene->start;
			@lineSplitted[4] = $gene->end;
			$outPutAnnotationGene = @lineSplitted[0] . "\t" . @lineSplitted[1] . "\t" . @lineSplitted[2] . "\t" . @lineSplitted[3] . "\t" . @lineSplitted[4] . "\t";
			$outPutAnnotationGene .= @lineSplitted[5] . "\t" . @lineSplitted[6] . "\t" . @lineSplitted[7] . "\t" . @lineSplitted[8];
		}
		# SAVE THE RECORDS REMAPPED
		system("echo '" . $outPutAnnotationGene . "' >> newCoordinatesAPIFile.gff3");
	}
}
print("\n#######");
print("\nDONE");
print("\n#######\n");
