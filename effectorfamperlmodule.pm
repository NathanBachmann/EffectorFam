#################
#
#A list of perl subroutines for EffectorFAM cgi script
#
################

# revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}

# get_frame
#
# A subroutine to translate a frame of DNA

sub get_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
        return translate_dna ( substr ( $seq, $start - 1, $end -$start + 1) );
}

#translate_dna
#
#translate the dna sequences

sub translate_dna {

	my($seq) = @_;
	my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'tca' => 'S',
    'TCC' => 'S',    # Serine
    'tcc' => 'S',
    'TCG' => 'S',    # Serine
    'tcg' => 'S',
    'TCT' => 'S',    # Serine
    'tct' => 'S',
    'TTC' => 'F',    # Phenylalanine
    'ttc' => 'F',
    'TTT' => 'F',    # Phenylalanine
    'ttt' => 'F',
    'TTA' => 'L',    # Leucine
    'tta' => 'L',
    'TTG' => 'L',    # Leucine
    'ttg' => 'L',
    'TAC' => 'Y',    # Tyrosine
    'tac' => 'Y',
    'TAT' => 'Y',    # Tyrosine
    'tat' => 'Y',
    'TAA' => '*',    # Stop
    'taa' => '*',
    'TAG' => '*',    # Stop
    'tag' => '*',
    'TGC' => 'C',    # Cysteine
    'tgc' => 'C',
    'TGT' => 'C',    # Cysteine
    'tgt' => 'C',
    'TGA' => '*',    # Stop
    'tga' => '*',
    'TGG' => 'W',    # Tryptophan
    'tgg' => 'W',
    'CTA' => 'L',    # Leucine
    'cta' => 'L',
    'CTC' => 'L',    # Leucine
    'ctc' => 'L',
    'CTG' => 'L',    # Leucine
    'ctg' => 'L',
    'CTT' => 'L',    # Leucine
    'ctt' => 'L',
    'CCA' => 'P',    # Proline
    'cca' => 'P',
    'CCC' => 'P',    # Proline
    'ccc' => 'P',
    'CCG' => 'P',    # Proline
    'ccg' => 'P',
    'CCT' => 'P',    # Proline
    'cct' => 'P',
    'CAC' => 'H',    # Histidine
    'cac' => 'H',
    'CAT' => 'H',    # Histidine
    'cat' => 'H',
    'CAA' => 'Q',    # Glutamine
    'caa' => 'Q',
    'CAG' => 'Q',    # Glutamine
    'cag' => 'Q',
    'CGA' => 'R',    # Arginine
    'cga' => 'R',
    'CGC' => 'R',    # Arginine
    'cgc' => 'R',
    'CGG' => 'R',    # Arginine
    'cgg' => 'R',
    'CGT' => 'R',    # Arginine
    'cgt' => 'R',
    'ATA' => 'I',    # Isoleucine
    'ata' => 'I',
    'ATC' => 'I',    # Isoleucine
    'atc' => 'I',
    'ATT' => 'I',    # Isoleucine
    'att' => 'I',
    'ATG' => 'M',    # Methionine
    'atg' => 'M',
    'ACA' => 'T',    # Threonine
    'aca' => 'T',
    'ACC' => 'T',    # Threonine
    'acc' => 'T',
    'ACG' => 'T',    # Threonine
    'acg' => 'T',
    'ACT' => 'T',    # Threonine
    'act' => 'T',
    'AAC' => 'N',    # Asparagine
    'aac' => 'N',
    'AAT' => 'N',    # Asparagine
    'aat' => 'N',
    'AAA' => 'K',    # Lysine
    'aaa' => 'K',
    'AAG' => 'K',    # Lysine
    'aag' => 'K',
    'AGC' => 'S',    # Serine
    'agc' => 'S',
    'AGT' => 'S',    # Serine
    'agt' => 'S',
    'AGA' => 'R',    # Arginine
    'aga' => 'R',
    'AGG' => 'R',    # Arginine
    'agg' => 'R',
    'GTA' => 'V',    # Valine
    'gta' => 'V',
    'GTC' => 'V',    # Valine
    'gtc' => 'V',
    'GTG' => 'V',    # Valine
    'gtg' => 'V',
    'GTT' => 'V',    # Valine
    'gtt' => 'V',
    'GCA' => 'A',    # Alanine
    'gca' => 'A',
    'GCC' => 'A',    # Alanine
    'gcc' => 'A',
    'GCG' => 'A',    # Alanine
    'gcg' => 'A',
    'GCT' => 'A',    # Alanine
    'gct' => 'A',
    'GAC' => 'D',    # Aspartic Acid
    'gac' => 'D',
    'GAT' => 'D',    # Aspartic Acid
    'gat' => 'D',
    'GAA' => 'E',    # Glutamic Acid
    'gaa' => 'E',
    'GAG' => 'E',    # Glutamic Acid
    'gag' => 'E',
    'GGA' => 'G',    # Glycine
    'gga' => 'G',
    'GGC' => 'G',    # Glycine
    'ggc' => 'G',
    'GGG' => 'G',    # Glycine
    'ggg' => 'G',
    'GGT' => 'G',    # Glycine
    'ggt' => 'G',
    );
    my $i;
    my $codon;
    my $protein;
    
    
    for (my $i=0; $i < (length($seq) - 2) ; $i += 3)
	{
		$codon = substr($seq, $i, 3);
		if (exists $genetic_code{$codon})
		{
			$protein .= $genetic_code{$codon};
		}
	}
	
	return $protein;
}

#frames2fasta
#
#Converts ORFs to fasta format

sub frames2fasta
{
	my($seq, $start, $frame) = @_;
	my $head;
	my $element;
	my $field;
	my $orf;
	my $last = 0;
	my $first = 0;
	my $base;
	my $count = 1;
	
	my @part = split(/\*/, $seq);
	foreach $element (@part)
	{
		my $size = length($element);
		
		$base = ($size * 3) + 3;
		
		$first = $last + 1 + $start;
		$last = $last + $base;
		
		#prints the orf if they are larger than 150 amino acids
		if ($size > '150')
		{
			$field = ">RF_$frame\_ORF_$count $first..$last\n$element\n";
			$orf .= $field;
			$count++;
		}
		
	}

	return $orf;
}

#frames2fasta_revcom
#
#Coverts ORF on reverse strand to fasta format

sub frames2fasta_revcom
{
	my($aa, $start, $total, $frame) = @_;
	my $element;
	my $high = 0;
	my $low = 0;
	my $field;
	my $orf;
	my $base;
	my $count = "1";
	
	my @part = split(/\*/, $aa);
	foreach $element (@part)
	{
		my $size  = length($element);
		
		$base = ($size * 3) + 2;
		
		$high = $total - $start;
		$low = $high - $base;
		
		#check is the low number is negative number or zero
		if (($low =~ /-\d/) or ($low =~ /^0$/))
		{
			$low = 1;
		}

		
		#prints the ORFs if they are larger than 150 amino acids
		if ($size > '150')
		{
			$field = ">RF_$frame\_ORF_$count $high..$low\n$element\n";
			$orf .= $field;
			$count++;
		}
		
		$total = ($low + $start) - 1;
	}	
	
	return $orf;
}

#dir_size
#
#Calculates current size of the directory

sub dir_size 
{
	my $dir = shift;
	my $dir_size = 0;
	opendir (DIR, $dir) or die "Unable to open $dir\n";
	foreach (readdir DIR)
	{
		$dir_size += -s "$dir/$_";
	}
	return $dir_size;
}

1;