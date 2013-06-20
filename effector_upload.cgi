#!/usr/bin/perl -w
#Written by Nathan Bachmann (nathan_bachmann@hotmail.com) under the supervision of Dr Scott Beatson
#(s.beatson@uq.edu.au) on the 13/02/12
use strict;
use CGI;
use CGI::Carp qw ( fatalsToBrowser );  
use File::Basename;
use effectorfamperlmodule;

$CGI::POST_MAX = 1024 * 10000; #limits uplaod to 10 Mb
use constant MAX_DIR_SIZE => 100 * 1_048_576; #limits total directory size to 100 Mb

#variable list
my $model;
my $result;
my @list;
my $line;
my $line2;
my $line4;
my $line5;
my $name;
my $hmm;
my @part;
my $gene;
my $data;
my $bitScore;
my $eValue;
my $N;
my $eValueDom;
my %hash;
my %hash2;
my $line3;
my @piece;
my $desc;
my @element;
my @element2;
my $model2;
my $bit;
my $tag;
my %HoH;
my $key;
my $key2;
my @array;
my $hi_score;
my @effectorfam_res;
my @field;
my $check;
my $aa;
my $seq;
my $size;
my $protein;
my $cds;
my $header1;
my $header2;
my $anno;

#identify the upload directory
my $upload_dir = "/var/www/EffectorFam";
my $img_location = "http://effectorfam-scmb.biosci.uq.edu.au";

#read the form Variables
my $query = new CGI;
my $filename = $query->param("file");
my $txt = $query->param("seq");
my $count = 1;

#check for existing data
foreach (1..10)
{
	if (-e "$upload_dir/upload/enterdata.$count.faa")
	{
		$count++;
		redo;
	}
	else
	{
		last;
	}
}
#safety limits
my $safe_filename_characters = "a-zA-Z0-9_.-"; #ensure specific characters are used in the name

#prints the table headers in the browser
print "Content-type: text/html\n\n";

print <<Endofheader;
<html>
<head>
	<title>Effectorfam Results</title>
</head>

<body>
<img src="$img_location/img/banner_res.gif" height="205" width="1000">
</div>

Endofheader

#checks uplaod file size
if (defined($filename))
{
	print "<p>Your upload file is $filename</p>";
}
else
{
	print "<p>Error 02: Upload file may be too big (Max uplaod is 10 Mb)</p>";
	exit;
}

#checks to male the sure the directory is under 100 Mb
if (dir_size($upload_dir) + $ENV{CONTENT_LENGTH} > MAX_DIR_SIZE)
{
	print "<p>Upload directory is full</p>";
	exit;
}	

#make the filename safe by removing the pth using the filepaser routine from Basename
my ($name, $path, $extension) = fileparse ($filename, '\..*');
$filename = $name . $extension;

$filename =~ tr/ /_/; #replace spaces with underscore
$filename =~ s/[^$safe_filename_characters]//g; #remove any characters that are not in the safe list

#Double check that the filename is safe
if ($filename =~ /^([$safe_filename_characters]+)$/)
{
	$filename = $1;
	#getting the file handle using CG.pm upload function
	my $upload_filehandle = $query->upload("file");
	
	#save the file or prints error if it doesn't work
	open ( UPLOADFILE, ">$upload_dir/upload/enterdata.$count.faa") or die "$!";
	binmode UPLOADFILE; #write the file in binary reduce chances of corruption

	while ( <$upload_filehandle> )
	{
		print UPLOADFILE;
	}

	close UPLOADFILE;
	$filename = "enterdata.$count.faa";
}
elsif ($txt =~ /^>+/)
{
	open (SAVEFILE, ">$upload_dir/upload/enterdata.$count.faa") or die "Couldn't open save file\n";
	print SAVEFILE $txt;
	close SAVEFILE;
	$filename = "enterdata.$count.faa";	
}
else
{
	die "filename contains invalid characters\n";
}


#check to see if the input file is nucleic acid
open (IN2, "$upload_dir/upload/$filename") or die "$!";

foreach $line5 (<IN2>)
{
	chomp $line5;
	if ($line5 =~ ">")
	{
		next;
	}
	elsif ($line5 =~ /[bdefhiklmnpqrsvwyz]/i)
	{
		$aa = 'T';
	}
	else
	{
		$aa = 'F';
		$seq .= $line5;
	}
}

#if sequence is nucleic acid then translate it
if ($aa =~ 'F')
{
	open (NEWFILE, ">$upload_dir/upload/moddata.$count.faa") or die "$!";
	$size = length($seq);

	#translate reading frame 1
	$protein = get_frame($seq, 1);
	$cds = frames2fasta($protein, 0, 1);
	print NEWFILE "$cds";

	#translate reading frame 2
	$protein = get_frame($seq, 2);
	$cds = frames2fasta($protein, 1, 2);
	print NEWFILE "$cds";

	#translate reading frame 3
	$protein = get_frame($seq, 3);
	$cds = frames2fasta($protein, 2, 3);
	print NEWFILE "$cds";

	#get reverse complement of input sequence
	my $rev_seq = revcom($seq);

	#translate reading frame 4
	$protein = get_frame($rev_seq, 1);
	$cds = frames2fasta_revcom($protein, 0, $size, 4);
	print NEWFILE "$cds";

	#translate reading frame 5
	$protein = get_frame($rev_seq, 2);
	$cds = frames2fasta_revcom($protein, 1, $size, 5);
	print NEWFILE "$cds";

	#translate reading frame 6
	$protein = get_frame($rev_seq, 3);
	$cds = frames2fasta_revcom($protein, 2, $size, 6);
	print NEWFILE "$cds";
	
	close NEWFILE;
}
elsif ($aa =~ 'T')
{
	`cp $upload_dir/upload/enterdata.$count.faa $upload_dir/upload/moddata.$count.faa`;
}
	
#perform HMM search
open (IN, "$upload_dir/db/hmm_list.txt"); #opens a list of HMM names

@list = <IN>;

foreach $line3 (@list)	#runs hmmsearch from HMMer3
{
	chomp $line3;
	@piece = split(/\t/, $line3);
	$model = $piece[0];
	$anno = $piece[1];
	$desc = $piece[2];
	$hash{$model} = $desc;
	$hash2{$model} = $anno;
	$result = `hmmsearch $upload_dir/db/$model.hmm $upload_dir/upload/moddata.$count.faa`;
	open (OUT, ">$upload_dir/out/$model.$count.out");
	print OUT "$result";
	close OUT;
}

close IN;

#Parser HMM results

foreach $line (@list)
{
	@element = split(/\t/, $line);
	$model2 = $element[0];
	open (IN2, "$upload_dir/out/$model2.$count.out") or die "Couldn't open $model2\n";
	while ($line2 = <IN2>)
	{
		chomp $line2;
		if ($line2 =~ /^#\squery\sHMM\sfile:\s+(\S+)$/)
		{
			$name = $1;
			@part = split(/\//, $name);
			$hmm = $part[-1];
		}
		elsif ($line2 =~ /Scores for complete sequences/)
		{
			$check = "yes";
		}
		elsif ($line2 =~ /^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(.+)$/)
		{
			$gene = $9;
			$data = $10;
			$bitScore = $2;
			$eValue = $1;
			$N = $8;
			$eValueDom = $4;
			if (($eValue < 0.00001) and ($eValueDom < 0.00001) and ($bitScore > 100))
			{
				push (@effectorfam_res, "$gene\t$data\t$model2\t$hash2{$model2}\t$hash{$model2}\t$bitScore\t$eValue\t$N");
			}
			else
			{
				next;
			}
		}
		else
		{
			next;
		}
	}
	close IN2;
}

#checks the HMMer ran correctly
if ($check =~ /yes/)
{
	$check = "good";
}
else
{
	print "<p> Error 01: You may have invalid characters in your input sequences.<p>";
	print "<p> Support: nathan.bachmann\@uqconnect.edu.au<p>";
	exit;
}


#This section seraches for any proteins that are identfied by multiple HMMs and record which one has the hightest bitscore

foreach (@effectorfam_res)
{
	$line4 = $_;
	@element2 = split(/\t/, $line4);
	$bit = @element2[4];
	$tag = @element2[0];
	$HoH{$tag}{$bit} = $line4;
}	

#creates output file
open (OUT2, ">$upload_dir/upload/db_res.$count.tab");

if ($aa =~ 'T')
{
	$header1 = "Protein";
	$header2 = "Description";
}
elsif ($aa =~ 'F')
{
	$header1 = "ORF";
	$header2 = "Coordinates";
}

print OUT2 "$header1\t$header2\tHMM\tAnnotation\tHMM description\tBit Score\tE Value\tNo. of Domains\n";

#draws the table headers for output page
print <<EndofTable;
<table border=1>
	<tr>
		<th>$header1</th>
		<th>$header2</th>
		<th>HMM</th>
		<th>Annotation</th>
		<th>HMM description</th>
		<th>Bitscore</th>
		<th>E value</th>
		<th>No. of Domains</th>
		<th>Raw Output</th>
	</th>
EndofTable



#finds the duplicates
foreach $key (sort keys %HoH)
{
	@array=();
	foreach $key2 (sort { $b <=> $a } keys %{$HoH{$key}})
	{
		push (@array, $key2);
	}
	$hi_score = $array[0];
	print OUT2 "$HoH{$key}{$hi_score}\n";
	@field = split(/\t/, $HoH{$key}{$hi_score});
	print "<tr><td><b>$field[0]</b></td><td width=\"140\" style=\"white-space:pre-warp\">$field[1]</td><td>$field[2]</td><td width=\"140\" style=\"white-space:pre-warp\">$field[3]</td><td width=\"140\" style=\"white-space:pre-warp\">$field[4]</td><td>$field[5]</td><td>$field[6]</td><td>$field[7]</td><td><a href=\"$img_location/out/$field[2].$count.out\">Alignment</a></td></tr>";
}

close OUT2;
	
print "</table>\n";
print "<p>View text version of results: <a href=\"$img_location/upload/db_res.$count.tab\">File</a></p>";
print "</body><html>\n";
