#!/usr/bin/perl -w
#Written by Nathan Bachmann (nathan_bachmann@hotmail.com) under the supervision of Dr Scott Beatson
#(s.beatson@uq.edu.au) on the 13/02/12
use strict;
use CGI;
use CGI::Carp qw ( fatalsToBrowser );  
use File::Basename;
use effectorfamperlmodule;

$CGI::POST_MAX = 1024 * 10000; #limits upload to 10 Mb
use constant MAX_DIR_SIZE => 100 * 1_048_576; #total uploads set to 100 Mb

#variable list
my $line1;
my $line2;
my $line3;
my $line4;
my $line5;
my $aa;
my $seq;
my $size;
my $protein;
my $cds;
my @part;
my $model;
my $desc;
my $anno;
my %hash;
my %hash2;
my $result;
my $gene;
my $check;
my $check2 = "0";
my $check3;
my $data;
my $eValue;
my $bitScore;
my $eValueDom;
my $N;
my $hmm;
my $aliFrom;
my $aliTo;
my @field;
my @field2;
my $cood;
my $hmm2;
my $anno2;
my $aliFrom2;
my @element;
my $start;
my $end;
my $newStart;

my $header1;
my $header2;

#identify the upload directory
my $upload_dir = "/var/www/EffectorFam";
my $img_location = "http://smms-steel.biosci.uq.edu.au/EffectorFam/";

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


#print table headers
print "Content-type: text/html\n\n";

print <<Endofheader;
<html>
<head>
	<title>EffectorFAM Results</title>
</head>

<body>
<img src="$img_location/img/banner_res.gif" height="205" width="1000">
</div>

Endofheader

#check size
if (defined($filename))
{
	print "<p>EffectorFAM search complete</p>";
	#print "<p>Your upload file is $filename</p>";
}
else
{
	print "<p>Error 2: Upload file may be too big (Max upload is 10 Mb)</p>";
	exit;
}

#check to make sure the directory is under 100MB
if (dir_size($upload_dir) + $ENV{CONTENT_LENGTH} > MAX_DIR_SIZE)
{
	print "<p>Upload Directory is full</p>";
	exit;
}

#make the filename safe by removing the path using the filepaser routine from Basename
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


#prints the table headers in the browser


#check to see if the input file is nucleic acid
open (IN2, "$upload_dir/upload/$filename") or die "$!";

foreach $line1 (<IN2>)
{
	chomp $line1;
	if ($line1 =~ ">")
	{
		next;
	}
	elsif ($line1 =~ /[efilpqz]/i)
	{
		$aa = 'T';
	}
	else
	{
		$aa = 'F';
		$seq .= $line1;
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
#Runs the input sequence file against the effectorFAM HMMs
$result = `hmmscan $upload_dir/db/Efam $upload_dir/upload/moddata.$count.faa`;
open (OUT, ">$upload_dir/out/Efam.$count.out");
print OUT "$result";
close OUT;


#open the HMM list and stores HMM description in a hash
open (IN, "$upload_dir/db/hmm_list.txt");

while ($line2 = <IN>)
{
	chomp $line2;
	@part = split(/\t/, $line2);
	$model = $part[0];
	$desc = $part[2];
	$anno = $part[1];
	$hash{$model} = $desc;
	$hash2{$model} = $anno;
}

close IN;
	

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

print OUT2 "$header1\t$header2\tHMM\tAnnotation\tHMM description\tBit Score\tE Value\tNo. of Domains\t Alignment from\tAlignment to\n";

#Parser HMM results
open (IN3, "$upload_dir/out/Efam.$count.out") or die "Couldn't open Efam.$count.out\n";
while ($line3 = <IN3>)
{
	chomp $line3;
	if ($line3 =~ /^Query:\s+(\S+)/) # get protein name
	{
		$gene = $1;
		$check = "0"; #Used to ensure only best hit is record
		$check2 = "1"; #checks if HMM work
		$check3 = "0";
	}
	if ($line3 =~ /^Description:\s(.+)/) #get gene name 
	{
		$data = $1;
	}
	if ($line3 =~ /^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\w+)\s+(.+)$/) #matches main results line
	{
		$eValue = $1;
		$bitScore = $2;
		$eValueDom = $4;
		$N = $8;
		$hmm = $9;
		if (($eValue < 0.00001) and ($eValueDom < 0.00001) and ($bitScore > 100) and ($check =~ /0/))
		{
			print OUT2 "$gene\t$data\t$hmm\t$hash2{$hmm}\t$hash{$hmm}\t$bitScore\t$eValue\t$N\t";
			$check = "1";
		}
		else
		{
			next;
		}
	}
	if ($line3 =~ /^\s+\d\s!\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) #matches the alignment line
	{
		$aliFrom = $8;
		$aliTo = $9;
		if (($3 < 0.00001) and ($1 > 100) and ($check3 =~ /0/))
		{
			print OUT2 "$aliFrom\t$aliTo\n";
			$check3 = "1";
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

close IN3;	
close OUT2;


#checks the HMMer ran correctly
if ($check2 =~ /1/)
{
	$check2 = "good";
}
else
{
	print "<p> Error 01: You may have invalid characters in your input sequences.<p>";
	print "<p> Support: nathan.bachmann\@uqconnect.edu.au<p>";
	exit;
}


#Efam2gen
#
#Creates an artemis compatible tab file with effector protein marked as CDS features 
if ($aa =~ 'F')
{
	open (IN5, "$upload_dir/upload/db_res.$count.tab") or die "Couldn't open effectorFam result file number: $count\n";
	open (OUT3, ">$upload_dir/upload/Efam_annotation.$count.tab");
	while ($line5 = <IN5>)
	{
		chomp $line5;
		if ($line5 =~ /^ORF/)
		{
			next;
		}
		else
		{
			chomp $line5;
			@field2 = split(/\t/, $line5);
			$cood = $field2[1];
			$hmm2 = $field2[2];
			$anno2 = $field2[3];
			$aliFrom2 = $field2[8] * 3;
			@element = split(/\.\./, $cood);
			$start = $element[0];
			$end = $element[1];
		
			if ($end > $start)
			{
				$newStart = ($start + $aliFrom2) - 3;
				#writes the EffectorFAM results in genbank format
				print OUT3 "     misc_feature    $newStart..$end\n";
				print OUT3 "                     /note=\"$hmm2\"\n";
				print OUT3 "                     /note=\"$anno2\"\n";
				print OUT3 "                     /color=2\n";
			}
			if ($end < $start)
			{
				$newStart = ($start - $aliFrom2) + 3;
				#writes the EffectorFAM results in genbank format
				print OUT3 "     misc_feature    complement($newStart..$end)\n";
				print OUT3 "                     /note=\"$hmm2\"\n";
				print OUT3 "                     /note=\"$anno2\"\n";
				print OUT3 "                     /color=2\n";
			}
			else
			{
				next;
			}
		}
	}

	close IN5;
	close OUT3;
}

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
	</th>
EndofTable


#prints the results to the browser window
open (IN4, "$upload_dir/upload/db_res.$count.tab");

while ($line4 = <IN4>)
{
	chomp $line4;
	if ($line4 =~ /Alignment from/)
	{
		next;
	}
	else
	{
		@field = split(/\t/, $line4);
		print "<tr><td><b>$field[0]</b></td><td width=\"140\" style=\"white-space:pre-warp\">$field[1]</td><td>$field[2]</td><td width=\"140\" style=\"white-space:pre-warp\">$field[3]</td><td width=\"140\" style=\"white-space:pre-warp\">$field[4]</td><td>$field[5]</td><td>$field[6]</td><td>$field[7]</td></tr>";
	}
}

close IN4;
	
print "</table>\n";
print "<p>View text version of results: <a href=\"$img_location/upload/db_res.$count.tab\">File</a></p>";

#if upload sequence is nucleotide generate a genbank table file
if ($aa =~ 'F')
{
	print "<p>Get Artemis annotation file: <a href=\"$img_location/upload/Efam_annotation.$count.tab\">File</a></p>";
}

print "</body><html>\n";
