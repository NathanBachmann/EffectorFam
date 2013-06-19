#!/usr/bin/perl -wT
#Written by Nathan Bachmann (nathan_bachmann@hotmail.com) under the supervision of Dr Scott Beatson
#(s.beatson@uq.edu.au) on the 13/02/12
use strict;
use CGI;
use CGI::Carp qw ( fatalsToBrowser );

my $img_location = "http://smms-steel.biosci.uq.edu.au/EffectorFam";

print "Content-type: text/html\n\n";

print <<Endofhtml
<html>
<head>
<title>EffectorFAM</title>
</head>
<table cellspacing="0" cellpadding="0" border="0"
 bgcolor="black" id="shell" height="250" width="800">
<tr height="50"><td colspan="2" bgcolor="white">
<table title="Banner" id="banner" border="0">
<tr><td><img src="$img_location/img/banner_v2.gif" height="205" width="1000"></td></tr>
</table>
</td></tr>
<tr height="200"><td bgcolor="white" valign="top">
<table id="navigation" title="Navigation" border="0" width="200">
<tr><td><b>Current Satus:</b></tr></td>
<tr><td>Beta version 0.1</tr></td>
<tr><td></tr></td>
<tr><td><b>Example Data:</b></tr></td>
<tr><td><a href=\"$img_location/example_data/exampledata.faa\">Sequences</a></tr></td>
</table>
</td><td bgcolor="white">
<table title="Content" id="content" border="0">
<tr><td><form action="effector_upload.cgi" method="post" 
enctype="multipart/form-data">
EffectorFAM is a database of Hidden Markov Models (HMMs) designed specifically for identifying and predicting effector proteins that are translocated via the Type III Secretion System. 
<p>This portal allows users to submit protein or nucleotide sequences in FASTA format, which is then compared to the HMMs to generate a list of potential effectors.</p>
<p>If a nucleotide sequences is uploaded it will translate all ORF longer then 150 bp across all six reading frames. Note this may take a minute to calculate.</p>
<br>
<p>Input sequence(s) in fasta format:
<TEXTAREA ROWS=10 COLS=80 NAME="seq" WRAP="virtual"></TEXTAREA></p>
<p> Upload a file in fasta format: <input name="file" type="file"></p>
<p><input name="Submit" value="Find" type="submit"><input type="reset" value="Clear"></p>
<br>
<p><b>For support and inquiries:</b> nathan.bachmann\@uqconnect.edu.au</p>
</form></td></tr>
</table>
</td></tr></table>
</br>
</html>

Endofhtml
