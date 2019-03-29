EffectorFam
===========

EffectorFAM is a database of Hidden Markov Models (HMMs) designed specifically
for identifying and predicting effector proteins that are translocated via the
Type III Secretion System.

This web based tool allows users to submit protein or nucleotide sequences in
FASTA format, which is then compared to the HMMs to generate a list of 
potential effectors.

If given nucleotide sequences EffectorFam will translate all ORFs longer then
150 bp across all six reading frames. Note this may take a minute to calculate.

Written by Nathan Bachmann (nathan_bachmann@hotmail.com) under the supervision       
of Dr Scott Beatson (s.beatson@uq.edu.au)


Basic install instructions (on a LAMP Ubuntu 12.04.2 LTS box)::

    $ cd /var/www
    $ sudo git clone https://github.com/NathanBachmann/EffectorFam.git
    $ cd EffectorFam
    $ mkdir cgi-bin
    $ vi /etc/apache2/sites-enabled/000-default
    $ #Updating so looks like this:
         #ScriptAlias /cgi-bin/ /var/www/EffectorFam/cgi-bin/
         #<Directory "/var/www/EffectorFam/cgi-bin">
    $ vi effector_upload.cgi
    $ #Updating the variables on lines 55 and 56 respectively (to current site)
    $ sudo apt-get install hmmer
    $ sudo service apache2 restart
    $ Browse over to: To be updated

Official live site
------------------

Please see: To be updated
