#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd;

# set up module, internal program paths for JUMPga
# should be run in the /programs dir
my $pipeDir = getcwd;
my @t=split /\//,$pipeDir;
my $quotedPath=join("\\/",@t);#print $quotedPath;

# Update path one by one
# main program: (A) -I $path/g (B) my $path=$pipeDir
system(qq(find $pipeDir -name 'JUMPg*.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/g/'));
system(qq(find $pipeDir -name 'JUMPg*.pl' | xargs perl -pi -e 's/my \\\$path=\\"[0-9a-zA-Z\\/\\.\\_\\-]+\\"/my \\\$path=\\"$quotedPath\\"/'));
# jump_c.pl: $code_path="$path/c/customizedDB";
system(qq(find $pipeDir -name 'jump_c.pl' | xargs perl -pi -e 's/my \\\$code_path=\\"[0-9a-zA-Z\\/\\.\\_\\-]+\\"/my \\\$code_path=\\"$quotedPath\\/c\\/customizedDB\\"/'));
# c/customizedDB/MFM_junction/junction_seq_translation.pl: -I $path/g
system(qq(find $pipeDir -name 'junction_seq_translation.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/g/'));
# trinity_tr6.pl: -I $path/c
system(qq(find $pipeDir -name 'trinity_tr6.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/c/'));
# formatFas.pl: -I $path/c
system(qq(find $pipeDir -name 'formatFas.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/c/'));
# builddb.pl: -I $path/d
system(qq(find $pipeDir -name 'builddb.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/d/'));
# DatabaseGeneration.pm
system(qq(find $pipeDir -name 'DatabaseGeneration.pm' | xargs perl -pi -e 's/my \\\$JUMPscript=\\"[0-9a-zA-Z\\/\\.\\_\\-]+\\"/my \\\$JUMPscript=\\"$quotedPath\\/s\\/jump.pl\\"/'));
# jump.pl: (A) -I $path/s (B) my $library = "$path/s"
system(qq(find $pipeDir -name 'jump.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/s/'));
system(qq(find $pipeDir -name 'jump.pl' | xargs perl -pi -e 's/my \\\$library = \\"[0-9a-zA-Z\\/\\.\\_\\-]+\\"/my \\\$library = \\"$quotedPath\\/s\\"/'));
# idsum_beta_pit.pl
system(qq(find $pipeDir -name 'idsum_beta_pit.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/f/'));
# spectrumQC.pl
system(qq(find $pipeDir -name 'spectrumQC.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/f/'));
# rundtas
system(qq(find $pipeDir -name 'rundtas.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/s/'));
system(qq(find $pipeDir -name 'rundtas.pl' | xargs perl -pi -e 's/my \\\$runsearch_shell=\\"[0-9a-zA-Z\\/\\.\\_\\-]+\\"/my \\\$runsearch_shell=\\"$quotedPath\\/rundtas\\/runsearch_shell\.pl\\"/'));
# runsearch_shell.pl
system(qq(find $pipeDir -name 'runsearch_shell.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/s/'));
# jump_g_postAnnotations.pl
system(qq(find $pipeDir -name 'jump_g_postAnnotations.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/g/'));
system(qq(find $pipeDir -name 'jump_g_postAnnotations.pl' | xargs perl -pi -e 's/my \\\$code_path=\\"[0-9a-zA-Z\\/\\.\\_\\-]+\\"/my \\\$code_path=\\"$quotedPath\\/pa\\"/'));
# pa/junction/consolidate_junction_peptides.pl
system(qq(find $pipeDir -name 'consolidate_junction_peptides.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/g/'));
# pa/ref/uniPro2ucsc.pl
system(qq(find $pipeDir -name 'uniPro2ucsc.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/g/'));
# pa/transcript_6FT/annoPep.pl: -I $path/g
system(qq(find $pipeDir -name 'annoPep.pl' | xargs perl -pi -e 's/\\-I [0-9a-zA-Z\\/\\.\\_\\-]+/\\-I $quotedPath\\/g/'));

=head
system("find $pipeDir -name '*.pm' -o -name '*.pl' -o -name '*.params' | xargs perl -pi -e 's/Release date: [0-9\\/]+/Release date: $date/'");
system("find $pipeDir -name '*.pm' -o -name '*.pl' -o -name '*.params' | xargs perl -pi -e 's/Release version: version [0-9.]+/Release version: version $versionNo/'");
system("find $pipeDir -name '*.pm' -o -name '*.pl' -o -name '*.params' | xargs perl -pi -e 's/Version [0-9.]+/Version $versionNo/g'");
system("find $pipeDir -name '*.pm' -o -name '*.pl' -o -name '*.params' | 
                xargs perl -pi -e 's/data1\\/pipeline\\/release\\/version[0-9.]+/data1\\/pipeline\\/release\\/version$versionNo/g'");
=cut
