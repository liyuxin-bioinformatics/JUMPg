#!/usr/bin/perl

use strict;
use warnings;

if (scalar(@ARGV)!=3)
{
        die "perl annotateIntronic.pl reads.alg.out.updated reads_vs_geneLocus.tab pep2readPos.txt\n";
}


