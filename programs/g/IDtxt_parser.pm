#!/usr/bin/perl
package IDtxt_parser;
use strict;
sub new
{
        my($class) = @_;
        my $self = {};
        bless ($self,$class);
        return $self;
}

#---------------------------------------------------------------------------------------------------------
1;

sub printIDtxt {
	shift @_;
	my ($scanhash,$output)=@_;

	open(OUT,">$output");
	foreach my $outfile (keys %$scanhash) {
		print OUT $scanhash->{$outfile}->{IDtxt_line},"\n";
	}
	close OUT;
}

sub parseTab
{
	shift @_;
	my ($input,$colkey)=@_;
	my $hash;

	open(IN,$input) || die "cannot open ID.txt: $input\n";
	my $header=<IN>; chomp($header);
	my @colnames=split /\t/,$header;
	while (<IN>) {
		chomp;
		my @t=split /\t/,$_;
		for (my $i=0; $i<=$#t; $i++) {
			next if ($i==$colkey);
			$hash->{$t[$colkey]}->{$colnames[$i]}=$t[$i];
		}
	}
	close IN;

	return $hash;
}

sub parse_IDtxt
{
	#shift @_;
	my ($self,$idtxt)=@_;

	my ($scanhash,$peptidehash,$proteinhash);

	open(IN,"$idtxt") || die "cannot open ID.txt: $idtxt\n";
	while(<IN>)
	{
		chomp;
		next if (/^Database/ || /^Peptide/);
		my ($Peptide,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos,$precursor_peak_intensity_percentage)=split(/\;/,$_);
		my @dirs=split /\//,$Outfile;
		$Outfile=$dirs[$#dirs];
		my $pep=$self->intpep($Peptide);			# pure peptide sequence
		my $nomod=$self->nomod_pep($Peptide);
		#my $pep=intpep($Peptide);			# pure peptide sequence
		my $newscan=(defined($scanhash->{$Outfile}))?0:1; 	# new scan or duplicated scan for shared proteins?

		# scanhash
		$scanhash->{$Outfile}->{IDtxt_line}=$_;
		$scanhash->{$Outfile}->{measuredMH}=$measuredMH;
		$scanhash->{$Outfile}->{calcMH}=$calcMH;
		$scanhash->{$Outfile}->{ppm}=$ppm;
		$scanhash->{$Outfile}->{XCorr}=$XCorr;
		$scanhash->{$Outfile}->{dCn}=$dCn;
		$scanhash->{$Outfile}->{Ions}=$Ions;
		$scanhash->{$Outfile}->{tryptic}=$tryptic;
		$scanhash->{$Outfile}->{precursor_peak_intensity_percentage}=$precursor_peak_intensity_percentage;
		$scanhash->{$Outfile}->{peptide}=$pep;
		$scanhash->{$Outfile}->{nomod}=$nomod;
		$scanhash->{$Outfile}->{pos}=$pos;
		$scanhash->{$Outfile}->{full_peptide}=$Peptide;
		$scanhash->{$Outfile}->{proteins}->{$Protein}='';

		# peptidehash
		$peptidehash->{$pep}->{outfiles}->{$Outfile}='';
		$peptidehash->{$pep}->{full_peptide}=$Peptide;
		$peptidehash->{$pep}->{proteins}->{$Protein}='';
		if (!defined($peptidehash->{$pep}->{XCorr}) || $XCorr>$peptidehash->{$pep}->{XCorr})
		{
			$peptidehash->{$pep}->{top_score_outfile}=$Outfile;
			
			$peptidehash->{$pep}->{measuredMH}=$measuredMH;
			$peptidehash->{$pep}->{calcMH}=$calcMH;
			$peptidehash->{$pep}->{ppm}=$ppm;
			$peptidehash->{$pep}->{XCorr}=$XCorr;
			$peptidehash->{$pep}->{dCn}=$dCn;
			$peptidehash->{$pep}->{Ions}=$Ions;
			$peptidehash->{$pep}->{tryptic}=$tryptic;
			$peptidehash->{$pep}->{precursor_peak_intensity_percentage}=$precursor_peak_intensity_percentage;
		}

		# proteinhash
		# future development
	}
	close IN;

	return ($scanhash,$peptidehash);
}

sub intpep
{
	shift @_;
	my ($pep)=@_;

	if ($pep =~ m/\./)
	{
		chop($pep); chop($pep); $pep=reverse($pep);
		chop($pep); chop($pep); $pep=reverse($pep);
	}

	return $pep;
}

sub nomod_pep
{
	my ($self,$pep)=@_;
	my $nomod=$self->intpep($pep);
	$nomod =~ s/[\@\#\*\^\~\$\%]//g;
	return $nomod;
}
