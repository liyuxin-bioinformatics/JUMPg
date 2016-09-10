#!/usr/bin/perl -wT
package PrimarySeq;
use strict;
use vars qw(@NAMES @TABLES %IUPAC_DNA $CODONS $TRCOL);

BEGIN
{
	#our (@NAMES, @TABLES, %IUPAC_DNA,$CODONS, $TRCOL);

	 @NAMES =			#id
	(
	 'Standard',		#1
	 'Vertebrate Mitochondrial',#2
	 'Yeast Mitochondrial',# 3
	 'Mold, Protozoan, and CoelenterateMitochondrial and Mycoplasma/Spiroplasma',#4
	 'Invertebrate Mitochondrial',#5
	 'Ciliate, Dasycladacean and Hexamita Nuclear',# 6
	 '', '',
	 'Echinoderm Mitochondrial',#9
	 'Euplotid Nuclear',#10
	 '"Bacterial"',# 11
	 'Alternative Yeast Nuclear',# 12
	 'Ascidian Mitochondrial',# 13
	 'Flatworm Mitochondrial',# 14
	 'Blepharisma Nuclear',# 15
	 'Chlorophycean Mitochondrial',# 16
	 '', '',  '', '',
	 'Trematode Mitochondrial',# 21
	 'Scenedesmus obliquus Mitochondrial', #22
	 'Thraustochytrium Mitochondrial' #23
	 );

     @TABLES =
	qw(
	   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
	   FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
	   FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   '' ''
	   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
	   FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
	   FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   '' '' '' ''
	   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG   
	   FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	   );

 %IUPAC_DNA = ( A => [qw(A)],
             C => [qw(C)],
             G => [qw(G)],
             T => [qw(T)],
             U => [qw(U)],
             M => [qw(A C)],
             R => [qw(A G)],
             W => [qw(A T)],
             S => [qw(C G)],
             Y => [qw(C T)],
             K => [qw(G T)],
             V => [qw(A C G)],
             H => [qw(A C T)],
             D => [qw(A G T)],
             B => [qw(C G T)],
             X => [qw(G A T C)],
             N => [qw(G A T C)]
             );

	my @nucs = qw(t c a g);
	my $x = 0;
	($CODONS, $TRCOL) = ({}, {});
	for my $i (@nucs) 
	{
		for my $j (@nucs) 
		{
			for my $k (@nucs) 
			{
				my $codon = "$i$j$k";
				$CODONS->{$codon} = $x;
				$TRCOL->{$x} = $codon;
				$x++;
			}
		}
   	 }
}


sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);

  return $self;
}


############################################################

1;

sub printFas {
	my ($self, $seqhash, $output) = @_;

	open(OUT,">$output");
	foreach my $id (keys %$seqhash) {
		print OUT ">$id\n$seqhash->{$id}\n";
	}
	close OUT;
}

sub parseFas
{
	my ($self, $fasFile) = @_;

	my $seq; # hash pointer
	my $id=''; # sequence id

	open(IN,$fasFile) || die "Cannot open $fasFile!!!";
	while(<IN>)
	{
		#chomp;
		#s/[\r]*[\n]*//gm;
		s/\r?\n//;
		if (/^>/)
		{
			s/>//;
			my @t=split /[\s\t]+/,$_;
			$id=shift @t;
			#$ss='';
			$seq->{$id}='';
		}
		else
		{
			#$ss.=$_;
			$seq->{$id}.=$_;
		}
	}
	close IN;

	return $seq;
}

sub translate
{
    my ($self, $seq, $frame) = @_;
    $frame ||= 0;
    if( $frame ) {$seq = substr ($seq,$frame); }

    my $id=1;
    my ($partial) = 0;
    $partial = 2 if length($seq) % 3 == 2;

    $seq = lc $seq;
    $seq =~ tr/u/t/;
    my $protein = "";
        if ($seq =~ /[^actg]/ )
        { #ambiguous chars
                for (my $i = 0; $i < (length($seq) - 2 ); $i+=3)
                {
                        my $triplet = substr($seq, $i, 3);
                        if (exists $CODONS->{$triplet})
                        {
                                $protein .= substr($TABLES[$id-1], $CODONS->{$triplet},1);
                        }
                        else
                        {
                                $protein .= $self->_translate_ambiguous_codon($triplet);
                        }
                }
        }
        else
        { # simple, strict translation
                for (my $i = 0; $i < (length($seq) - 2 ); $i+=3)
                {
            my $triplet = substr($seq, $i, 3);
            if (exists $CODONS->{$triplet})
                        {
                $protein .= substr($TABLES[$id-1], $CODONS->{$triplet}, 1);
                        }
                        else
                        {
                $protein .= 'X';
                        }
        }
    }

    if ($partial == 2)
        { # 2 overhanging nucleotides
                my $triplet = substr($seq, ($partial -4)). "n";
                if (exists $CODONS->{$triplet})
                {
                        my $aa = substr($TABLES[$id-1], $CODONS->{$triplet},1);
                        $protein .= $aa;
                }
                else
                {
                        $protein .= $self->_translate_ambiguous_codon($triplet, $partial);
                }
    }

    return $protein;
}


sub _translate_ambiguous_codon
{
    my ($self, $triplet, $partial) = @_;       #$self
    $partial ||= 0;
    my $id = 1;
    my $aa;
    my @codons =$self-> _unambiguous_codons($triplet);
    my %aas =();

    foreach my $codon (@codons)
    {
        $aas{substr($TABLES[$id-1],$CODONS->{$codon},1)} = 1;
    }

    my $count = scalar keys %aas;
    if ( $count == 1 )
    {
        $aa = (keys %aas)[0];
    }
    elsif ( $count == 2 )
    {
        if ($aas{'D'} and $aas{'N'}) {      $aa = 'B';  }
        elsif ($aas{'E'} and $aas{'Q'}) {           $aa = 'Z';  }
        else {      $partial ? ($aa = '') : ($aa = 'X');        }
    }
    else {      $partial ? ($aa = '') :  ($aa = 'X');    }

    return $aa;
}


sub _unambiguous_codons
{
    my ($self,$value) = @_;           #why no $self
    my @nts = ();
    my @codons = ();
    my ($i, $j, $k);

    @nts = map { $IUPAC_DNA{uc $_} }  split(//, $value);

    for my $i (@{$nts[0]})
    {
        for my $j (@{$nts[1]})
        {
            for my $k (@{$nts[2]})
            {
                push @codons, lc "$i$j$k";
            }
        }
    }

    return @codons;
}


sub translate_3frames
{
    my ($self,$seq)=@_;

    my @seqs;
    my $f = 0;
    while ($f != 3) {
        my $trl = $self->translate($seq, $f);
        push @seqs, $trl;
        $f++;
    }

    return @seqs;
}

sub translate_6frames
{
    my ($self, $seq)=@_;
    my @seqs1 = $self->translate_3frames($seq);
    $seq=$self->revcom($seq);
    my @seqs2 = $self->translate_3frames($seq);
    return @seqs1, @seqs2;
}

sub revcom
{
    my ($self, $seq)=@_;
    $seq  =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    $seq = reverse $seq;
    return $seq;
}

sub random_sequence
{
	my ($self, $length)=@_;
	my $seq=''; my @bp=qw(A T C G);
	for (my $i=0; $i<$length; $i++)
	{
		$seq.=$bp[int(rand()*4)];
	}
	return $seq;
}
=head
sub parseFas
{
	my ($self, $fas)=@_;
	my $seqs; # hash for sequences
	my ($id,$seq)=('',''); # individual record
	open(IN,"$fas") or  die "Cannot open $fas\n";
	while(<IN>)
	{
		chomp;
		if (/^>/)
		{
			if ($id ne '')
			{
				$seqs->{$id}=$seq;
				($id,$seq)=('','');
			}

			s/^>//;
			$id=(split /\s+/,$_)[0];
			$seq='';
		}
		else { $seq.=$_; }
	}
	if ($id ne '') {$seqs->{$id}=$seq;}
	close IN;

	return $seqs;
}
=cut
