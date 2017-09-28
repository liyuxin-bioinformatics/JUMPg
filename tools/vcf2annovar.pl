#!/usr/bin/perl

my $vcf=shift @ARGV;

open(IN,$vcf) || die "Cannot open vcf file $vcf\n";
while(<IN>)
{
	chomp;
	s/^\s*//;
	next if (/^#/);

	my ($chr,$pos,$id,$ref,$alts,$score,$filter,$infor)=split /\t/,$_;
	next unless ($filter eq 'PASS'); # only use 'filtered' mutations

	my @altss=split /\,/,$alts;
	foreach my $alt (@altss)
	{
		my ($start,$end);
		# SNV
		if ( length($ref)==1 and length($alt)==1 ){
			$start=$end=$pos;
		}
		# insertion
		# vcf: chr1    20809349        .       C       CGTGT
		# ann: chr1    20809350        20809350       -		GTGT
		elsif ( length($ref)==1 and length($alt)>1 ){
			$start=$end=$pos+1;
			$ref='-';
			$alt=shiftHeadLetter($alt);
		}
		# deletion
		# vcf: chr17   17039561        .       CCAG    C
		# ann: chr17   17039562        17039564        CAG     -
		elsif ( length($ref)>1 and length($alt)==1 ){
			$start=$pos+1;
			$end=$pos+length($ref)-1;
			$alt='-';
			$ref=shiftHeadLetter($ref); # 
		}
		# other special case? 
		# chr1    168170853       .       GACAC   G,GACACAC 
		# micro satellites: less chance to be translated; ignore in this version
		else { next; }

		# print
		print "$chr\t$start\t$end\t$ref\t$alt\n";
	}
}
close IN;

#--------------------------------------
sub shiftHeadLetter
{
	my ($s)=@_;
	$s=reverse($s);
	chop($s);
	return reverse($s);
}
