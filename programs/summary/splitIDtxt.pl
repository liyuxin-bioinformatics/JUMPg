#!/usr/bin/env perl

use strict;
use warnings;

if (scalar(@ARGV)!=5)
{
        die "Usage: perl splitIDtxt.pl idsum_dir mut_ids jun_ids rdt_ids output_dir\n";
}

my $sum_dir=$ARGV[0];
my $mut_ids=$ARGV[1];
my $jun_ids=$ARGV[2];
my $rdt_ids=$ARGV[3];
my $output_dir=$ARGV[4];

my (%id2type,%types,%ids,%lines,%typeCounts);
# %id2type{cu_1234}=mut|jun|6FT
unless ($mut_ids eq 'na') {assignType(\%id2type,$mut_ids,'mutation');}
unless ($jun_ids eq 'na') {assignType(\%id2type,$jun_ids,'junction');}
#unless ($rdt_ids eq 'na') {assignType(\%id2type,$rdt_ids,'6FT');}
unless ($rdt_ids eq 'na') {assignType(\%id2type,$rdt_ids,'transcript_6FT');}
# %types{$peptide}=ref|cu_1234
# %typeCounts{ref_proteins|mut|jun|6FT}{peptide|PSM}
typeHash("$sum_dir/publications/id_uni_pep.txt",\%types,\%id2type,\%typeCounts);
# only for non-ref IDs:  %ids{cu_1234}='GSTO1|uc001kya.3|p.A140D|chr10:106022789-106022789.C_A'
idHash("$sum_dir/publications/id_all_pep.txt",\%ids);
# %lines{$type}{$line}=''
splitIDtxt("$sum_dir/ID.txt",\%types,\%ids,\%lines);
printSplitIDtxt(\%lines,$output_dir,\%typeCounts);
#----------------------------------
sub assignType
{
	my ($id2type,$input,$type)=@_;
	open(IN,$input) or die "cannot open ID.txt: $input\n";
	while(<IN>)
	{
		chomp;
		my $id=(split /\t/,$_)[0];
		$$id2type{$id}=$type;
	}
	close IN;
}

sub printSplitIDtxt
{
	my ($lines,$output_dir,$typeCounts)=@_;

	unless (-e $output_dir) { system(qq(mkdir $output_dir)); }
	system(qq(cp $sum_dir/publications/id_uni_pep.txt $output_dir/accepted_peptides.txt));

	foreach my $type (keys %{$lines})
	{
		my $output="$output_dir/$type\_ID.txt";
		#$output =~ s/\|/_/;
		open(OUT,">$output");
		foreach my $line (keys %{$$lines{$type}})
		{
			print OUT "$line\n";
		}
		close OUT;
	}

	open(OUT,">$output_dir/JUMPg_summary.txt");
	print OUT "Type\tPSMs\tPeptides\n";
	#print  "Type\tPSMs\tPeptides\n";
	foreach my $type (keys %{$typeCounts})
	{
		print OUT "$type\t$$typeCounts{$type}{PSM}\t$$typeCounts{$type}{peptide}\n";
		#print "$type\t$$typeCounts{$type}{PSM}\t$$typeCounts{$type}{peptide}\n";
	}
	close OUT;
}

sub splitIDtxt
{
	my ($input,$types,$ids,$lines)=@_;

	open(IN,$input) or die "cannot open ID.txt: $input\n";
	#open(REF,">ref_proteins_ID.txt");
	while(<IN>)
	{
		next if (/^Database/);
		next if (/^Peptide/);
		chomp;
		my ($peptide,$accession,$Outfile)=split /\;/,$_;
		$peptide = (split /\./,$peptide)[1];
		next if ($accession =~ /Decoy/);
		if (!defined($$types{$peptide}))
		{
			print "WARNING: $peptide in ID.txt not found in id_uni_pep.txt (skipped row)\n";
			next;
		}

		#if ($$types{$peptide} eq 'ref_proteins' and $accession =~ m/^sp/) {print REF "$_\n";}
		if ($$types{$peptide} eq 'ref_proteins' and $accession !~ m/^cu/) # ref peptide, ref accession
		{
			$$lines{$$types{$peptide}}{$_}='';
		}
		#elsif ($accession =~ m/$$types{$peptide}/)
		else
		{
			if (!defined($$ids{$accession}))
			{
				print "WARNING: $accession ($peptide) in ID.txt not found in id_all_pep.txt (skipped row)\n"; 
				next;
			}
			my @t=split /\;/,$_;
			my $accessionType=$id2type{$accession};
			$t[1]=$$ids{$accession}; # non-ref peptide, get full ID name
			# make sure the accession type and peptide type are matched
			my $peptideType=$$types{$peptide};
			next unless ($accessionType eq $peptideType);
			$$lines{$peptideType}{join(";",@t)}='';
		}
		#else {next;}
	}
	close IN;
}

sub idHash
{
	my ($input,$ids)=@_;

	open(IN,$input) or die "cannot open id_uni_pep.txt: $input\n";
	my $line=<IN>;
	$line=<IN>;
	$line=<IN>;
	$line=<IN>;
	while(<IN>)
	{
		chomp;
		my ($peptide,$group,$accession,$annotation,$gn)=split /\t/,$_;
		$peptide = (split /\./,$peptide)[1];

		my ($left,$central,$right)=split /\|/,$accession;
		if ($left eq 'cu')
		{
			$$ids{$accession}=(split /\s+/,$annotation)[0];
		}
	}
	close IN;
}

sub typeHash
{
	# %types{$peptide}=ref_proteins|cu_1234
	my ($input,$types,$id2type,$typeCounts)=@_;

	open(IN,$input) or die "cannot open id_uni_pep.txt: $input\n";
	my $line=<IN>;
	$line=<IN>;
	$line=<IN>;
	$line=<IN>;
	while(<IN>)
	{
		chomp; my $type='';
		my ($peptide,$group,$accession,$annotation,$gn,$PSMs,$run)=split /\t/,$_;
		$peptide = (split /\./,$peptide)[1];

		#my ($left,$central,$right)=split /\|/,$accession;
		#if ($left eq 'cu')
		if (defined($$id2type{$accession}))
		{
			#my ($type,$id)=split /\_/,$central;
			$$types{$peptide}=$$id2type{$accession};
			$type=$$id2type{$accession};
		}
		else { $$types{$peptide}='ref_proteins'; $type='ref_proteins'; }

		if (!defined($$typeCounts{$type}{'peptide'}))
		{
			$$typeCounts{$type}{'peptide'}=$$typeCounts{$type}{'PSM'}=0;
		}
		$$typeCounts{$type}{'peptide'}++;
		$$typeCounts{$type}{'PSM'}+=$PSMs;
	}
	close IN;
}
