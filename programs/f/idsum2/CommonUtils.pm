#!/usr/bin/perl 

## Release date: 11/01/2015
## Release version: version 12.1.0
## Module name: idsum2::CommonUtils

package idsum2::CommonUtils;
use strict;
use Data::Dumper;

use Sys::Hostname;
use Socket;
#Changed to Monoisotopic 1/9/05 DMD (5 decimal DMD 1/5/07)
my %residMass =
(
A => 71.0371137847,
B => 114.5962,
C => 103.0091847847,
D => 115.0269430238,
E => 129.0425930880,
F => 147.0684139130,
G => 57.0214637206,
H => 137.0589118585,
I => 113.0840639771,
K => 128.0949630140,
L => 113.0840639771,
M => 131.0404849130,
N => 114.0429274411,
O => 114.1472,
P => 97.0527638489,
Q => 128.0585775053,
R => 156.1011110236,
S => 87.0320284043,
T => 101.0476784684,
V => 99.0684139130,
W => 186.0793129499,
X => 113.1594,
Y => 163.0633285326,
Z => 128.6231
);
#(A => 71.03711, B => 114.5962, C => 103.00919, D => 115.02694, E => 129.04259, F => 147.06841,
                #G => 57.02146, H => 137.05891, I => 113.08406, K => 128.09496, L => 113.08406, M => 131.04048,
                #N => 114.04293, O => 114.1472, P => 97.05276, Q => 128.05858, R => 156.10111, S => 87.03203,
                #T => 101.04768, V => 99.06841, W => 186.07931, X => 113.1594, Y => 163.06333, Z => 128.6231);
my %phobicity =(A => 16, B => 114.5962, C => 250, D => -249, E => -150, F => 500,
                G => -331, H => -463, I => 441, K => -500, L => 476, M => 332,
                N => -379, O => 114.1472, P => -492, Q => -276, R => -277, S => -285,
                T => -108, V => 302, W => 488, X => 113.1594, Y => 200, Z => 128.6231);
my %enzymes = (trypsin => "[KR\-](?!P)", chymotrypsin => "[WYF\-](?!P)", 
               gluc_nahpo => "[ED\-](?!P)", gluc => "[E\-](?!P)", 
               lysc => "[K\-](?!P)", argc => "[R\-](?!P)", aspc => "[D]");
my $NormColor = "RED";
my $ShareColor = "BLUE";
my $H = 1.00728;

sub new {
  my($class) = @_;
  my $self = {};
  bless ($self,$class);
                                                                                                                                                             
  return $self;
}

sub server{
##################### changed by xusheng ##################
#	my $server = "proteox.genetics.emory.edu";
###### changed by yanji
        #my $server = "10.4.14.5";
	#my $server = "spiders02.stjude.org";
	#my $server = "10.4.15.169";
	#my $server = "spiderscluster.stjude.org";
	
	# Tim Modified
	my($addr)=inet_ntoa((gethostbyname(hostname))[4]);
	#print "$addr\n";
	my $server = $addr;
###### end of change
	return $server;
}

sub get_monoresidhash{
	shift ;
	my ($hash) = @_;
	%$hash = %residMass;
}

sub get_critical_values{
  shift @_;
  my ($hash) = @_;
	# R22
  my ($n, $x, $xx) = (6, 95, 99);
  $$hash{14}{$n}{$x} = .990; $$hash{14}{$n}{$xx} = .998; $n++;  $$hash{14}{$n}{$x} = .909; $$hash{14}{$n}{$xx} = .970; $n++;
  $$hash{14}{$n}{$x} = .846; $$hash{14}{$n}{$xx} = .922; $n++;  $$hash{14}{$n}{$x} = .787; $$hash{14}{$n}{$xx} = .873; $n++;
  $$hash{14}{$n}{$x} = .734; $$hash{14}{$n}{$xx} = .826; $n++;  $$hash{14}{$n}{$x} = .688; $$hash{14}{$n}{$xx} = .781; $n++;
  $$hash{14}{$n}{$x} = .648; $$hash{14}{$n}{$xx} = .740; $n++;  $$hash{14}{$n}{$x} = .616; $$hash{14}{$n}{$xx} = .705; $n++;
  $$hash{14}{$n}{$x} = .590; $$hash{14}{$n}{$xx} = .674; $n++;  $$hash{14}{$n}{$x} = .568; $$hash{14}{$n}{$xx} = .647; $n++;
  $$hash{14}{$n}{$x} = .548; $$hash{14}{$n}{$xx} = .624; $n++;  $$hash{14}{$n}{$x} = .531; $$hash{14}{$n}{$xx} = .605; $n++;
  $$hash{14}{$n}{$x} = .516; $$hash{14}{$n}{$xx} = .589; $n++;  $$hash{14}{$n}{$x} = .503; $$hash{14}{$n}{$xx} = .575; $n++;
  $$hash{14}{$n}{$x} = .491; $$hash{14}{$n}{$xx} = .562; $n++;  $$hash{14}{$n}{$x} = .480; $$hash{14}{$n}{$xx} = .551; $n++;
  $$hash{14}{$n}{$x} = .470; $$hash{14}{$n}{$xx} = .541; $n++;  $$hash{14}{$n}{$x} = .461; $$hash{14}{$n}{$xx} = .532; $n++;
  $$hash{14}{$n}{$x} = .452; $$hash{14}{$n}{$xx} = .524; $n++;  $$hash{14}{$n}{$x} = .445; $$hash{14}{$n}{$xx} = .516; $n++;
  $$hash{14}{$n}{$x} = .438; $$hash{14}{$n}{$xx} = .508; $n++;  $$hash{14}{$n}{$x} = .432; $$hash{14}{$n}{$xx} = .501; $n++;
  $$hash{14}{$n}{$x} = .426; $$hash{14}{$n}{$xx} = .495; $n++;  $$hash{14}{$n}{$x} = .419; $$hash{14}{$n}{$xx} = .489; $n++;
  $$hash{14}{$n}{$x} = .414; $$hash{14}{$n}{$xx} = .483;
	%{$$hash{15}} = %{$$hash{14}}; %{$$hash{16}} = %{$$hash{14}}; %{$$hash{17}} = %{$$hash{14}}; %{$$hash{18}} = %{$$hash{14}}; 
	%{$$hash{19}} = %{$$hash{14}}; %{$$hash{20}} = %{$$hash{14}}; %{$$hash{21}} = %{$$hash{14}}; %{$$hash{22}} = %{$$hash{14}}; 
	%{$$hash{23}} = %{$$hash{14}}; %{$$hash{24}} = %{$$hash{14}}; %{$$hash{25}} = %{$$hash{14}}; %{$$hash{26}} = %{$$hash{14}}; 
	%{$$hash{27}} = %{$$hash{14}}; %{$$hash{28}} = %{$$hash{14}}; %{$$hash{29}} = %{$$hash{14}}; %{$$hash{30}} = %{$$hash{14}}; 

	# R21
	$n = 5;
  $$hash{11}{$n}{$x} = .987; $$hash{11}{$n}{$xx} = .998; $n++;  $$hash{11}{$n}{$x} = .913; $$hash{11}{$n}{$xx} = .970; $n++;
  $$hash{11}{$n}{$x} = .828; $$hash{11}{$n}{$xx} = .919; $n++;  $$hash{11}{$n}{$x} = .763; $$hash{11}{$n}{$xx} = .868; $n++;
  $$hash{11}{$n}{$x} = .710; $$hash{11}{$n}{$xx} = .816; $n++;  $$hash{11}{$n}{$x} = .664; $$hash{11}{$n}{$xx} = .760; $n++;
  $$hash{11}{$n}{$x} = .625; $$hash{11}{$n}{$xx} = .713; $n++;  $$hash{11}{$n}{$x} = .592; $$hash{11}{$n}{$xx} = .675; $n++;
  $$hash{11}{$n}{$x} = .565; $$hash{11}{$n}{$xx} = .649; $n++;  $$hash{11}{$n}{$x} = .544; $$hash{11}{$n}{$xx} = .627; $n++;
  $$hash{11}{$n}{$x} = .525; $$hash{11}{$n}{$xx} = .607; $n++;  $$hash{11}{$n}{$x} = .509; $$hash{11}{$n}{$xx} = .580; $n++;
  $$hash{11}{$n}{$x} = .495; $$hash{11}{$n}{$xx} = .573; $n++;  $$hash{11}{$n}{$x} = .482; $$hash{11}{$n}{$xx} = .559; $n++;
  $$hash{11}{$n}{$x} = .469; $$hash{11}{$n}{$xx} = .547; $n++;  $$hash{11}{$n}{$x} = .460; $$hash{11}{$n}{$xx} = .536; $n++;
  $$hash{11}{$n}{$x} = .450; $$hash{11}{$n}{$xx} = .526; $n++;  $$hash{11}{$n}{$x} = .441; $$hash{11}{$n}{$xx} = .517; $n++;
  $$hash{11}{$n}{$x} = .434; $$hash{11}{$n}{$xx} = .509; $n++;  $$hash{11}{$n}{$x} = .427; $$hash{11}{$n}{$xx} = .501; $n++;
  $$hash{11}{$n}{$x} = .420; $$hash{11}{$n}{$xx} = .493; $n++;  $$hash{11}{$n}{$x} = .414; $$hash{11}{$n}{$xx} = .486; $n++;
	$$hash{11}{$n}{$x} = .407; $$hash{11}{$n}{$xx} = .479; $n++;  $$hash{11}{$n}{$x} = .402; $$hash{11}{$n}{$xx} = .472; $n++;  
	$$hash{11}{$n}{$x} = .396; $$hash{11}{$n}{$xx} = .466; $n++;	$$hash{11}{$n}{$x} = .391; $$hash{11}{$n}{$xx} = .460; $n++;
	%{$$hash{12}} = %{$$hash{11}}; %{$$hash{13}} = %{$$hash{11}};
	
  
	#R11
	$n = 4;                                                                                                                                                  
  $$hash{8}{$n}{$x} = .977; $$hash{8}{$n}{$xx} = .998; $n++;  $$hash{8}{$n}{$x} = .863; $$hash{8}{$n}{$xx} = .937; $n++;
  $$hash{8}{$n}{$x} = .748; $$hash{8}{$n}{$xx} = .839; $n++;  $$hash{8}{$n}{$x} = .673; $$hash{8}{$n}{$xx} = .782; $n++;
  $$hash{8}{$n}{$x} = .615; $$hash{8}{$n}{$xx} = .725; $n++;  $$hash{8}{$n}{$x} = .570; $$hash{8}{$n}{$xx} = .677; $n++;
  $$hash{8}{$n}{$x} = .534; $$hash{8}{$n}{$xx} = .639; $n++;  $$hash{8}{$n}{$x} = .505; $$hash{8}{$n}{$xx} = .606; $n++;
  $$hash{8}{$n}{$x} = .481; $$hash{8}{$n}{$xx} = .580; $n++;  $$hash{8}{$n}{$x} = .461; $$hash{8}{$n}{$xx} = .558; $n++;
  $$hash{8}{$n}{$x} = .445; $$hash{8}{$n}{$xx} = .539; $n++;  $$hash{8}{$n}{$x} = .430; $$hash{8}{$n}{$xx} = .522; $n++;
  $$hash{8}{$n}{$x} = .417; $$hash{8}{$n}{$xx} = .508; $n++;  $$hash{8}{$n}{$x} = .406; $$hash{8}{$n}{$xx} = .495; $n++;
  $$hash{8}{$n}{$x} = .396; $$hash{8}{$n}{$xx} = .484; $n++;  $$hash{8}{$n}{$x} = .386; $$hash{8}{$n}{$xx} = .473; $n++;
  $$hash{8}{$n}{$x} = .379; $$hash{8}{$n}{$xx} = .464; $n++;  $$hash{8}{$n}{$x} = .371; $$hash{8}{$n}{$xx} = .455; $n++;
  $$hash{8}{$n}{$x} = .364; $$hash{8}{$n}{$xx} = .446; $n++;  $$hash{8}{$n}{$x} = .357; $$hash{8}{$n}{$xx} = .439; $n++;
  $$hash{8}{$n}{$x} = .352; $$hash{8}{$n}{$xx} = .432; $n++;  $$hash{8}{$n}{$x} = .346; $$hash{8}{$n}{$xx} = .426; $n++;
	$$hash{8}{$n}{$x} = .341; $$hash{8}{$n}{$xx} = .420; $n++;  $$hash{8}{$n}{$x} = .337; $$hash{8}{$n}{$xx} = .414; $n++;  
	$$hash{8}{$n}{$x} = .332; $$hash{8}{$n}{$xx} = .409; $n++;	$$hash{8}{$n}{$x} = .328; $$hash{8}{$n}{$xx} = .404; $n++;
	$$hash{8}{$n}{$x} = .324; $$hash{8}{$n}{$xx} = .399; $n++;
	%{$$hash{9}} = %{$$hash{8}}; %{$$hash{10}} = %{$$hash{8}}; 
	
	#R10
	$n = 3;                                                                                                                                                  
  $$hash{3}{$n}{$x} = .970; $$hash{3}{$n}{$xx} = .994; $n++;	$$hash{3}{$n}{$x} = .829; $$hash{3}{$n}{$xx} = .926; $n++;
  $$hash{3}{$n}{$x} = .710; $$hash{3}{$n}{$xx} = .821; $n++;	$$hash{3}{$n}{$x} = .625; $$hash{3}{$n}{$xx} = .740; $n++;
  $$hash{3}{$n}{$x} = .568; $$hash{3}{$n}{$xx} = .680; $n++;	$$hash{3}{$n}{$x} = .526; $$hash{3}{$n}{$xx} = .634; $n++;
  $$hash{3}{$n}{$x} = .493; $$hash{3}{$n}{$xx} = .598; $n++;	$$hash{3}{$n}{$x} = .466; $$hash{3}{$n}{$xx} = .568; $n++;
  $$hash{3}{$n}{$x} = .444; $$hash{3}{$n}{$xx} = .542; $n++;	$$hash{3}{$n}{$x} = .426; $$hash{3}{$n}{$xx} = .522; $n++;
  $$hash{3}{$n}{$x} = .410; $$hash{3}{$n}{$xx} = .503; $n++;	$$hash{3}{$n}{$x} = .396; $$hash{3}{$n}{$xx} = .488; $n++;
  $$hash{3}{$n}{$x} = .384; $$hash{3}{$n}{$xx} = .475; $n++;	$$hash{3}{$n}{$x} = .374; $$hash{3}{$n}{$xx} = .463; $n++;
  $$hash{3}{$n}{$x} = .365; $$hash{3}{$n}{$xx} = .452; $n++;	$$hash{3}{$n}{$x} = .356; $$hash{3}{$n}{$xx} = .442; $n++;
  $$hash{3}{$n}{$x} = .349; $$hash{3}{$n}{$xx} = .433; $n++;	$$hash{3}{$n}{$x} = .342; $$hash{3}{$n}{$xx} = .425; $n++;
  $$hash{3}{$n}{$x} = .337; $$hash{3}{$n}{$xx} = .418; $n++;	$$hash{3}{$n}{$x} = .331; $$hash{3}{$n}{$xx} = .411; $n++;
  $$hash{3}{$n}{$x} = .326; $$hash{3}{$n}{$xx} = .404; $n++;	$$hash{3}{$n}{$x} = .321; $$hash{3}{$n}{$xx} = .399; $n++;
	$$hash{3}{$n}{$x} = .317; $$hash{3}{$n}{$xx} = .393; $n++;  $$hash{3}{$n}{$x} = .312; $$hash{3}{$n}{$xx} = .388; $n++;  
	$$hash{3}{$n}{$x} = .308; $$hash{3}{$n}{$xx} = .384; $n++;	$$hash{3}{$n}{$x} = .305; $$hash{3}{$n}{$xx} = .380; $n++;
	$$hash{3}{$n}{$x} = .301; $$hash{3}{$n}{$xx} = .376; $n++;	$$hash{3}{$n}{$x} = .298; $$hash{3}{$n}{$xx} = .372; $n++;
	%{$$hash{4}} = %{$$hash{3}}; %{$$hash{5}} = %{$$hash{3}}; %{$$hash{6}} = %{$$hash{3}}; %{$$hash{7}} = %{$$hash{3}};
}

sub create_dbHash  #returns a hash that contains name, annotation, and sequence for all proteins in given database
{
	shift @_;
	my ($db,$dbHash) = @_;
	open (DB, "<$db") || die "Cannot open database: $db $!\n";
                                                                                                                                                             
#  my %dbHash;
  	my $inProt = 0;
  	my ($name, $annotation, $sequence);
	my $number = 0;
  	while (<DB>)
	{
    		chomp($_);
    		if (/^>/)
		{
      			$inProt = 1;
			$_=~s/\#\#//g;
      			$_ =~ s/^>([a-zA-Z0-9\.\_\-\|\:]+)[\s\,\;]*(.*)//;
			$name = $1;
      			$annotation = $2;
      			if (!defined($annotation))
			{
				print $_;
        			$annotation = "";
      			}
########### comment out by xusheng on 11/6/2011 ########
#			print "\r     Gathering sequences from the database for: $name                    ";
      			$$dbHash{$name}{'protein'} = $name;
      			$$dbHash{$name}{'annotation'} = $annotation;
			my $species;
			$annotation =~ s/\s+\[MASS\=\d+\].*\Z//o;
			$annotation =~ s/\s+\Z//;
			if ($annotation =~ /\[[\w\.\-\s]+\]\Z/)
			{
				$annotation =~ s/\[([\w\.\-\s]+)\]\Z//o;
				$species = $1;
			} else 
			{
				$species = "Not specified";
			}
			$$dbHash{$name}{'species'} = $species; 
			$$dbHash{$name}{'number'} = $number;
			$number++;
    		} 
		elsif ($inProt)
		{
      			$sequence = $_;
	      		$sequence =~ s/\s//g;
      			$sequence =~ s/[\*\#\@\%\&\^\~\$\_]//g;
			next if ($sequence eq "");
      			$$dbHash{$name}{'sequence'} .= $sequence;
    		}
  	}
	print "\n";
  	return ($number);
}

sub get_peptides_emPAI{
  shift @_;
  my ($sequence, $min_len, $enzyme, $max_len) = @_;
  $max_len = 1000 if (!defined($max_len));
  $sequence = "-".$sequence."-";
  $enzyme = lc $enzyme;
  my $pattern = $enzymes{$enzyme} || die "$enzyme is not in Enzyme Hash\n";;
  my @parts = split(/($pattern)/, $sequence);
  my ($i, $arraysize) = (2, scalar(@parts));
  my @tryp_array;
  while ($i < $arraysize){
    if ((length($parts[$i])>=$min_len-1) && (length($parts[$i])<$max_len)){
        my $Nterm = get_Nterm(\@parts, $i);
        my $Cterm = get_Cterm(\@parts, $i);
        push (@tryp_array, "$Nterm$parts[$i]$Cterm");
    }
    $i += 2;
  }
  return \@tryp_array;
}

sub get_peptides{
	shift @_;
	my ($sequence, $min_len, $enzyme) = @_;
	$sequence = "-".$sequence."-";
	$enzyme = lc $enzyme;
	my $pattern = $enzymes{$enzyme} || die "$enzyme is not in Enzyme Hash\n";;
	my @parts = split(/($pattern)/, $sequence);
	my ($i, $arraysize) = (2, scalar(@parts));
	my @tryp_array;
	while ($i < $arraysize){
		if (length($parts[$i])>$min_len){
				my $Nterm = get_Nterm(\@parts, $i);
				my $Cterm = get_Cterm(\@parts, $i);
				push (@tryp_array, "$Nterm$parts[$i]$Cterm");
		}
		$i += 2;
	}
	return \@tryp_array;
}

sub get_Nterm{  #return 2 amino acids on N terminal: both are preceding the internal sequence
        my ($peps, $i) = @_;
        my $Nterm = "";
        my $offset = 1;
                                                                                                                                                             
        while(length($Nterm) < 3){
                if (defined($$peps[$i-$offset])){
                        $Nterm = $$peps[$i-$offset].$Nterm;
                                                                                                                                                             
                }
                last if $$peps[$i-$offset] =~ /-/;
                $offset++;
        }
				my $length = length($Nterm);
				if ($length>=3){
        $Nterm =~ s/.*(\w\w)\Z/$1\./;
				} elsif ($length==2){
					$Nterm =~ s/([\-\w]{2})/$1\./;
				} elsif ($length==1){
					$Nterm .= "\.";
				}
                                                                                                                                                             
        return $Nterm;
}
                                                                                                                                                             
sub get_Cterm{  #returns 2 amino acids on C terminal:  one w/n the internal sequence and one trailing
        my ($peps, $i) = @_;
        my $Cterm = "";
        my $offset = 1;
                                                                                                                                                             
        while(length($Cterm) < 4){
                if (defined($$peps[$i+$offset])){
                        $Cterm .= $$peps[$i+$offset];
                                                                                                                                                             
                }
								if (!defined($$peps[$i+$offset])){
									for (@$peps){
										print "test $_\n";
									}
									exit;
								}
                last if ($$peps[$i+$offset] =~ /\-/);
                $offset++;
        }
				my $length = length($Cterm);
				if ($length>=4){
					$Cterm =~ s/(\w)([\w\-]{2}).*/$1\.$2/;
				} elsif ($length==3) {
					$Cterm =~ s/(\w)(.*)/$1\.$2/;
				} elsif ($length==2) {
					$Cterm =~ s/(\w)([\w\-])/$1\.$2/;
				} else {
					$Cterm = "\.".$Cterm;
				}
                                                                                                                                                             
        return $Cterm;
}

sub get_MW{  #returns molecular weight of given sequence
	shift @_;
  my ($aacids) = @_;
  my $MW = 18.0105+$H;
	my $prev_aa = "";
  for (@$aacids){
     $MW += $residMass{$_} if defined($residMass{$_});
		 if ($_ eq "\#" || $_ eq "\*" || $_ eq "\@"){ #added DMD 4/20/06 to calculate with modifications
			 my $aa = $prev_aa.$_;
			 $MW += $residMass{$aa} if defined($residMass{$aa});
			 $MW -= $residMass{$prev_aa} if defined($residMass{$prev_aa});
		 }
		 $prev_aa = $_;
  }
  return $MW;
}

sub get_bestiso{ # Added DMD 4/21/06
	shift @_;
	my ($mw, $charge) = @_;
	
	my $M = ($mw/$charge)+$H;
	my ($first, $second);
	if ($mw > 5200){
		$first = $M+(3*$H)/$charge;
		$second = $M+(4*$H)/$charge;
	} elsif ($mw > 4500){
		$first = $M+(3*$H)/$charge;
		$second = $M+(2*$H)/$charge;
	} elsif ($mw > 3800){
		$first = $M+(2*$H)/$charge;
		$second = $M+(3*$H)/$charge;
	} elsif ($mw > 3100){
		$first = $M+(2*$H)/$charge;
		$second = $M+($H)/$charge;
	} elsif ($mw > 2400){
		$first = $M+($H)/$charge;
		$second = $M+(2*$H)/$charge;
	} elsif ($mw > 1700){
		$first = $M+($H)/$charge;
		$second = $M;
	} elsif ($mw <= 1700) {
		$first = $M;
		$second = $M+($H)/$charge;
	}
	
	return ($first, $second);
}

sub get_Pho{  #returns hydrophobicity of given sequence
	shift @_;
  my ($aacids) = @_;
  my $Pho = 0;
                                                                                                                                                           
  for (@$aacids){
    $Pho += $phobicity{$_} if defined($phobicity{$_});
  }
  return $Pho/100;
}

sub ispepseq{  #returns 1 if given item is a peptide sequence
	shift @_;
        my ($pepseq) = @_;
                                                                                                                                                             
        if ($pepseq =~ m/[a-zA-Z\-\s][\.][a-zA-Z\@\%\&\^\~\$\#\*]+[\.][a-zA-Z\-\s]/){
                return 1;
        } else {
                return 0;
        }
}

sub istryptic{  #returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
	shift @_;
        my ($Sequence) = @_;

        my $Nterm = 0;
        my $Cterm = 0;
        if ($Sequence =~ m/^[RK]\.[A-OQ-Z].*/ || $Sequence =~ m/^[\-]\..*/) {$Nterm = 1}; # Fixed bug DMD 3/20/2007
        #if (($Sequence =~ m/.\.[A-Z\#\@\*]*[KR][\#\@\*]*\.[A-OQ-Z]\Z/) || ($Sequence =~ m/.\.[A-Z\#\@\*]*\.\-\Z/)) {$Cterm = 1};
	if (($Sequence =~ m/.\.[A-Z\#\@\%\&\^\~\$\*\^\~\$]*[KR][\#\@\%\&\^\~\$\*\^\~\$]*\.[A-OQ-Z]\Z/) || ($Sequence =~ m/.\.[A-Z\#\@\%\&\^\~\$\*\^\~\$]*\.\-\Z/)) {$Cterm = 1};
        if ($Nterm && $Cterm) {  
                return 3;
        } elsif ($Nterm || $Cterm){
                return 1 if $Nterm;
		return 2 if $Cterm;
        } else {
                return 0;
        }
                                                                                                                                                             
}

sub istryptic1{  #returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
        shift @_;
        my ($Sequence) = @_;
	 $Sequence = "\." . $Sequence;
	 $Sequence .= "\.";

        my $Nterm = 1;
        my $Cterm = 1;
        if ($Sequence =~ m/\.[K|R][A-OQ-Z].*/) {$Nterm = 0}; # Fixed bug DMD 3/20/2007
        if ($Sequence =~ m/.\.[A-Z\#\@\%\&\^\~\$\*]*[KR][\#\@\%\&\^\~\$\*]*\./) {$Cterm = 0};
        if ($Nterm && $Cterm) {
                return 3;
        } elsif ($Nterm || $Cterm){
                return 1 if $Nterm;
                return 2 if $Cterm;
        } else {
                return 1;
        }
    
}

sub get_Peppos{  #returns all positions of the given peptide within the given sequence
	shift @_;
        my ($Pepseq, $Proseq) = @_;
        my $OrigPepseq = $Pepseq;
        my @Positions;
                                                                                                                                                             
        $Pepseq =~ s/\-//g;
        $Pepseq =~ s/([A-Z\s]*)\.([A-Z][A-Z]*)\.([A-Z]*)/$1$2$3/;
        my $first = $1;
        my $last = $3;
        my $fpos = index($Proseq, $Pepseq);
        while ($fpos != -1){
                $fpos += (1 + length($first));
                my $lpos = $fpos + length($Pepseq) - 1 - length($last);
                my $position = ("$fpos-$lpos");
                push (@Positions, $position);
                $fpos = index($Proseq, $Pepseq, $fpos);
        }
                                                                                                                                                             
        return \@Positions;
}

sub put_Spacer{  #places a spacer within the sequence for better display
	shift @_;
        my ($Seq) = @_;
				$Seq = uc $Seq;
        #Split Sequence into single chars
        my @Chars = split("", $Seq);
        my $Sequence = "";
                                                                                                                                                             
        #Tracking Variables
        my $charnum = 1;
        my $Norm = 1;
        my $closebrack = 1;
        my $inRed = 0;
        for (@Chars){
                if(/[A-Z0-9\s\=\/\#]/){
                        if ($Norm){
                                if ($charnum == 10){
                                        if (!$closebrack){
                                                if ($inRed){
                                                        $Sequence .= "$_<\/FONT> <FONT COLOR = \"$NormColor\">";
                                                } else {
                                                        $Sequence .= "$_<\/FONT> <FONT COLOR = \"$ShareColor\">";
                                                }
                                                $charnum = 1;
                                        } else {
                                                $Sequence .= "$_ ";
                                                $charnum = 1;
                                        }
                                } else {
                                        $Sequence .= $_;
                                        $charnum++;
                                }
                        } else {
                                $Sequence .= $_;
                                if ($_ eq "D"){
                                        $inRed = 1;
                                }
                        }
                } elsif (/</){
                        $Norm = 0;
                        $Sequence .= $_;
                } elsif (/>/){
                        $Norm = 1;
                        $Sequence .= $_;
                        if ($closebrack){
                                $closebrack = 0;
                        } else {
                                $inRed = 0;
                                $closebrack = 1;
                        }
                }
        }
                                                                                                                                                             
        return $Sequence;
}

sub put_trypColors{  #colors the found peptides: blue for shared and red for non-shared
	shift @_;
        my ($NewSeq, $Shared) = @_;
                                                                                                                                                             
        #extraction of lowercase
        my @extracts;
        my $TestSeq = $NewSeq;
        while ($TestSeq =~ s/([a-z][a-z]*)//){
                my $lower = $1;
                push (@extracts, $lower);
        }
                                                                                                                                                             
        #save original extracts for later use
        my @OrigEx = @extracts;
        #putting in colors
        for my $pep (@$Shared){
                for (@extracts){
                        $_ =~ s/$pep/\U$pep/ig;
                }
        }
        for (@extracts){
                $_ =~ s/([A-Z])/<FONT COLOR = \"$ShareColor\">$1<\/FONT>/g;
                $_ =~ s/([a-z])/<FONT COLOR = \"$NormColor\">\U$1<\/FONT>/g;
        }
                                                                                                                                                             
        #putting back into $NewSeq
        my $i = 0;
        for (@OrigEx){
                $NewSeq =~ s/$_/$extracts[$i]/g;
                $i++;
        }
        return $NewSeq;
}

sub put_Colors{  #colors the found peptides: blue for shared and red for non-shared
        shift @_;
        my ($NewSeq, $Shared) = @_;
                                                                                                                                                             
        #extraction of lowercase
        my @extracts;
        my $TestSeq = $NewSeq;
        while ($TestSeq =~ s/([a-z][a-z]*)//){
                my $lower = $1;
                push (@extracts, $lower);
        }
                                                                                                                                                             
        #save original extracts for later use
        my @OrigEx = @extracts;
        #putting in colors
        for my $pep (@$Shared){
                for (@extracts){
                        $_ =~ s/$pep/\U$pep/ig;
                }
        }
                                                                                                                                                             
        for (@extracts){
                $_ =~ s/([A-Z][A-Z]*)/<FONT COLOR = $ShareColor>$1<\/FONT>/g;
                $_ =~ s/([a-z][a-z]*)/<FONT COLOR = $NormColor>\U$1<\/FONT>/g;
        }
                                                                                                                                                             
        #putting back into $NewSeq
        my $i = 0;
        for (@OrigEx){
                $NewSeq =~ s/$_/$extracts[$i]/ig;
                $i++;
        }
                                                                                                                                                             
        return $NewSeq;
}

sub calc_stdev{
  shift @_;
  my ($sum, $sum2, $count) = @_;
                                                                                                                                                             
  my $stdev = 0;
	if ($count>1){
		if ($sum2 > $sum**2/$count){
  		$stdev = sqrt(($sum2-(($sum*$sum)/$count))/($count-1));
		}
	}
                                                                                                                                                             
  return $stdev;
}

sub calc_CV{
  shift @_;
  my ($mean, $stdev) = @_;
                                                                                                                                                             
  my $CV = 0;
  $CV = ($stdev/$mean)*100 if ($mean > 0);
                                                                                                                                                             
  return $CV;
}

sub get_line{ # added March 24, 2006 DMD
	shift @_;
	my ($x1, $y1, $x2, $y2) = @_;
	
	my $slope = ($y2-$y1)/($x2-$x1);
  my $intercept = $y1 - $slope * $x1;
	
	return ($slope, $intercept);
}

sub moving_average{
	shift @_;
  my ($array, $binsize, $sdcutoff, $cyclenum) = @_;
  $binsize = sprintf("%.0f", $binsize/2)*2;

  my $begin = ($binsize/2) - 1;
  my $end = scalar(@$array)-($binsize/2)-1;
  my $initial = $begin;
  my $i = 0;
  while ($begin <= $end){
    # Calculate initial mean and sd
    my ($sum, $sum2, $n) = (0,0,0);
    for (my $j=$i; $j<=$i+$binsize-1; $j++){
      $sum += $$array[$j]{'value'};
      $sum2 += $$array[$j]{'value'}**2;
      $n++;
    }
    my $mean = $sum/$n;
    my $stdev = calc_stdev($sum, $sum, $sum2, $n);

    # Remove outliers
    my $cycle = 1;
    while ($cycle <= $cyclenum){
      ($sum, $sum2, $n) = (0,0,0);
      for (my $j=$i; $j<=$i+$binsize-1; $j++){
        next if (abs($$array[$j]{'value'}-$mean) > $sdcutoff*$stdev);
        $sum += $$array[$j]{'value'};
        $sum2 += $$array[$j]{'value'}**2;
        $n++;
      }
      $mean = $sum/$n;
      $stdev = calc_stdev($sum, $sum, $sum2, $n);
      $cycle++;
    }
    $$array[$begin]{'norm'} = $mean;
    $$array[$begin]{'sd'} = $stdev;
    $begin++;
    $i++;
  }

  # Record values for beginning and end
  for (my $j=0; $j<$initial; $j++){
    $$array[$j]{'norm'} = $$array[$initial]{'norm'};
    $$array[$j]{'sd'} = $$array[$initial]{'sd'};
  }
  for (my $j=$end+1; $j<scalar(@$array); $j++){
    $$array[$j]{'norm'} = $$array[$end]{'norm'};
    $$array[$j]{'sd'} = $$array[$end]{'sd'};
  }
}

1;
