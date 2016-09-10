#!/usr/bin/perl

## Release date: 05/01/2015
## Release version: version 12.0.0
## Module name: Utils::Params

use strict;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

package Utils::Params;

my $ERROR;

#################################################################
#CLASS METHODS
#################################################################
sub new {
  my ($class,%args) = @_;
  my $self = {};
  $self->{PATH} = $args{-path};
  bless $self;
  
  if($args{CREATE_DEFAULT}){
    Spiders::Params->make_default($args{path});
  }
  
  eval {
    $self->parse_param();
  };
  $ERROR = $@ if($@);
  return $self;
}


#################################################################
#PUBLIC OBJECT METHODS
#################################################################
sub error {
  return $ERROR;
}

sub path {
  my $self = shift;
  $self->{PATH} = shift if(@_);
  return $self->{PATH};
}

sub enzymes {
  my $self = shift;
  $self->{ENZYMES} = shift if(@_);
  return $self->{ENZYMES};
}

sub current_enzyme {
  my($self) = @_;
  my $num = $self->get('enzyme_number');
  return $self->enzyme_by_number($num);
}

sub enzyme_by_number {
  my ($self,$num) = @_;
  while(my($key,$val) = each %{$self->{ENZYMES}}){
    return $key if($val eq $num);
  }
  return 'Unknown';
}

sub get_parameters {
	my $self=shift;

	my($quiet) = @_;

	my $runsearch_params;
	chomp(my $pwd = `pwd`);

  #if sequest.params exists in this directory, use it.
  #if not, ask for a path and copy that sequest.params into the current directory.
  #while giving the option to exit.

	if(-e './runsearch.params'){$runsearch_params = "$pwd/runsearch.params";}
	else {
    
    #if($quiet){die "Not found: sequest.params\n";}
		if($quiet){die "Not found: sequest.params\n";}
		else{
      #print "\nNot found: sequest.params\nSpecify it's location, or hit <Enter> to quit.\n";
			print "\nNot found: runsearch.params\nSpecify it's location, or hit <Enter> to quit.\n";
      
			while(1){
				print "Path to runsearch.params: ";
				chomp(my $sequest_params = <STDIN>);
				if(-e $sequest_params){
					copy($sequest_params,"./runsearch.params");
					last;
				}
				exit if($sequest_params eq "");
			}
		}
	}

	my $p = $self->new(path=>"$pwd/runsearch.params");
	return $p;
}

sub get_dynamic_modifications
{
######### dynamic modification
	my ($self,$parameter) = @_;

	my $modif;
	my $largest_modif=0;
	foreach my $param (keys %$parameter)
	{
		if($param =~/dynamic\_(\w+)/ and $parameter->{$param}>0)
		{
			$modif->{$1} = $parameter->{$param};
			if($parameter->{$param} > $largest_modif)
			{
				$largest_modif = $parameter->{$param};		
			}
		}
	}
	return ($modif,$largest_modif);
}

sub get_static_modification
{
	my ($self,$parameter) = @_;
	my $static_modif;
	foreach my $param (keys %$parameter)
	{
		if($param =~/add\_(\w+)\_/ and $parameter->{$param}>0)
		{
			$static_modif->{$1} = $parameter->{$param};
		}
	}
	return $static_modif;
}
	
sub parse_param {
  my($self) = @_;
  my $path = $self->{PATH};

  if(open(P,"< $path")){
    my($line);
    my $phash = {};
    
    while($line = <P>){
      
      my $linehash = {};
      my $comments = "";
      
      if( $line =~ s/\s*([;\#].*)$// ) {$comments = $1;}
      
      $linehash->{Comments} = $comments;
      
      if($line =~ /^(.+?)\s*=\s*(.+)$/){
        my ($key,$data) = ($1,$2);
        $data =~ s/\s+$//o;
#        $linehash->{data} = $data;
        $phash->{$key} = $data;
      }      
   }
   
   
    
    close P;
    $self->{PARAMETERS} = $phash;
	return $phash;
  }
  return 0;
  
}

sub create_sequest_params {
  my($self) = @_;
  my $phash = $self->{PARAMETERS};
  my %phash_table = %$phash;

  open FILE, ">sequest.params" or die $!;
    print FILE "[SEQUEST]\n";
    foreach (keys %phash_table) {
        my $output = "$_ = $phash_table{$_}->{data}\n";
        $output =~ s/\;//g;
        print FILE $output;
    }
  close FILE;
}

sub get_full_text {
  my($self) = @_;
  local $/ = undef;
  my $path = $self->{PATH} || return;
  open(P,"< $path") || return;
  my $p = <P>;
  close P;
  return $p;
}

sub write_text {
  my ($self) = @_;
  #figure out something safe;
}

##################################################################

sub reconstruct {
  my ($self,$newpath) = @_;
  my $path = $self->{PATH};
  open(P,"< $path") || die "Couldn't open $path: $!";
  my @lines = <P>;
  close P;
  
  my $mask = "%-40s";
  
  open(NEW,"> $newpath") || die "Couldnt't open $newpath: $!";
  my $line;
  my @newlines = ();
  my $parameters = $self->{PARAMETERS};
  
  foreach $line (@lines){
    if($line =~ /^(.+?)\s*=\s*(.+)$/ && exists $parameters->{$1}){
      printf(NEW $mask,"$1 = $parameters->{$1}->{data}");
      print NEW $parameters->{$1}->{Comments}."\n";
    }
    else{
      print NEW $line;
    }
  }
  close NEW;
  
}


sub reconstruct_new {
  my ($self,$newpath) = @_;
  my $path = $self->{PATH};
  open(P,"< $path") || die "Couldn't open $path: $!";
  my @lines = <P>;
  close P;

  my $mask = "%-40s";

  open(NEW,"> $newpath") || die "Couldnt't open $newpath: $!";
  my $line;
  my @newlines = ();
  my $parameters = $self->{PARAMETERS};

  foreach $line (@lines){
        next if($line =~/duc/);
        if($line=~/SEQUEST_ENZYME_INFO/)
        {
                @lines=();
        }
    if($line =~ /^(.+?)\s*=\s*(.+)$/ && exists $parameters->{$1}){
      printf(NEW $mask,"$1 = $parameters->{$1}->{data}");
      print NEW $parameters->{$1}->{Comments}."\n";
    }
    else{
      print NEW $line;
    }
  }
  close NEW;

}



#-----------------------------------------------------------------

sub rewrite {
  my($self) = @_;
  return $self->reconstruct($self->{PATH});
}

#################################################################

sub update {
  my($self,$key,$newval) = @_;
  my $parameters = $self->{PARAMETERS};
  
  if(exists $parameters->{$key}){
    $parameters->{$key}->{data} = $newval;
  }
}

sub get {
  my($self,$key) = @_;
  
  my $parameters = $self->{PARAMETERS};
  if($key eq 'current_enzyme'){return $self->current_enzyme();}
  if(my $p  = $parameters->{$key}){return $p->{data}}
  
  return 'Unknown';
}

#################################################################

sub set_default_dir {
  my ($self,$defaultdir) = @_;
  if(opendir(DIR,$defaultdir)){
    $defaultdir =~ s/\/$//o;
    $defaultdir .= '/';
    $self->{defaultdir} = $defaultdir;
    
    my @defaults = map {"$defaultdir".$_}
    grep {/\.params/} readdir(DIR);
    
    $self->{defaults} = \@defaults;
  }
}

sub get_defaults {
  my($self) = @_;
  return $self->{defaults};
}

sub generateSeqeust
{
	my ($self, $p,$out)=@_;

	my $e;
	my @tmp=split(/\s+/,$$p{enzyme_info});
	if ( $$p{digestion} eq 'full' ) { $e="$tmp[0] 1 1 $tmp[1] $tmp[2]"; }
	elsif ( $$p{digestion} eq 'partial' ) { $e="$tmp[0] 2 1 $tmp[1] $tmp[2]"; }
	$e =~ s/Tryptic/Trypsin/;

	my $term_diff_search_options='';
	if (defined($$p{dynamic_N_term})) { $term_diff_search_options.=$$p{dynamic_N_term}; }
	else { $term_diff_search_options.='0.000'; }
	$term_diff_search_options.=' ';
	if (defined($$p{dynamic_C_term})) { $term_diff_search_options.=$$p{dynamic_C_term}; }
	else { $term_diff_search_options.='0.000'; }

	my @ks=(keys %$p);
	my $diff_search_options='';
	foreach my $k (@ks)
	{
		if ($k =~ m/^dynamic_(\w)$/)
		{
			$diff_search_options.=" $$p{$k} $1";
		}
	}

	if($$p{peptide_tolerance_units}==1)
	{
		$$p{peptide_tolerance_units}=0;
	}
	
	open(OUT,">$out");
	print OUT <<"EOF";
[SEQUEST]
first_database_name = $$p{database_name}
peptide_mass_tolerance = $$p{peptide_tolerance}   # +/- peptide mass tolerance
peptide_mass_units = $$p{peptide_tolerance_units} # 0 = amu, 1 = mmu, 2 = ppm
ion_series = 0 0 0 $$p{ion_series}                 
fragment_ion_tolerance = $$p{frag_mass_tolerance}  # Da
num_output_lines = 10
num_results = 250
num_description_lines = 5
show_fragment_ions = 0
print_duplicate_references = 1
enzyme_info = $e
max_num_differential_per_peptide = $$p{max_modif_num}
diff_search_options =$diff_search_options
term_diff_search_options = $term_diff_search_options
mass_type_parent = 1
mass_type_fragment = 1
normalize_xcorr = 0
remove_precursor_peak = 0
ion_cutoff_percentage = 0.0000
max_num_internal_cleavage_sites = $$p{max_mis_cleavage}
protein_mass_filter = 0 0
match_peak_count = 0
match_peak_allowed_error = 1
match_peak_tolerance = 1
digest_mass_range = $$p{min_peptide_mass} $$p{max_peptide_mass}
add_Cterm_peptide = $$p{add_Cterm_peptide}
add_Cterm_protein = 0.0000
add_Nterm_peptide = $$p{add_Nterm_peptide}
add_Nterm_protein = 0.0000
add_G_Glycine = $$p{add_G_Glycine}
add_A_Alanine = $$p{add_A_Alanine}
add_S_Serine = $$p{add_S_Serine}
add_P_Proline = $$p{add_P_Proline}
add_V_Valine = $$p{add_V_Valine}
add_T_Threonine = $$p{add_T_Threonine}
add_C_Cysteine = $$p{add_C_Cysteine}
add_L_Leucine = $$p{add_L_Leucine}
add_I_Isoleucine = $$p{add_I_Isoleucine}
add_X_LorI = $$p{add_X_LorI}
add_N_Asparagine = $$p{add_N_Asparagine}
add_O_Ornithine = $$p{add_O_Ornithine}
add_B_avg_NandD = $$p{add_B_avg_NandD}
add_D_Aspartic_Acid = $$p{add_D_Aspartic_Acid}
add_Q_Glutamine = $$p{add_Q_Glutamine}
add_K_Lysine = $$p{add_K_Lysine}
add_Z_avg_QandE = $$p{add_Z_avg_QandE}
add_E_Glutamic_Acid = $$p{add_E_Glutamic_Acid}
add_M_Methionine = $$p{add_M_Methionine}
add_H_Histidine = $$p{add_H_Histidine}
add_R_Arginine = $$p{add_R_Arginine}
add_Y_Tyrosine = $$p{add_Y_Tyrosine}
add_W_Tryptophan = $$p{add_W_Tryptophan}
add_J_user_amino_acid = $$p{add_J_user_amino_acid}
add_U_user_amino_acid = $$p{add_U_user_amino_acid}
max_num_differential_AA_per_mod = $$p{max_modif_num}
EOF

close(OUT);	
}






#-----------------------------------------------------------------
sub make_default {
  my ($class,$path) = @_;
  open(OUT,"> $path") || die "Couldn't create $path";

############## changed by xusheng on 10/18/2011 ###################
# due to new format of sequest parameter file
  (my $text = <<'EOF') =~ s/^\s*//gm;
[SEQUEST]
first_database_name = /var/www/html/core/core_data/database/static/Pseudorev_MOUSE.fasta.hdr
second_database_name = 
peptide_mass_tolerance = 50
peptide_mass_units = 2
ion_series = 0 1 1 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
fragment_ion_tolerance = 0.0000
num_output_lines = 10
num_results = 250
num_description_lines = 3
show_fragment_ions = 0
print_duplicate_references = 0
enzyme_info = Trypsin 1 1 KR P
max_num_differential_per_peptide = 3
diff_search_options = 15.9949146221 M     
term_diff_search_options = 0.0000 0.0000
nucleotide_reading_frame = 0
mass_type_parent = 1
mass_type_fragment = 1
normalize_xcorr = 0
remove_precursor_peak = 0
ion_cutoff_percentage = 0.0000
max_num_internal_cleavage_sites = 2
protein_mass_filter = 0 0
match_peak_count = 0
match_peak_allowed_error = 1
match_peak_tolerance = 1.0000
partial_sequence = 
sequence_header_filter = 
digest_mass_range = 600.0000 35000.0000
add_Cterm_peptide = 0.0000
add_Cterm_protein = 0.0000
add_Nterm_peptide = 0.0000
add_Nterm_protein = 0.0000
add_G_Glycine = 0.0000
add_A_Alanine = 0.0000
add_S_Serine = 0.0000
add_P_Proline = 0.0000
add_V_Valine = 0.0000
add_T_Threonine = 0.0000
add_C_Cysteine = 57.02146374
add_L_Leucine = 0.0000
add_I_Isoleucine = 0.0000
add_X_LorI = 0.0000
add_N_Asparagine = 0.0000
add_O_Ornithine = 0.0000
add_B_avg_NandD = 0.0000
add_D_Aspartic_Acid = 0.0000
add_Q_Glutamine = 0.0000
add_K_Lysine = 0.0000
add_Z_avg_QandE = 0.0000
add_E_Glutamic_Acid = 0.0000
add_M_Methionine = 0.0000
add_H_Histidine = 0.0000
add_F_Phenylalanine = 0.0000
add_R_Arginine = 0.0000
add_Y_Tyrosine = 0.0000
add_W_Tryptophan = 0.0000
add_J_user_amino_acid = 0.0000
add_U_user_amino_acid = 0.0000
max_num_differential_AA_per_mod = 
EOF
  
  print OUT $text;
  close OUT;
}
1;
