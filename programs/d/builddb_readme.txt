                ##################################
                ## builddb version 11.0.0	##
                ## UPDATED: 2014-11-10          ##
                ##################################

builddb is a tool designed to build a database for JUMP or SEQUEST search and
a corresponding protein inference table (PIT) which is necessary for the
protein grouping during jump -f running.
It requires a parameter file, builddb.params, in which some parameters for the
database generation and a jump search parameter file (e.g. jump.params) should
be specified.
The followings are the descriptions of parameters

1. input_database
First, a user should specify at least one .fasta file which will be used as a
template for the database. Note that the full (absolute) path of the .fasta
file needs to be defined as following.
	input_database1 = /home/usr/MOUSE.fasta
Multiple .fasta files can also be used to generate a database as following.
	input_database1 = /home/usr/HUMAN.fasta
	input_database2 = /home/usr/MOUSE.fasta
	input_database3 = /home/usr/ECOLI.fasta
Please don't forget numbering for the input databases (i.e. .fasta files).

2. output_prefix
It is the prefix for the new database and PIT files. The suffix will be
automatically generated according to the conditions for the database (those
conditions should be defined in a jump search parameter file like
jump.params).
	output_prefix = mouse
Then, the newly generated database will have a name such as
"mouse_ft_mc2_c0.fasta.hdr (or mdx)"

3. include_contaminants
If a user wants to include contaminant sequences, this parameter needs to be
set to 1. Otherwise, set to 0.
If the user uses input_databases which already include contaminant sequences,
this parameter should be set to 0 to avoid duplicating contaminants.

4. input_contaminants
The path of a .fasta file containing contaminant sequences. Note that the full
(absolute) path should be specified.

5. decoy_generation
If a user wants to include decoy sequences of proteins (and contaminants, if
exist), this parameter needs to be set to 1. Otherwise, 0.

6. decoy_generation_method
A user can select the method of generating decoy sequences.
"decoy_generation_method = 1" will generate decoys by simply reversing protein
sequences.
"decoy_generation_method = 2" will generate decoys by reversing protein
sequences and then swapping every K and R with its preceding amino acid.

7. jump.params
A user SHOULD specify a jump search parameter file (e.g. jump.params). From
the file, the information of search engine (whether SEQUEST or JUMP),
static/dynamic modifications and so on will be obtained and be directly used
to generate a new database. Please put the full (absolute) path of the jump
search parameter file.
Although "database_name" needs to be defined in a jump.params file, it will be
ignored in the builddb. So, users don't need to care about it.

8. list_protein_abundance (optional)
This parameter is used to generate a PIT file. The protein abundance
information will be used to group and sort proteins in the PIT. Multiple
protein abundance information can be used as following.
	list_protein_abundance1 = /home/usr/mouse_abundance.txt
	list_protein_abundance2 = /home/usr/rat_abundance.txt
The file for protein abundance information should be a tab-delimited text
format with the following rules.
	Column1 = Uniprot accession of a protein (e.g. P12345)
	Column2 = abundance of the protein (numeric value)
	Column3, 4, 5, ... = any information/description (will be ignored)
Also, the full (absolute) path needs to be specified

9. list_... (optional)
These parameters are also used to generate a PIT file. A user can specify
protein annotation information such as the lists of kinases, transcription
factors and so forth. The file for a protein annotation should have the
following format (tab-delimited text)
	Column1 = Uniprot accession of a protein
	Column2, 3, 4, ...  = any information/description (will be ignored)
When a user adds any protein annotation information, the prefix "list_" should
be kept. The letters after "list_" will be used as a name of the
user-specified annotation.
	list_oncogene = /home/usr/oncogene_list.txt
Then, "oncogene" item will be added to the right end of a PIT file.
