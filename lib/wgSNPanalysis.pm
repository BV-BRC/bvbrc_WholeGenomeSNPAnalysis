package wgSNPanalysis;

use File::Slurp;
use IPC::Run;
use Cwd qw(abs_path getcwd);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Path qw(rmtree make_path);
use File::Spec;
use strict;
use Data::Dumper;
use File::Basename;
use File::Temp;
use JSON;
use Text::CSV qw(csv);
use Getopt::Long::Descriptive;
use P3DataAPI;
use YAML;
use Template;
use Module::Metadata;
use warnings;
use strict;

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(work_dir staging_dir output_dir app params));

sub new
{
    my($class) = @_;
    my $self = {};
    bless $self, $class;

    $self->{api} = P3DataAPI->new();
    
    return $self;
}

sub run_WholeGenomeSNPanalysis
{

    my($self, $app, $app_def, $raw_params, $params) = @_;
    my $begin_time = time();
    print STDERR "Processed parameters for application " . $app->app_definition->{id} . ": ", Dumper($params);
    #
    # Set up work and staging directories.
    # 
    $self->app($app);
    $self->params($params);

    my $tmp_dir = $params->{_tmpdir};
    $tmp_dir //= getcwd . "/tmp.$$";
    
    my $work_dir = "$tmp_dir/work";
    my $staging_dir = "$tmp_dir/staging";
    my $raw_fasta_dir = "$staging_dir/raw_fastas";
    my $clean_fasta_dir = "$tmp_dir/clean_fastas";
    my $output_dir = "$tmp_dir/output";


    $self->{work_dir} = $work_dir;
    $self->{staging_dir} = $staging_dir;
    $self->{raw_fasta_dir} = $raw_fasta_dir;
    $self->{clean_fasta_dir} = $clean_fasta_dir;
    $self->{output_dir} = $output_dir;
    

    eval {
        make_path($work_dir, $staging_dir, $raw_fasta_dir, $clean_fasta_dir, $output_dir);
        run($app, $params, $staging_dir, $output_dir, $work_dir, $raw_fasta_dir, $clean_fasta_dir);
        };
    my $err = $@;
    if ($err)
    {
	warn "Run failed: $@";
    }
    save_output_files($app, $output_dir);
    my $end_time = time();
    printf STDERR "TOTAL TIME ELAPSED %.2f\n", $end_time - $begin_time ;
}


sub run
{
    my($app, $params, $staging_dir, $output_dir, $work_dir, $raw_fasta_dir, $clean_fasta_dir) = @_;

    #
    # set up for genome groups
    #
    if ($params->{input_genome_type} ne 'genome_group' and
	$params->{input_genome_type} ne 'genome_fasta')
    {
	    die "Gene set type $params->{input_genome_type} not supported";
    }

    # FUTURE DEV NOTE: If we support other input options
    # The fastas need to end up in the raw_fasta_dir - NB March 2025
    if ($params->{input_genome_type} eq 'genome_group')
    {
        my $genome_group_path = $params->{input_genome_group};
        my $group_name = basename($genome_group_path);
        my $api = P3DataAPI->new;
        # my $group_genome_ids = $api->retrieve_patric_ids_from_genome_group($genome_group_path);
        my @group_genome_ids = $api->retrieve_patric_ids_from_genome_group($genome_group_path);
        # $api ->retrieve_contigs_in_genomes($group_genome_ids, $raw_fasta_dir, "%s");
        $api ->retrieve_contigs_in_genomes(@group_genome_ids, $raw_fasta_dir, "%s");
        print Dumper @group_genome_ids;
        # my @genome_metadata_fields = (
        #   "genome_name", "genome_id", "ncbi_taxon_id", "organism_name", "taxon_lineage_ids", "taxon_lineage_names",
        #   "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "genome_status", "strain", 
        #   "serovar", "isolation_source", "isolation_comments", "collection_date", "collection_year", "season", 
        #   "isolation_country", "geographic_group", "geographic_country", "geographic_location", "host_name", "host_common_name",
        #   "host_gender", "host_age", "lab_host", "passage", "sequencing_center", "sequencing_status", "sequencing_platform", "chromosomes", 
        #   "plasmids", "contigs", "genome_length", "gc_content", "contig_l50", "contig_n50", "cds", "genome_quality", 
        #   "antimicrobial_resistance", "antimicrobial_resistance_evidence", "bioproject_accession", "biosample_accession", "assembly_accession",
        #   "sra_accession", "genbank_accession", "refseq_accession", "comments", "date_inserted", "date_modified", "public", "owner", "members", "",
        #   "accession", "patric_id", "public", "genome_status", "state");
        my @genome_metadata_fields = (
        "genome_id", "genome_name", "strain", "genbank_accessions", "subtype", "lineage", "clade", "h1_clade_global", "h1_clade_us", "h5_clade", "host_group", 
        "host_common_name", "host_scientific_name", "collection_year", "geographic_group", "isolation_country", "genome_status", "state_province", "state");
        my @genome_group_metadata = $api->retrieve_genome_metadata(@group_genome_ids, \@genome_metadata_fields);
        my $json_string = encode_json(@genome_group_metadata);
        print Dumper @genome_group_metadata;
        # write metdata to json
        my $top = getcwd;
        write_file("$top/genome_metadta.json", JSON::XS->new->pretty->canonical->encode(\@genome_group_metadata));
    }

    
    # Prep the config and get ready to run Snakemake
    print STDERR "Starting the config json....\n";
    my $json_string = encode_json($params);

    #
    # Create json config file for the execution of this workflow.
    # If we are in a production deployment, we can find the workflows
    # by looking in $KB_TOP/workflows/app-name
    # Otherwise they are in the module directory; this is indicated
    # by the value of $KB_MODULE_DIR (note this is set for both
    # deployed and dev-container builds; the deployment case
    # is determined by the existence of $KB_TOP/workflows)
    #

    my %config_vars;
    # temp

    my $wf_dir = "$ENV{KB_TOP}/workflows/$ENV{KB_MODULE_DIR}";
    if (! -d $wf_dir)
    {
	$wf_dir = "$ENV{KB_TOP}/modules/$ENV{KB_MODULE_DIR}/workflow";
    }
    -d $wf_dir or die "Workflow directory $wf_dir does not exist";

    #
    # Find snakemake. We need to put this in a standard location in the runtime but for now
    # use this.
    #
    my @sm_dirs = ("$ENV{KB_RUNTIME}/bin", "$ENV{KB_RUNTIME}/artic-ncov2019/bin");
    my $snakemake;
    for my $dir (@sm_dirs)
    {
	my $p = "$dir/snakemake";
	if (-x $p)
	{
	    $snakemake = $p;
	    last;
	}
    }
    if (!$snakemake)
    {
	die "Cannot find snakemake in @sm_dirs\n";
    }
    warn "Snakemake found at $snakemake\n";
    system($snakemake, "--version");

    $config_vars{cores} = $ENV{P3_ALLOCATED_CPU} // 2;
    $config_vars{snakemake} = $snakemake;
    $config_vars{workflow_dir} = $wf_dir;
    $config_vars{input_data_dir} = $staging_dir;
    $config_vars{output_data_dir} = $output_dir;
    $config_vars{work_data_dir} = $work_dir;
    $config_vars{clean_data_dir} = $clean_fasta_dir;
    $config_vars{raw_fasta_dir} = $raw_fasta_dir;

    # add the params to the config file
    $config_vars{params} = $params;

    # write config to current working directory to avoid finding tmp dir
    my $top = getcwd;
    write_file("$top/config.json", JSON::XS->new->pretty->canonical->encode(\%config_vars));

    print STDERR "Check point 1: Starting snakemake....\n";

    my @cmd = (
        $snakemake,
        "--cores",
        $config_vars{cores},
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going",
        "--snakefile",
        "$wf_dir/snakefile/prep_ksnp_snakefile"
    );
    print STDERR "Run: @cmd\n";

    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "Snakemake prep_ksnp_snakefile command failed $?: @cmd";
    }

        my @cmd2 = (
        $snakemake,
        "--cores",
        $config_vars{cores},
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going",
        "--snakefile",
        "$wf_dir/snakefile/run_ksnp_snakefile"
    );
    print STDERR "Run: @cmd2\n";

    my $ok2 = IPC::Run::run(\@cmd2);
    if (!$ok2)
    {
     die "Snakemake prep_ksnp_snakefile command failed $?: @cmd2";
    }
}

#
# Run preflight to estimate size and duration.
#
sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $time = 60 * 60 * 12;
    my $pf = {
            	cpu => 1,
            	memory => "32G",
            	runtime => $time,
                storage => 0,
              };

    return $pf;
}

sub save_output_files
{
    my($app, $output) = @_;

    my %suffix_map = (
            phyloxml => 'phyloxml',
            tre => 'nwk',
            NJ => 'txt',
            ML => 'txt',
            SNPs_all => 'tsv',
            SNPs_all_matrix => 'txt',
            core_SNPs => 'tsv',
            5 => 'tsv',
            nonCore_SNPs => 'tsv',
            core_SNPs => 'tsv',
            core_SNPs_matrix => 'txt',
            vcf => 'vcf',
            fasta => 'contigs',
            html => 'html');
    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

    if (opendir(D, $output))
    {
	while (my $p = readdir(D))
	{
	    next if ($p =~ /^\./);
	    my @cmd = ("p3-cp", "--recursive", @suffix_map, "$output/$p", "ws:" . $app->result_folder);
	    print STDERR "saving files to workspace... @cmd\n";
	    my $ok = IPC::Run::run(\@cmd);
	    if (!$ok)
	    {
		warn "Error $? copying output with @cmd\n";
	    }
	}
    closedir(D);
    }
}

1;