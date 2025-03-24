use Carp::Always;
use Bio::KBase::AppService::AppScript;
use File::Slurp;
use IPC::Run;
use Cwd qw(abs_path getcwd);
use File::Path qw(rmtree make_path);
use strict;
use Data::Dumper;
use File::Basename;
use File::Temp;
use JSON::XS;
use Getopt::Long::Descriptive;
# Get perl module;
use lib '/home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_WholeGenomeSNPAnalysis/lib/';  # use NB dev

use wgSNPanalysis qw(my_function);

my $wgSNPanalysis = new WholeGenomeSNPanalysis();

my $app = Bio::KBase::AppService::AppScript->new(sub { $wgSNPanalysis->run_WholeGenomeSNPanalysis(@_); },
						 sub { $wgSNPanalysis->preflight(@_); });

$app->run(\@ARGV);