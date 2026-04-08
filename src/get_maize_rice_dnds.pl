#!/usr/bin/env perl
use lib '/c/Users/julia/src/bioperl-1.6.924';
use lib '/c/Users/julia/src/ensembl/modules';
use lib '/c/Users/julia/ensembl-compara/modules';
use strict;
use warnings;

use Bio::EnsEMBL::Registry;

# Load Ensembl Plants registry
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'plants.ensembl.org',
    -user => 'anonymous'
);

# Get Compara adaptor
my $compara_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Multi', 'compara');

my $genome_db_adaptor = $compara_dba->get_GenomeDBAdaptor();
my $homology_adaptor  = $compara_dba->get_HomologyAdaptor();

# Get species
my $maize = $genome_db_adaptor->fetch_by_registry_name("zea_mays");
my $rice  = $genome_db_adaptor->fetch_by_registry_name("oryza_sativa");

# Fetch ALL orthologs (not just 1:1)
my $homologies = $homology_adaptor->fetch_all_by_GenomeDB_pair(
    $maize,
    $rice
);

print join("\t", qw(maize_gene rice_gene dn ds dnds type)), "\n";

foreach my $homology (@$homologies) {

    my @members = @{$homology->get_all_Members()};

    # Ensure correct species assignment
    my ($gene1, $gene2);
    if ($members[0]->genome_db->name eq 'zea_mays') {
        ($gene1, $gene2) = @members;
    } else {
        ($gene2, $gene1) = @members;
    }

    my $maize_id = $gene1->stable_id;
    my $rice_id  = $gene2->stable_id;

    my $dn   = $homology->dn;
    my $ds   = $homology->ds;
    my $dnds = $homology->dnds_ratio;
    my $type = $homology->description;

    print join("\t",
        $maize_id,
        $rice_id,
        defined $dn ? $dn : "NA",
        defined $ds ? $ds : "NA",
        defined $dnds ? $dnds : "NA",
        $type
    ), "\n";
}