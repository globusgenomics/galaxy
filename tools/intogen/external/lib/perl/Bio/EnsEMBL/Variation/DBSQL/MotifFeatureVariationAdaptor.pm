=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::MotifFeatureVariationAdaptor

=head1 SYNOPSIS
    my $reg = 'Bio::EnsEMBL::Registry';
  
    $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
    my $mfa  = $reg->get_adaptor('human', 'funcgen', 'MotifFeature');
    my $va   = $reg->get_adaptor('human', 'variation', 'Variation');
    my $vfa  = $reg->get_adaptor('human', 'variation', 'VariationFeature');
    my $mfva = $reg->get_adaptor('human','variation','MotifFeatureVariation');

    # fetch MotifFeature by dbID
    my $mf_dbID = 928174;
    my $mf      = $mfa->fetch_by_dbID($mf_dbID);

    # fetch all MotifFeatureVariations falling in the MotifFeature
    my $mfvs = $mfva->fetch_all_by_MotifFeatures([$mf]);

    # fetch by VariationFeatures
    my $variation_name = 'rs191666497';
    my $v              = $va->fetch_by_name($variation_name);
    my $vfs            = $vfa->fetch_all_by_Variation($v);
    
    # fetch all MotifFeatureVariations for the list of VariationFeatures
    $mfvs = $mfva->fetch_all_by_VariationFeatures($vfs);
     
=head1 DESCRIPTION

This adaptor allows you to fetch MotifFeatureVariation objects either by the MotifFeature
the associated VariationFeature falls in, or by VariationFeature directly. Storing
MotifFeatureVariation objects in a variation schema database is also supported. In the
database there will be a separate row for each alternative allele of a MotifFeatureVariation, 
but the methods here will fetch all alleles associated with the MotifFeatureVariation
at once.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::MotifFeatureVariationAdaptor;

use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

use base qw(Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor);

=head2 store

  Arg [1]    : Bio::EnsEMBL::Variation::MotifFeatureVariation $mfv
  Description: Store the MotifFeatureVariation in the database
  Status     : At risk

=cut

sub store {
    my ($self, $mfv, $rf) = @_;
    my $dbh = $self->dbc->db_handle;
    my $sth = $dbh->prepare_cached(q{
        INSERT DELAYED INTO motif_feature_variation (
            variation_feature_id,
            feature_stable_id,
            motif_feature_id,
            allele_string,
            somatic,
            consequence_types,
            motif_name,
            motif_start,
            motif_end,
            motif_score_delta,
            in_informative_position
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?)
    });

    for my $allele (@{ $mfv->get_all_alternate_MotifFeatureVariationAlleles }) {
        $sth->execute(
            $mfv->variation_feature->dbID,
            $rf->stable_id,
            $mfv->feature->dbID,
            $allele->allele_string,
            $mfv->variation_feature->is_somatic,
            (join ',', map { $_->SO_term } @{ $allele->get_all_OverlapConsequences }),
            $mfv->motif_feature->display_label,
            $allele->motif_start,
            $allele->motif_end,
            $allele->motif_score_delta,
            $allele->in_informative_position,
        );
    }
}

=head2 fetch_all_by_MotifFeatures

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::MotifFeature
  Description: Fetch all germline MotifFeatureVariations associated with the
               given list of MotifFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : Stable

=cut

sub fetch_all_by_MotifFeatures {
    my ($self, $motif_features) = @_;
    my $regulatory_features = $self->_associated_regulatory_features($motif_features);
    return $self->fetch_all_by_MotifFeatures_with_constraint($regulatory_features, $motif_features, 'somatic = 0');
}

=head2 fetch_all_somatic_by_MotifFeatures

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::MotifFeature
  Description: Fetch all somatic MotifFeatureVariations associated with the
               given list of MotifFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : Stable

=cut

sub fetch_all_somatic_by_MotifFeatures {
    my ($self, $motif_features) = @_;
    my $regulatory_features = $self->_associated_regulatory_features($motif_features);
    return $self->fetch_all_by_MotifFeatures_with_constraint($regulatory_features, $motif_features, 'somatic = 1');
}

=head2 fetch_all_by_MotifFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::MotifFeature
  Arg [2]    : listref of SO terms
  Description: Fetch all germline MotifFeatureVariations associated with the
               given list of MotifFeatures and with consequences from given list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : At risk

=cut

sub fetch_all_by_MotifFeatures_SO_terms {
    my ($self, $motif_features, $terms) = @_;
    my $regulatory_features = $self->_associated_regulatory_features($motif_features);
    my $constraint = $self->_get_consequence_constraint($terms);
    return $self->fetch_all_by_MotifFeatures_with_constraint($regulatory_features, $motif_features, $constraint.' AND somatic = 0');
}

=head2 fetch_all_somatic_by_MotifFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::MotifFeature
  Arg [2]    : listref of SO terms
  Description: Fetch all somatic MotifFeatureVariations associated with the
               given list of MotifFeatures and with consequences from given list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : At risk

=cut

sub fetch_all_somatic_by_MotifFeatures_SO_terms {
    my ($self, $motif_features, $terms) = @_;
    my $regulatory_features = $self->_associated_regulatory_features($motif_features);
    my $constraint = $self->_get_consequence_constraint($terms);
    return $self->fetch_all_by_MotifFeatures_with_constraint($regulatory_features, $motif_features, $constraint.' AND somatic = 1');
}

sub _associated_regulatory_features {
    my ($self, $motif_features) = @_;
    my $rfa = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
                    -species => $self->db->species, 
                    -type    => 'RegulatoryFeature',
                );			
    my @regulatory_features;
    for my $mf (@$motif_features) {
        my $rf = $rfa->fetch_all_by_attribute_feature($mf)->[0];
        push @regulatory_features, $rf;
    }
    return \@regulatory_features;
}

=head2 fetch_all_by_VariationFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeature
  Arg [2]    : listref of SO terms
  Description: Fetch all germline MotifFeatureVariations associated with the
               given list of VariationFeatures and with consequences from given list of SO terms
  Returntype : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : At risk

=cut

sub fetch_all_by_VariationFeatures_SO_terms {
    my ($self, $vfs, $motif_features, $terms, $without_children, $included_so) = @_;
    my $constraint = $self->_get_consequence_constraint($terms, $without_children, $included_so);
    if (!$constraint) {
      return [];
    }
    return $self->SUPER::fetch_all_by_VariationFeatures_with_constraint($vfs, $motif_features, $constraint);
}

=head2 count_all_by_VariationFeatures_SO_terms

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Arg [2]    : listref of SO terms
  Description: Count MotifFeatureVariations associated with given
               VariationFeatures and with consequences from given list of SO terms
  Returntype : int
  Status     : At risk

=cut

sub count_all_by_VariationFeatures_SO_terms {
    my ($self, $vfs, $motif_features, $terms, $included_so) = @_;
    my $constraint = $self->_get_consequence_constraint($terms, 1, $included_so);
    if (!$constraint) {
      return 0;
    }
    return $self->SUPER::count_all_by_VariationFeatures_with_constraint($vfs, $motif_features, $constraint);
}

=head2 fetch_all_by_MotifFeatures_with_constraint

  Arg [1]    : listref of Bio::EnsEMBL::Funcgen::MotifFeature
  Arg [2]    : extra SQL constraint for the query
  Description: Fetch all MotifFeatureVariations associated with the
               given list of MotifFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : At risk

=cut

sub fetch_all_by_MotifFeatures_with_constraint {
    my ($self, $regulatory_features, $motif_features, $constraint) = @_;
    my $mf_ids = join ',', map {$_->dbID} @$motif_features;
    my $features = $self->SUPER::fetch_all_by_Features_with_constraint($regulatory_features, $constraint . " AND motif_feature_id IN ($mf_ids)");
    for my $feature (@$features) {
        $feature->{feature} = undef;
        $feature->{feature} = $feature->motif_feature;
    }   
    return $features;
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    #warn $sth->sql;
    my (
        $motif_feature_variation_id,
        $variation_feature_id, 
        $feature_stable_id,
        $motif_feature_id, 
        $allele_string,
        $somatic,
        $consequence_types,
        $motif_name,
        $motif_start,
        $motif_end,
        $motif_score_delta,
        $in_informative_position,
    );
    
    $sth->bind_columns(
        \$motif_feature_variation_id,
        \$variation_feature_id, 
        \$feature_stable_id,
        \$motif_feature_id,
        \$allele_string,
        \$somatic,
        \$consequence_types,
        \$motif_name,
        \$motif_start,
        \$motif_end,
        \$motif_score_delta,
        \$in_informative_position,
    );
    
    my %mfvs;

    while ($sth->fetch) {
        my ($ref_allele, $alt_allele)   = split /\//, $allele_string;
        my $key = $variation_feature_id.'_'.$feature_stable_id.'_'.$motif_feature_id;
        my $mfv = $mfvs{$key};
      
        unless ($mfv) {
            $mfv = Bio::EnsEMBL::Variation::MotifFeatureVariation->new_fast({
                _variation_feature_id        => $variation_feature_id,
                _feature_stable_id           => $feature_stable_id,
                regulatory_feature_stable_id => $feature_stable_id,
                motif_feature_id             => $motif_feature_id,
                motif_name                   => $motif_name,
                adaptor                      => $self,
            });
            
            $mfvs{$key} = $mfv;
            
            my $ref_allele = Bio::EnsEMBL::Variation::MotifFeatureVariationAllele->new_fast({
                is_reference            => 1,
                variation_feature_seq   => $ref_allele,
                motif_feature_variation => $mfv,
                motif_start             => $motif_start,
                motif_end               => $motif_end,
                dbID                    => $motif_feature_variation_id,
            });

            $mfv->add_MotifFeatureVariationAllele($ref_allele);
        }

        my $overlap_consequences = [ map { $OVERLAP_CONSEQUENCES{$_} } split /,/, $consequence_types ];
        my $allele = Bio::EnsEMBL::Variation::MotifFeatureVariationAllele->new_fast({
            is_reference                => 0,
            variation_feature_seq       => $alt_allele,
            motif_feature_variation     => $mfv, 
            overlap_consequences        => $overlap_consequences,
            motif_start                 => $motif_start,
            motif_end                   => $motif_end,
            motif_score_delta           => $motif_score_delta,
            in_informative_position     => $in_informative_position,
            dbID                        => $motif_feature_variation_id,
        });
        $mfv->add_MotifFeatureVariationAllele($allele);
    }
    return [values %mfvs];
}

sub _tables {
    return (
        ['motif_feature_variation', 'mfv']
    );
}

sub _columns {
    return 
    qw(
        motif_feature_variation_id 
        variation_feature_id 
        feature_stable_id
        motif_feature_id 
        allele_string
        somatic 
        consequence_types 
        motif_name
        motif_start
        motif_end
        motif_score_delta
        in_informative_position    
    );
}

1;
