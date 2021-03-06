# $Id: Path.pm,v 1.1.2.2 2003/03/27 10:07:56 lapp Exp $
#
# BioPerl module for Path
#
# Cared for by Hilmar Lapp <hlapp at gmx.net> 
#
# (c) Hilmar Lapp, hlapp at gmx.net, 2003.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2003.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Path - a path for an ontology term graph

=head1 SYNOPSIS

  $path = Bio::Ontology::Path->new( -identifier     => "16847",
                                    -subject_term   => $subj,
                                    -object_term    => $obj,
                                    -predicate_term => $pred,
                                    -distance       => 3 );

=head1 DESCRIPTION

This is a basic implementation of Bio::Ontology::PathI.

Essiantially this is a very thin extension of
L<Bio::Ontology::Relationship>. It basically adds a method distance().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

 Hilmar Lapp <hlapp@gmx.net>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::Path;
use vars qw( @ISA );
use strict;
use Bio::Ontology::PathI;
use Bio::Ontology::Relationship;

@ISA = qw( Bio::Ontology::Relationship
           Bio::Ontology::PathI );




=head2 new

 Title   : new
 Usage   : $rel = Bio::Ontology::Path->new(-identifier   => "16847",
                                           -subject_term => $subject,
                                           -object_term  => $object,
                                           -predicate_term => $type );
                                           -distance     => 3 );
 Function: Creates a new Bio::Ontology::Path.
 Returns : A new Bio::Ontology::Path object.
 Args    : -identifier     => the identifier of this relationship [scalar]
           -subject_term   => the subject term [Bio::Ontology::TermI]
           -object_term    => the object term [Bio::Ontology::TermI]  
           -predicate_term => the predicate term [Bio::Ontology::TermI]
           -distance       => the distance between subject and object

=cut

sub new {

    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );
   
    my ( $distance ) = 
	$self->_rearrange( [qw( DISTANCE)
			    ], @args );
   
    $distance      && $self->distance($distance);
                                                    
    return $self;
    
} # new



=head2 init

 Title   : init()
 Usage   : $rel->init();   
 Function: Initializes this Path to all undef.
 Returns : 
 Args    :

=cut

sub init {
    my $self = shift;
    
    $self->SUPER::init(@_);
    $self->{ "_distance" } = undef;
   
} # init


=head2 distance

 Title   : distance
 Usage   : $obj->distance($newval)
 Function: Get/set the distance between the two terms connected
           by this path.

           Note that modifying the distance may not be meaningful. The
           implementation here is not connected to any graph engine,
           so changing an existing value may simply render the
           attribute's value wrong.

 Example : 
 Returns : value of distance (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub distance{
    my $self = shift;

    return $self->{'_distance'} = shift if @_;
    return $self->{'_distance'};
}

=head2 to_string

 Title   : to_string()
 Usage   : print $rel->to_string();
 Function: to_string method for Path.
 Returns : A string representation of this Path.
 Args    :

=cut

sub to_string {
    my( $self ) = @_;
    
    my $s = $self->SUPER::to_string();
    $s .= "-- Distance:\n";
    $s .= $self->distance() if defined($self->distance());
    $s .= "\n";
    
    return $s;
    
} # to_string



1;
