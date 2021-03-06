# $Id: QueryI.pm,v 1.1 2002/11/20 08:39:03 lstein Exp $
#
# BioPerl module for Bio::DB::QueryI.pm
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::DB::QueryI - Object Interface to queryable sequence databases

=head1 SYNOPSIS

   # using Bio::DB::Query::GenBank as an example
   my $query_string = 'Oryza[Organism] AND EST[Keyword]';
   my $query = Bio::DB::Query::GenBank->new(-db=>'nucleotide',
                                            -query=>$query_string);
   my $count = $query->count;
   my @ids   = $query->ids;

   # get a genbank database handle
   $gb = new Bio::DB::GenBank;
   my $stream = $db->get_Stream_by_query($query);
   while (my $seq = $stream->next_seq) {
      ...
   }

   # initialize the list yourself
   my $query = Bio::DB::Query::GenBank->new(-ids=>['X1012','CA12345']);

=head1 DESCRIPTION

This interface provides facilities for managing sequence queries such
as those offered by Entrez.  A query object is created by calling
new() with a database-specific argument list. From the query object
you can either obtain the list of IDs returned by the query, or a
count of entries that would be returned.  You can pass the query
object to a Bio::DB::RandomAccessI object to return the entries
themselves as a list or a stream.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via email or the
web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Lincoln Stein

Email lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::QueryI;
use strict;
use Bio::Root::RootI;

use vars qw(@ISA $VERSION);

@ISA = qw(Bio::Root::RootI);
$VERSION = '0.1';

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::QueryI->new(@args);
 Function: constructor
 Returns : QueryI object
 Args    : -query       a query string
           -ids         a list of ids as an arrayref

Create new QueryI object.  You may initialize with either a query
string or with a list of ids.  If both ids and a query are provided,
the former takes precedence.

Subclasses may recognize additional arguments.

=cut

=head2 count

 Title   : count
 Usage   : $count = $db->count;
 Function: return count of number of entries retrieved by query
 Returns : integer
 Args    : none

Returns the number of entries that are matched by the query.

=cut

sub count   {
  my $self = shift;
  my @ids  = $self->ids;
  scalar @ids;
}

=head2 ids

 Title   : ids
 Usage   : @ids = $db->ids([@ids])
 Function: get/set matching ids
 Returns : array of sequence ids
 Args    : (optional) array ref with new set of ids

=cut

sub ids     {
  my $self = shift;
  $self->throw_not_implemented;
}

=head2 query

 Title   : query
 Usage   : $query = $db->query([$query])
 Function: get/set query string
 Returns : string
 Args    : (optional) new query string

=cut

sub query   {
  my $self = shift;
  $self->throw_not_implemented;
}

1;
