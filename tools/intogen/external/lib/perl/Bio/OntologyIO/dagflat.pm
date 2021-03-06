# $Id: dagflat.pm,v 1.2.2.6 2003/06/30 05:04:06 lapp Exp $
#
# BioPerl module for Bio::OntologyIO::dagflat
#
# Cared for by Hilmar Lapp, hlapp at gmx.net
#
# (c) Christian M. Zmasek, czmasek@gnf.org, 2002.
# (c) Hilmar Lapp, hlapp at gmx.net, 2003.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
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

dagflat - a base class parser for GO flat-file type formats 

=head1 SYNOPSIS

  use Bio::OntologyIO;

  # do not use directly -- use via Bio::OntologyIO
  # e.g., the GO parser is a simple extension of this class
  my $parser = Bio::OntologyIO->new
	( -format       => "go",
          -defs_file    => "/home/czmasek/GO/GO.defs",
	  -files        => ["/home/czmasek/GO/component.ontology",
	                    "/home/czmasek/GO/function.ontology",
	                    "/home/czmasek/GO/process.ontology"] );

  my $go_ontology = $parser->next_ontology();

  my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
  my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );

=head1 DESCRIPTION

Needs Graph.pm from CPAN.

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

Christian M. Zmasek

Email: czmasek@gnf.org  or  cmzmasek@yahoo.com

WWW:   http://www.genetics.wustl.edu/eddy/people/zmasek/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head2 CONTRIBUTOR

 Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package  Bio::OntologyIO::dagflat;

use vars qw( @ISA );
use strict;

use Bio::Root::IO;
use Bio::Ontology::SimpleGOEngine;
use Bio::Ontology::Ontology;
use Bio::Ontology::TermFactory;
use Bio::OntologyIO;

use constant TRUE         => 1;
use constant FALSE        => 0;


@ISA = qw( Bio::OntologyIO );


=head2 new

 Title   : new
 Usage   : $parser = Bio::OntologyIO->new(
                             -format => "go",
                             -defs_file => "/path/to/GO.defs",
                             -files => ["/path/to/component.ontology",
                                        "/path/to/function.ontology",
                                        "/path/to/process.ontology"] );
 Function: Creates a new dagflat parser.
 Returns : A new dagflat parser object, implementing L<Bio::OntologyIO>.
 Args    : -defs_file  => the name of the file holding the term
                          definitions
           -files      => a single ontology flat file holding the
                          term relationships, or an array ref holding
                          the file names (for GO, there will usually be
                          3 files: component.ontology, function.ontology,
                          process.ontology)
           -file       => if there is only a single flat file, it may
                          also be specified via the -file parameter
           -ontology_name => the name of the ontology; if not specified the
                          parser will auto-discover it by using the term
                          that starts with a '$', and converting underscores
                          to spaces
           -engine     => the Bio::Ontology::OntologyEngineI object
                          to be reused (will be created otherwise); note
                          that every Bio::Ontology::OntologyI will
                          qualify as well since that one inherits from the
                          former.

=cut

# in reality, we let OntologyIO::new do the instantiation, and override
# _initialize for all initialization work
sub _initialize {
    my ($self, @args) = @_;
    
    $self->SUPER::_initialize( @args );

    my ( $defs_file_name,$files,$name,$eng ) =
	$self->_rearrange([qw( DEFS_FILE
			       FILES
			       ONTOLOGY_NAME
			       ENGINE)
			   ],
			  @args );
    
    $self->_done( FALSE );
    $self->_not_first_record( FALSE );
    $self->_term( "" );
    delete $self->{'_ontologies'};

    # ontology engine (and possibly name if it's an OntologyI)
    $eng = Bio::Ontology::SimpleGOEngine->new() unless $eng;
    if($eng->isa("Bio::Ontology::OntologyI")) {
	$self->ontology_name($eng->name());
	$eng = $eng->engine() if $eng->can('engine');
    }
    $self->_ont_engine($eng);
    
    # flat files to parse
    $self->defs_file( $defs_file_name );
    $self->{_flat_files} = $files ? ref($files) ? $files : [$files] : [];

    # ontology name (overrides implicit one through OntologyI engine)
    $self->ontology_name($name) if $name;

} # _initialize

=head2 ontology_name

 Title   : ontology_name
 Usage   : $obj->ontology_name($newval)
 Function: Get/set the name of the ontology parsed by this module. 
 Example : 
 Returns : value of ontology_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub ontology_name{
    my $self = shift;

    return $self->{'ontology_name'} = shift if @_;
    return $self->{'ontology_name'};
}


=head2 parse

 Title   : parse()
 Usage   : $parser->parse();   
 Function: Parses the files set with "new" or with methods
           defs_file and _flat_files.

           Normally you should not need to call this method as it will
           be called automatically upon the first call to
           next_ontology().

 Returns : [Bio::Ontology::OntologyEngineI]
 Args    :

=cut

sub parse {
    my $self = shift;
    
    # setup the default term factory if not done by anyone yet
    $self->term_factory(Bio::Ontology::TermFactory->new(
					     -type => "Bio::Ontology::Term"))
	unless $self->term_factory();

    # create the ontology object itself
    my $ont = Bio::Ontology::Ontology->new(-name => $self->ontology_name(),
					   -engine => $self->_ont_engine());

    # parse definitions
    while( my $term = $self->_next_term() ) {
        $self->_add_term( $term, $ont );
    }

    # set up the ontology of the relationship types
    foreach ($self->_part_of_relationship(), $self->_is_a_relationship()) {
	$_->ontology($ont);
    }

    # pre-seed the IO system with the first flat file if -file wasn't provided
    if(! $self->_fh) {
	$self->_initialize_io(-file => shift(@{$self->_flat_files()}));
    }

    while($self->_fh) {
	$self->_parse_flat_file($ont);
	# advance to next flat file if more are available
	if(@{$self->_flat_files()}) {
	    $self->close();
	    $self->_initialize_io(-file => shift(@{$self->_flat_files()}));
	} else {
	    last; # nothing else to parse so terminate the loop
	}
    }
    $self->_add_ontology($ont);
    
    # not needed anywhere, only because of backward compatibility
    return $self->_ont_engine();
    
} # parse

=head2 next_ontology

 Title   : next_ontology
 Usage   :
 Function: Get the next available ontology from the parser. This is the
           method prescribed by Bio::OntologyIO.
 Example :
 Returns : An object implementing Bio::Ontology::OntologyI, and undef if
           there is no more ontology in the input.
 Args    :


=cut

sub next_ontology{
    my $self = shift;

    # parse if not done already
    $self->parse() unless exists($self->{'_ontologies'});
    # return next available ontology
    return shift(@{$self->{'_ontologies'}}) if exists($self->{'_ontologies'});
    return undef;
}


=head2 defs_file

 Title   : defs_file
 Usage   : $parser->defs_file( "GO.defs" );
 Function: Set/get for the term definitions filename.
 Returns : The term definitions file name [string].
 Args    : On set, the term definitions file name [string] (optional).

=cut

sub defs_file {
    my $self = shift;

    if ( @_ ) {
	my $f = shift;
        $self->{ "_defs_file_name" } = $f;
	$self->_defs_io->close() if $self->_defs_io();
	if(defined($f)) {
            $self->_defs_io( Bio::Root::IO->new( -input => $f ) );
        }
    }
    return $self->{ "_defs_file_name" };
} # defs_file


=head2 _flat_files

 Title   : _flat_files
 Usage   : $files_to_parse = $parser->_flat_files();
 Function: Get the array of ontology flat files that need to be parsed.

           Note that this array will decrease in elements over the
           parsing process. Therefore, it's value outside of this
           module will be limited. Also, be careful not to alter the
           array unless you know what you are doing.

 Returns : a reference to an array of zero or more strings
 Args    : none

=cut

sub _flat_files {
    my $self = shift;

    $self->{_flat_files} = [] unless exists($self->{_flat_files});
    return $self->{_flat_files};
}


# INTERNAL METHODS
# ----------------

=head2 _defs_io

 Title   : _defs_io
 Usage   : $obj->_defs_io($newval)
 Function: Get/set the Bio::Root::IO instance representing the
           definition file, if provided (see defs_file()).
 Example : 
 Returns : value of _defs_io (a Bio::Root::IO object)
 Args    : on set, new value (a Bio::Root::IO object or undef, optional)


=cut

sub _defs_io{
    my $self = shift;

    return $self->{'_defs_io'} = shift if @_;
    return $self->{'_defs_io'};
}

sub _add_ontology {
    my $self = shift;
    $self->{'_ontologies'} = [] unless exists($self->{'_ontologies'});
    foreach my $ont (@_) {
	$self->throw(ref($ont)." does not implement Bio::Ontology::OntologyI")
	    unless ref($ont) && $ont->isa("Bio::Ontology::OntologyI");
	# the ontology name may have been auto-discovered while parsing
	# the file
	$ont->name($self->ontology_name) unless $ont->name();
	push(@{$self->{'_ontologies'}}, $ont);
    }
}

# This simply delegates. See SimpleGOEngine.
sub _add_term {
    my ( $self, $term, $ont ) = @_;

    $term->ontology($ont) if $ont && (! $term->ontology);
    $self->_ont_engine()->add_term( $term );
} # _add_term 



# This simply delegates. See SimpleGOEngine
sub _part_of_relationship {
    my $self = shift;

    return $self->_ont_engine()->part_of_relationship(@_);
} # _part_of_relationship 



# This simply delegates. See SimpleGOEngine
sub _is_a_relationship {
    my $self = shift;

    return $self->_ont_engine()->is_a_relationship(@_);
} # _is_a_relationship 



# This simply delegates. See SimpleGOEngine
sub _add_relationship {
    my ( $self, $parent, $child, $type, $ont ) = @_;

    # note the triple terminology (subject,predicate,object) corresponds to
    # (child,type,parent)
    $self->_ont_engine()->add_relationship( $child, $type, $parent, $ont );


} # _add_relationship


# This simply delegates. See SimpleGOEngine
sub _has_term {
    my $self = shift;

    return $self->_ont_engine()->has_term( @_ );
} # _add_term



# This parses the relationships files
sub _parse_flat_file {
    my $self = shift;
    my $ont  = shift;
    
    my @stack       = ();
    my $prev_spaces = -1;
    my $prev_term   = "";
    
    while( my $line = $self->_readline() ) {
        
        if ( $line =~ /^!/ ) {
            next;
        }
        
        my $current_term   = $self->_get_first_termid( $line );
        my @isa_parents    = $self->_get_isa_termids( $line );
        my @partof_parents = $self->_get_partof_termids( $line );
        my @syns           = $self->_get_synonyms( $line );
        my @sec_go_ids     = $self->_get_secondary_termids( $line );
        my @cross_refs     = $self->_get_db_cross_refs( $line );
        
        
        if ( ! $self->_has_term( $current_term ) ) {
            my $term =$self->_create_ont_entry($self->_get_name($line,
								$current_term),
					       $current_term );
            $self->_add_term( $term, $ont );
        }
        
        my $current_term_object = $self->_ont_engine()->get_terms( $current_term );
        
        $current_term_object->add_dblink( @cross_refs );
        $current_term_object->add_secondary_id( @sec_go_ids );
        $current_term_object->add_synonym( @syns );
        unless ( $line =~ /^\$/ ) {
            $current_term_object->ontology( $ont );
        }
        foreach my $parent ( @isa_parents ) {
            if ( ! $self->_has_term( $parent ) ) {
                my $term = $self->_create_ont_entry($self->_get_name($line,
								     $parent),
						    $parent );
                $self->_add_term( $term, $ont );
            }
            
            $self->_add_relationship( $parent,
                                      $current_term,
                                      $self->_is_a_relationship(),
				      $ont);
             
        }
        foreach my $parent ( @partof_parents ) {
            if ( ! $self->_has_term( $parent ) ) {
                my $term = $self->_create_ont_entry($self->_get_name($line,
								     $parent),
						    $parent );
                $self->_add_term( $term, $ont );
            }
           
            $self->_add_relationship( $parent,
                                      $current_term,
                                      $self->_part_of_relationship(),
				      $ont);
        }
        
        my $current_spaces = $self->_count_spaces( $line );
        
        if ( $current_spaces != $prev_spaces  ) {
       
            if ( $current_spaces == $prev_spaces + 1 ) {
                push( @stack, $prev_term ); 
            }
            elsif ( $current_spaces < $prev_spaces ) {
                my $n = $prev_spaces -  $current_spaces;
                for ( my $i = 0; $i < $n; ++$i ) {
                    pop( @stack );
                }
            }
            else {
		$self->throw( "format error (file ".$self->file.")" );
            } 
        }
        
        my $parent = $stack[ @stack - 1 ];
        
        # add a relationship if the line isn't the one with the root term
	# of the ontology (which is also the name of the ontology)
        if ( index($line,'$') != 0 ) {
	    if ( $line !~ /^\s*[<%]/ ) {
		$self->throw( "format error (file ".$self->file.")" );
	    }
	    my $reltype = ($line =~ /^\s*</) ?
		$self->_part_of_relationship() :
		$self->_is_a_relationship();
            $self->_add_relationship( $parent, $current_term, $reltype,
				      $ont);
        }
        
        $prev_spaces = $current_spaces;
        
        $prev_term = $current_term;
        
    } 
    return $ont;
} # _parse_relationships_file



# Parses the 1st term id number out of line.
sub _get_first_termid {
    my ( $self, $line ) = @_;
    
    if ( $line =~ /;\s*([A-Z]{1,8}:\d{3,})/ ) {
        return $1;
    }
    else {
        $self->throw( "format error: no term id in line \"$line\"" );
    }
    
} # _get_first_termid



# Parses the name out of line.
sub _get_name {
    my ( $self, $line, $termid ) = @_;
    
    if ( $line =~ /([^;<%]+);\s*$termid/ ) {
        my $name = $1;
	# remove trailing and leading whitespace
        $name =~ s/\s+$//;
        $name =~ s/^\s+//;
	# remove leading dollar character; also we default the name of the
	# ontology to this name if preset to something else
	if(index($name,'$') == 0) {
	    $name = substr($name,1);
	    # replace underscores by spaces for setting the ontology name
	    $self->ontology_name(join(" ",split(/_/,$name))) 
		unless $self->ontology_name();
	}
        return $name;
    }
    else {
        return undef;
    }
} # _get_name   


# Parses the synonyms out of line.
sub _get_synonyms {
    my ( $self, $line ) = @_;
   
    my @synonyms = ();
   
    while ( $line =~ /synonym\s*:\s*([^;^<^%]+)/g ) {
        my $syn = $1;
        $syn =~ s/\s+$//;
        $syn =~ s/^\s+//;
        push( @synonyms, $syn );
    }
    return @synonyms;
    
} # _get_synonyms



# Parses the db cross refs out of line.
sub _get_db_cross_refs {
    my ( $self, $line ) = @_;
   
    my @refs = ();
   
    while ( $line =~ /;([^;^<^%^:]+:[^;^<^%^:]+)/g ) {
        my $ref = $1;
        if ( $ref =~ /synonym/ || $ref =~ /[A-Z]{1,8}:\d{3,}/ ) {
            next;
        }
        $ref =~ s/\s+$//;
        $ref =~ s/^\s+//;
        push( @refs, $ref );
    }
    return @refs;
    
} 


# Parses the secondary go ids out of a line
sub _get_secondary_termids {
    my ( $self, $line ) = @_;
    my @secs = ();
   
    while ( $line =~ /,\s*([A-Z]{1,8}:\d{3,})/g ) {
        my $sec = $1;
        push( @secs, $sec );
    }
    return @secs;
    
} # _get_secondary_termids 



# Parses the is a ids out of a line
sub _get_isa_termids {
    my ( $self, $line ) = @_;
    
    my @ids = ();
    
    $line =~ s/[A-Z]{1,8}:\d{3,}//;
    
    while ( $line =~ /%[^<^,]*?([A-Z]{1,8}:\d{3,})/g ) {
        push( @ids, $1 );
    }
    return @ids; 
} # _get_isa_termids



# Parses the part of ids out of a line
sub _get_partof_termids {
    my ( $self, $line ) = @_;
    
    my @ids = ();
    
    $line =~ s/[A-Z]{1,8}:\d{3,}//;
    
    while ( $line =~ /<[^%^,]*?([A-Z]{1,8}:\d{3,})/g ) {
        push( @ids, $1 );
    }
    return @ids; 
    
    
} # _get_partof_termids




# Counts the spaces at the beginning of a line in the relationships files
sub _count_spaces {
    my ( $self, $line ) = @_;
     
    if ( $line =~ /^(\s+)/ ) {
         return length( $1 );
    }
    else {
         return 0;
    }
} # _count_spaces




# "next" method for parsing the defintions file
sub _next_term {
    my ( $self ) = @_;

    if ( ($self->_done() == TRUE) || (! $self->_defs_io())) {
        return undef;
    }
    
    my $line      = "";
    my $termid    = "";
    my $next_term = $self->_term();
    my $def       = "";
    my $comment   = "";
    my @def_refs  = ();
    my $isobsolete;
    
    while( $line = ( $self->_defs_io->_readline() ) ) {
    
        if ( $line !~ /\S/ 
        ||   $line =~ /^\s*!/ ) {
            next;
        }
        elsif ( $line =~ /^\s*term:\s*(.+)/ ) {
	    $self->_term( $1 );
            last if $self->_not_first_record();
            $next_term = $1;
	    $self->_not_first_record( TRUE );
        }
        elsif ( $line =~ /^\s*[a-z]{1,8}id:\s*(.+)/ ) {
            $termid = $1;
        }
        elsif ( $line =~ /^\s*definition:\s*(.+)/ ) {
            $def = $1;   
	    $isobsolete = 1 if index($def,"OBSOLETE") == 0;
        }
        elsif ( $line =~ /^\s*definition_reference:\s*(.+)/ ) {
            push( @def_refs, $1 );  
        }
        elsif ( $line =~ /^\s*comment:\s*(.+)/ ) {
            $comment = $1;  
        }
    }
    $self->_done( TRUE ) unless $line; # we'll come back until done
    
    return $self->_create_ont_entry( $next_term, $termid, $def,
				     $comment, \@def_refs, $isobsolete);
} # _next_term





# Holds the GO engine to be parsed into
sub _ont_engine {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_ont_engine" } = $value;
    }
    
    return $self->{ "_ont_engine" };
} # _ont_engine




# Used to create ontology terms.
# Arguments: name, id
sub _create_ont_entry {
    my ( $self, $name, $termid, $def, $cmt, $dbxrefs, $obsolete ) = @_;

    if((!defined($obsolete)) && (index(lc($name),"obsolete") == 0)) {
	$obsolete = 1;
    }
    my $term = $self->term_factory->create_object(-name => $name,
						  -identifier => $termid,
						  -definition => $def,
						  -comment => $cmt,
						  -dblinks => $dbxrefs,
						  -is_obsolete => $obsolete);

    return $term;
} # _create_ont_entry



# Holds whether first record or not
sub _not_first_record {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_not_first_record" } = $value;
    }
    
    return $self->{ "_not_first_record" };
} # _not_first_record



# Holds whether done or not
sub _done {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_done" } = $value;
    }
    
    return $self->{ "_done" };
} # _done


# Holds a term.  
sub _term {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_term" } = $value;
    }
    
    return $self->{ "_term" };
} # _term  
  
  
