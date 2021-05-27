#!/usr/local/bin/perl -w

# lightweight fasta reader capabilities:
package Fasta_reader;

use strict;

sub new {
    my ($packagename, $fastaFile) = @_;
    my $self = { fastaFile => $fastaFile,
		 fileHandle => undef };
    bless ($self, $packagename);
    
    ## create filehandle
    my $filehandle = undef;
    open ($filehandle, $fastaFile) or die "Error: Couldn't open $fastaFile\n";
    $self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    my $orig_record_sep = $/;
    $/="\n>";
    my $filehandle = $self->{fileHandle};
    my $next_text_input = <$filehandle>;
    my $seqobj = undef;
    if ($next_text_input) {
	$next_text_input =~ s/^>|>$//g; #remove trailing > char.
	$next_text_input =~ tr/\t\n\000-\037\177-\377/\t\n/d; #remove cntrl chars
	my ($header, @seqlines) = split (/\n/, $next_text_input);
	my $sequence = join ("", @seqlines);
	$sequence =~ s/\s//g;
		
	$seqobj = Sequence->new($header, $sequence);
    }
    
    $/ = $orig_record_sep; #reset the record separator to original setting.
    
    return ($seqobj); #returns null if not instantiated.
}


#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}




##############################################
package Sequence;
use strict;

sub new {
    my ($packagename, $header, $sequence) = @_;
    
    ## extract an accession from the header:
    my ($acc, $rest) = split (/\s+/, $header, 2);
        
    my $self = { accession => $acc,
		 header => $header,
		 sequence => $sequence,
		 filename => undef };
    bless ($self, $packagename);
    return ($self);
}

####
sub get_accession {
    my $self = shift;
    return ($self->{accession});
}

####
sub get_header {
    my $self = shift;
    return ($self->{header});
}

####
sub get_sequence {
    my $self = shift;
    return ($self->{sequence});
}

#### 
sub get_FASTA_format {
    my $self = shift;
    my $header = $self->get_header();
    my $sequence = $self->get_sequence();
    $sequence =~ s/(\S{60})/$1\n/g;
    my $fasta_entry = ">$header\n$sequence\n";
    return ($fasta_entry);
}


####
sub write_fasta_file {
    my $self = shift;
    my $filename = shift;

    my ($accession, $header, $sequence) = ($self->{accession}, $self->{header}, $self->{sequence});
    $sequence =~ s/(\S{60})/$1\n/g;
   
    my $tempfile;
    if ($filename) {
	$tempfile = $filename;
    } else {
	my $acc = $accession;
	$acc =~ s/\W/_/g;
	$tempfile = "$acc.fasta";
    }
    
    open (TMP, ">$tempfile") or die "ERROR! Couldn't write a temporary file in current directory.\n";
    print TMP ">$header\n$sequence";
    close TMP;
    return ($tempfile);
}

1; #EOM


