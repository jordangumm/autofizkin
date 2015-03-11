#!/usr/bin/env perl

use common::sense;
use autodie;
use Getopt::Long;
use Hurwitz::Utils qw(take);
use List::MoreUtils qw(uniq);
use Pod::Usage;
use Readonly;
use Statistics::Descriptive::Discrete;

main();

# --------------------------------------------------
sub main {
    my $loc_file    =  '';
    my $in_file     = '-';
    my $out_file    = '-';
    my $show_mode   =   1;
    my $unique_file =  '';
    my ($help, $man_page);
    GetOptions(
        'l|loc=s'      => \$loc_file,
        'i|in:s'       => \$in_file,
        'o|out:s'      => \$out_file,
        'm|show-mode!' => \$show_mode,
        'u|unique:s'   => \$unique_file,
        'help'         => \$help,
        'man'          => \$man_page,
    ) or pod2usage(2);

    if ($help || $man_page) {
        pod2usage({
            -exitval => 0,
            -verbose => $man_page ? 2 : 1
        });
    }; 

    unless ($loc_file) {
        pod2usage('No location file');
    } 

    open my $loc_fh, '<', $loc_file;

    my $in;
    if ($in_file eq '-') {
        $in = \*STDIN;
    }
    else {
        open $in, '<', $in_file;
    }

    my $out_fh;
    if ($out_file eq '-') {
        $out_fh = \*STDOUT;
    }
    else {
        open $out_fh, '>', $out_file;
    }

    my %seen;
    if ($unique_file) {
        open my $tmp, '<', $unique_file;
        while (<$tmp>) {
            chomp;
            $seen{$_} = 1;
        }
        close $tmp;
    }

    READ:
    while (my $loc = <$loc_fh>) {
        chomp($loc);
        my ($read_id, $n_kmers) = split /\t/, $loc;

        next READ if $unique_file && defined $seen{ $read_id };

        if (my $mode = mode(take($n_kmers, $in))) {
            $seen{ $read_id }++ if $unique_file;
            say $out_fh ($show_mode) ? join("\t", $read_id, $mode) : $read_id;
        }
    }

    if ($unique_file) {
        open my $tmp, '>', $unique_file;
        say $tmp join "\n", keys %seen;
        close $tmp;
    }
}

# --------------------------------------------------
sub mode {
    my @vals = @_ or return;
    my $mode = 0;

    if (scalar @vals == 1) {
        $mode = shift @vals;
    }
    else {
        my @distinct = uniq(@vals);

        if (scalar @distinct == 1) {
            $mode = shift @distinct;
        }
        else {
            my $stats = Statistics::Descriptive::Discrete->new;
            $stats->add_data(@vals);
            return $stats->mode();
        }
    }

    return $mode;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

jellyfish-reduce.pl - Jellyfish output to read modes

=head1 SYNOPSIS

  kmerizer.pl -i input.fasta -k 20 -l input.loc | \
  jellyfish query -i subject.jf | \
  jellyfish-reduce.pl -l input.locs -o subject.modes

Required Arguments:

  -l|--loc        Path to the location file (shows readId/numKmers)

Options:

  -i|-in          Path to kmers/counts or '-' for STDIN (default)
  -o|--out        Path to output file or '-' for STDOUT (default)
  -u|--unique     Name of the file to read/write unique readIds 
  -m|--show-mode  Show the mode value (default true)
                  Use '--no-m' or '--no-show-mode' to negate

  --help          Show brief help and exit
  --man           Show full documentation

=head1 DESCRIPTION

Calculate mode of reads from Jellyfish output.

=head1 AUTHOR

Ken Youens-Clark E<lt>kyclark@email.arizona.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2015 kyclark

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut