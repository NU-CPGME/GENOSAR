#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "
taxonomy_to_sourmash.pl <taxonomy_table> </path/to/fasta/directory>

Converts ncbi datasets taxonomy table to sourmash table for 
use with lca index.

";

die $usage unless scalar(@ARGV) >= 2;

my ($taxfile, $fastadir) = (@ARGV);

## Read the taxonomy file
open (my $in, "<$taxfile") or die "ERROR: Can't open $taxfile: $!\n";
my %taxresults;
my %dups;
while (my $line = <$in>){
    chomp $line;
    my @tmp = split("\t", $line);
    next if $tmp[0] eq "Query";
    #accession,taxid,superkingdom,phylum,class,order,family,genus,species,strain
    my $strain = $tmp[2];
    my $species = $tmp[24];
    #next unless $species;
    #next unless $species =~ m/^[A-Z]/; ## remove those starting with lower-case letters
    my $strain_old = $strain;
    $strain =~ s/^(\S+ \S+).*/$1/;
    if ($strain =~ m/sp\.$/){
        $strain = $strain_old;
        $strain =~ s/^(\S+ \S+ \S+).*/$1/;
    }
    my $species_old = $species;
    $species =~ s/^(\S+ \S+).*/$1/;
    #my $string = join(",",@tmp[25,12,14,16,18,20,22]) . ",$strain";
    my $string = join(",",@tmp[12,14,16,18,20,22]) . ",$strain";
    #$taxresults{$strain} = $string;
    #next if ($species =~ m/sp\.$/);
    $taxresults{"$strain"} = $string;
    push @{$dups{$string}}, $strain;
    #print STDERR "$species($species_old) --> $string\n";
}
close ($in);

## Check for duplicates (debug)
foreach my $string (sort keys %dups){
    my @tmp = @{$dups{$string}};
    if (scalar @tmp > 1){
        my %hash;
        foreach my $slice (@tmp){
            $hash{$slice} = 1;
        }
        if (scalar(keys %hash) > 1){
            print STDERR "$string:\n";
            foreach my $key (sort keys %hash){
                print STDERR "\t$key\n";
            }

        }
    }

}


## Read the files next
my @files = glob("$fastadir/*");
my @output;
foreach my $file (@files){
    my $order = basename($file);
    $order =~ s/\..+$//;
    open (my $in, "<$file") or next;
    #print STDERR "file:$file order:$order\n";
    while (my $line = <$in>){
        chomp $line;
        last unless $line =~ m/^>/;
        $line =~ m/^>(\S+) ([^,]+)/;
        my ($acc, $strain) = ($1, $2);
        $acc =~ s/\..*//;
        my $strain_old = $strain;
        $strain =~ s/^(\S+ \S+).*/$1/;
        if ($strain =~ m/sp\.$/){
            $strain = $strain_old;
            $strain =~ s/^(\S+ \S+ \S+).*/$1/;
        }
        #print STDERR "species:$species strain:$strain\n";
        my $seen;
        if ($taxresults{"$strain"}){
            my $string = $taxresults{"$strain"};
            my $newstring = "$acc,$string";
            push @output, ([$order, $newstring]);
            $seen = 1;
        }
        unless ($seen){
            print STDERR "$order\t$acc\t$strain\n";
        }

        # my $seen;
        # my @tmp = split(" ", $strain);
        # for my $i (reverse(0 .. $#tmp)){
        #     my $test = join(" ", @tmp[0 .. $i]);
        #     #print STDERR "\t$test\n";
        #     if ($taxresults{$test}){
        #         my $string = $taxresults{$test};
        #         my $newstring = "$acc,$string,$test";
        #         push @output, ([$order, $newstring]);
        #         $seen = 1;
        #         last;
        #     }
        # }
        # if (!$seen){
        #     $strain =~ s/str\. //;
        #     $strain =~ s/strain //;
        #     my @tmp = split(" ", $strain);
        #     for my $i (reverse(0 .. $#tmp)){
        #         my $test = join(" ", @tmp[0 .. $i]);
        #         #print STDERR "\t$test\n";
        #         if ($taxresults{$test}){
        #             my $string = $taxresults{$test};
        #             my $newstring = "$acc,$string,$test";
        #             push @output, ([$order, $newstring]);
        #             $seen = 1;
        #             last;
        #         }
        #     }
        # }
        # unless ($seen){
        #     print STDERR "$order\t$acc\t$strain\n";
        # }
    }
    close ($in);
}

@output = sort{$a->[0] <=> $b->[0]} @output;
#print "accession,taxid,superkingdom,phylum,class,order,family,genus,species,strain\n";
print "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n";
for my $i (0 .. $#output){
    print "$output[$i][1],\n";
}

# 1: Query
# 2: Taxid
# 3: Tax name
# 4: Authority
# 5: Rank
# 6: Basionym
# 7: Basionym authority
# 8: Curator common name
# 9: Has type material
# 10: Group name
# 11: Domain/Realm name
# 12: Domain/Realm taxid
# 13: Kingdom name
# 14: Kingdom taxid
# 15: Phylum name
# 16: Phylum taxid
# 17: Class name
# 18: Class taxid
# 19: Order name
# 20: Order taxid
# 21: Family name
# 22: Family taxid
# 23: Genus name
# 24: Genus taxid
# 25: Species name
# 26: Species taxid
# 27: Scientific name is formal


