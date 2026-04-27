#!/usr/bin/env perl

my $version = "0.4.1";

## Changes from v0.4
## Will search for blast in PATH first
## Defaults to searching for schemes directory in the script directory
## Remove requirement for sequences.txt file pointing to location of allele sequences. Assumes they are in the same folder as the profiles.txt file

## Changes from v0.3.1
## Allow folder of all MLST schemes to be designated for alleles and profiles rather than manually entering each time

## Changes from v0.3
## Fixed bug where missing allele sequences were not being recorded
## Fixed bug where novel allele sequence files were not outputting correctly
## Remove makeblastdb dependency
## Add a pseudo-random identifier to the temp_alleles.fasta file so that the program can be run on multiple inputs in parallel
## Remove the temp_alleles.fasta file upon finishing

## Changes from v0.2.1
## Blast against all MLST targets at once, rather than one at a time. Should be MUCH faster.

## Changes from v0.2
## Modified to be more flexible and be able to do profiling from other PubMLST allele/profile combinations

## Changes from v0.1
## Automatically correct for issue in PA acs allele definitions where alleles 38 and 172 are both present in isolates with ST235.
## Make blast location modifiable from command line
## Make alignment slop user-definable

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions qw ( catfile path );

$|++;

my $blastdir = dirname(is_path("blastn"));
my $scriptpath = dirname(abs_path($0));
my $schemedir = "$scriptpath/schemes";
## Each organism / scheme should have its own directory in the schemedir.
## Directory names should have no spaces or special characters other than "_" or "-".
## Each organism / scheme directory should have the following files, named exactly as follows:
##   1. Allele sequence files. The file names should correspond to the allele names in the 'profiles.txt' file and have the suffix '.fas'.
##      For example:
##      acsA.fas
##      aroE.fas	
##      guaA.fas
##      mutL.fas
##      nuoD.fas
##      ppsA.fas
##      trpE.fas
##   2. 'profiles.txt' - tab-delimited text file of MLST definitions. Header line should start with "ST" followed by gene names matching the sequence file prefixes. 
##                       Subsequent columns are optional. Each line should contain an ST value and the corresponding gene allele sequence identifiers
##      For example:
##      ST	acsA	aroE	guaA	mutL	nuoD	ppsA	trpE	clonal_complex
##      1	1	1	1	1	1	1	1   CC1
##      2	6	130	1	3	4	4	17  CC2
##      3	9	22	110	13	15	6	5   CC1
##      4	17	137	1	3	11	12	111 CC3


my $usage = "
mlst_profiler.pl (version $version)

Determines the MLST profile from a fasta sequence.

Required:
  -o    Organism / scheme to use. To list all available organisms / schemes,
        run with the '-l' option.
        (See comments in the perl script for more info about how to add
        organisms / schemes)
    OR
  -p    File with MLST profiles. Format should be the same as profile files
        downloaded from pubmlst.org, i.e. a header line with the allele names
        starting with the sequence type ID, separated by tabs. Following lines
        list sequence types and allele IDs, separated by tabs.
        
        Example:
        ST	adk	atpA	dxr	glyA	recA	sodA	tpi
        1	1	1	1	10	1	3	5
        2	1	1	2	1	5	3	1
        3	1	1	2	1	1	1	1
        ...  
  -a    File with list of allele sequence files, one per line.
        Format should be: allele name (matching the column headers in the MLST
        profile above), tab, path to fasta file with allele sequences. Allele
        headers should be the allele name, underscore, allele number
        (for example: adk_1 is the first variant of the adk allele).
        File example:
        adk /path/to/adk.fas
        atpA    /path/to/atpA.fas
        ...

  -f    Query fasta sequence file. Can be gzipped.
    OR
  -F    list of query fasta sequence files, one path to sequence file per line.
        Can be gzipped.
        Optional: sequence ID can follow the path, separated by a tab.
          Example: /path/to/file.fasta<tab>genome1

Options:
  -s    slop, i.e. number of bases by which a potential novel allele can differ
        from its most closely related allele and not be considered potentially
        truncated
        (default: 2)
  -n    prefix for sequences of potential novel alleles. If this is given, will
        keep track of potential novel alleles not found in the given allele
        sequence files and output their sequences.
        (default: novel allele sequences will not be tracked)
  -c    MLST scheme gene targets
        (default: 7)

  -l    list available MLST profiles for use with the -o option
  -b    path to Blast+ binaries
        (default: '$blastdir')
  -d    path to organism/scheme directory. 
        (default: '$schemedir')

";

use Getopt::Std;
use vars qw( $opt_p $opt_a $opt_f $opt_F $opt_b $opt_s $opt_n $opt_o $opt_c $opt_d $opt_l );
getopts('p:a:f:F:b:s:n:o:d:l');

my $qfile       = $opt_f if $opt_f;
my $qlistfile   = $opt_F if $opt_F;
my $slopleng    = $opt_s ? $opt_s : 2;
$blastdir       = $opt_b ? $opt_b : $blastdir;
my $npref       = $opt_n if $opt_n;
my $targs       = $opt_c ? $opt_c : 7;
$schemedir      = $opt_d ? $opt_d : $schemedir;

die $usage unless ($opt_o or ($opt_p and $opt_a) or $opt_l);
if ($opt_o and $opt_p){
    print STDERR "WARNING: Both an organism / schema ($opt_o) and profile + sequence file paths (-p and -a) were given\n\tWill use the file paths provided.\n";
    $opt_o = "";
}

if ($blastdir){
    if (! -x "$blastdir/blastn"){
        die "ERROR: No excecutable blastn found at '$blastdir': please install in your PATH or provide path to directory for option '-b'\n";
    }
} else {
    die "ERROR: Blast not found in PATH and no blast directory provided with '-b'\n";
}

## Check the scheme directory for schemes
if ($opt_o or $opt_l){
    die "ERROR: option -d must be a directory\n" unless -d $schemedir; 
    opendir(my $dir, $schemedir);
    my %schemes;
    while (my $entry = readdir $dir){
        next unless -d "$schemedir/$entry";
        next if $entry =~ m/^\./;
        my $keep = 1;
        $keep = 0 unless -e "$schemedir/$entry/profiles.txt";
        $schemes{$entry} = 1 if $keep;
    }
    closedir $dir;
    if (scalar keys %schemes == 0){
        die "No MLST schemes found in directory '$schemedir'. Please see instructions for setting up scheme directory in the script or use -a and -p options to directly provide MLST profiles and sequences\n";
    }
    if ($opt_l){
        print "Available organisms / schemes:\n";
        foreach my $sc (sort keys %schemes){
            print "- $sc\n";
        }
        exit 0;
    }
    die "ERROR: $opt_o is not an available organism / scheme. Please use the '-l' option to list available organisms / schemes for -o option\n" unless $schemes{$opt_o};
    #$opt_a = "$schemedir/$opt_o/sequences.txt";
    $opt_p = "$schemedir/$opt_o/profiles.txt";
}
my $pfile = $opt_p;
my $afile = $opt_a if $opt_a;

die $usage unless ($opt_f or $opt_F);
die $usage if ($opt_f and $opt_F);

# read the profile matrix into a hash
my %profile;
my @allele_order;
my %seen_alleles;
my $header;
open (my $pin, "<$pfile") or die "ERROR: Can't open $pfile: $!\n";
my %non_alleles;
while (my $line = <$pin>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    my @tmp = split("\t", $line);
    if (!@allele_order){
        my $first = $tmp[0];
        die "ERROR: profile file '$pfile' does not appear to be formatted correctly (first line does not start with 'ST'). Please validate.\n" unless $first eq "ST";
        $header = "$first\t";
        for my $i (1 .. $#tmp){
            if ($i > $targs){
                $header .= "$tmp[$i]\t";
                next;
            }
            my $allele = $tmp[$i];
            $seen_alleles{$allele}++;
            push @allele_order, $allele;
        }
        chop $header;
        next;
    }
    my @vals;
    my $st = "$tmp[0]\t";
    for my $i (1 .. $#tmp){
        if ($i > $targs){
            $st .= "$tmp[$i]\t";
        }
        if ($non_alleles{$i}){
            $st .= "$tmp[$i]\t";
            next;
        }
        push @vals, $tmp[$i];
    }
    chop $st; #remove final tab character;
    my $val_string = join(",", @vals);
    die "ERROR: profile $val_string (ST $st) already exists as ST $profile{$val_string}\n" if $profile{$val_string};
    $profile{$val_string} = $st;
}
close ($pin);
my $non_count = scalar keys %non_alleles;

# Read in the allele file list and create a temporary sequence file with all alleles
my @file_list;
if ($afile){
    open (my $ain, "<$afile") or die "ERROR: Can't open $afile: $!\n";
    while (my $line = <$ain>){
        chomp $line;
        next if $line =~ m/^\s*$/;
        $line =~ s/\s*$//g;
        my ($id, $path) = split("\t", $line);
        $id =~ s/\s//g;
        push @file_list, ([$id,$path]);
    }
    close ($ain);
} else {
    foreach my $allele (@allele_order){
        push @file_list, ([$allele,"$schemedir/$opt_o/$allele.fas"]);
    }
}
my @chars = ("A".."Z", "a".."z", 0..9);
my $rstring;
$rstring .= $chars[rand @chars] for 1..8;
my $aout_file = "temp_alleles_$rstring.fasta";
open (my $aout, ">$aout_file");
my %afiles;
my $a_tot_count = 0;
foreach my $slice (@file_list){
    my ($id, $path) = @{$slice};
    unless (-e $path){
        die "ERROR: No sequence file found at '$path'\n";
    }
    open (my $tin, "<$path") or die "ERROR: Can't open $path: $!\n";
    while (my $tline = <$tin>){
        print $aout $tline;
        $a_tot_count++ if $tline =~ m/^>/;
    }
    close($tin);
    $afiles{$id} = 1;
}
close ($aout);
#my $status = system("$blastdir/makeblastdb -in temp_alleles.fasta -dbtype nucl >/dev/null");
#die "ERROR: makeblastdb exited with status $status\n" if $status;

## multiply the a_tot_count by 100, just to try to be sure to get all alignments
$a_tot_count = $a_tot_count * 100;

print "file\t", join("\t", @allele_order), "\t$header\n"; #header
my @qlist;
if ($opt_F){
    open (my $fin, "<$qlistfile") or die "ERROR: Can't open $qlistfile: $!\n";
    while (my $line = <$fin>){
        chomp $line;
        next if $line =~ m/^\s*$/;
        my ($path, $id) = split("\t", $line);
        $path =~ s/\s*$//;
        $path =~ s/^\s*//;
        unless ($id){
            $id = $path;
            $id =~ s/\/*[^\/]+\///g; #remove path
            $id =~ s/\.gz$//; #remove gz
            $id =~ s/\.[^.]+$//; #remove last suffix
        }
        $path =~ s/\|/\\\|/;
        push @qlist, ([$path, $id]);
    }
} else {
    my $id = $qfile;
    $id =~ s/\/*[^\/]+\///g; #remove path
    $id =~ s/\.gz$//; #remove gz
    $id =~ s/\.[^.]+$//; #remove last suffix
    $qfile =~ s/\|/\\\|/;
    push @qlist, ([$qfile, $id]);
}

my %novel_alleles;
my %nfiles;
foreach my $slice (@qlist){
    my ($file, $fileid) = @{$slice};
    #print STDERR "PROFILING $fileid...\n";
    
    # If the query file is gzipped, unzip it to a temporary file
    my $isgzip;
    if ($file =~ m/\.gz$/){
        system("gzip -cd $file > qtemp.fasta");
        $file = "qtemp.fasta";
        $isgzip = 1;
    }
    
    # Blast against the alleles to get the allele IDs
    #print STDERR "Blasting...\n";
    my %hits;
    open (my $bin, "$blastdir/blastn -query $file -subject $aout_file -max_target_seqs $a_tot_count -outfmt '6 qseqid qlen sseqid slen qstart qend sstart send pident' | ") or die "ERROR: blastn failed: $!\n";
    while (my $line = <$bin>){
        chomp $line;
        my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $pident) = split("\t", $line);
        (my $aid) = $sseqid =~ m/^(\S+)_\d+$/;
        push @{$hits{$aid}}, $line;
    }
    close ($bin);
    
    # Determine alleles.
    my @allele_vals;
    foreach my $allele (@allele_order){
        die "ERROR: No allele sequence file given for $allele\n" unless $afiles{$allele};
        #print STDERR "Determining $allele...\n";
        my $allele_num;
        my ($top_allele, $top_qid, $top_qstart, $top_qend, $top_pident, $top_line);
        my $trunc_hit;
        my $out_allele = "X";
        if ($hits{$allele}){
            foreach my $line (@{$hits{$allele}}){
                my ($qseqid, $qlen, $sseqid, $slen, $qstart, $qend, $sstart, $send, $pident) = split("\t", $line);
                my $aleng = abs($send - $sstart) + 1;
                (my $anum) = $sseqid =~ m/_(\d+)$/;
                if ($pident == 100){
                    if ($aleng == $slen){
                        if (!$allele_num){
                            ($allele_num, $top_allele, $top_qid, $top_qstart, $top_qend, $top_pident, $top_line) = ($anum, $anum, $qseqid, $qstart, $qend, $pident, $line);
                            next;
                        } else {
                            if ($anum eq $allele_num){
                                #print STDERR "Two identical alleles for $allele exist:\n\t$line\n";
                                next;
                            } else {
                                if ($allele_num =~ m/N/){
                                    ($allele_num, $top_allele, $top_qid, $top_qstart, $top_qend, $top_pident, $top_line) = ($anum, $anum, $qseqid, $qstart, $qend, $pident, $line);
                                } else {
                                    #print STDERR "At least two different alleles for $allele exist:\n\t$line\n";
                                    my @tmp = split("&", $allele_num);                          
                                    push @tmp, $anum;
                                    @tmp = grep {s/(^|\D)0+(\d)/$1$2/g,1} sort grep {s/(\d+)/sprintf"%09.9d",$1/ge,1} @tmp; #natural sort. See www.perlmonks.org/?node_id=68185
                                    $allele_num = join("&", @tmp);
                                }
                                next;
                            }
                        }
                    } else {
                        #print STDERR "Allele $allele is truncated\n\t$line\n" unless $allele_num;
                        $trunc_hit = 1;
                    }
                } else {
                    if ($top_qid){
                        next if ($qseqid eq $top_qid and $qstart <= $top_qend and $qend >= $top_qstart); #skip if this overlaps a previous hit
                        #print STDERR "Second sub-ideal hit found: $line\n"; ## Not worth keeping track of if there's already a top hit
                    } else {
                        next if $trunc_hit;
                        if (slop($slen,$aleng,$slopleng)){ #hit can differ by up to 2 bases from allele length;
                            ($top_qid, $top_qstart, $top_qend, $top_line) = ($qseqid, $qstart, $qend, $line);
                            #print STDERR "Potential novel allele found: $line\n";
                            
                            # get novel allele sequence
                            if (defined $npref){
                                open (my $nin, "<$file");
                                my $seq;
                                my $readseq;
                                while (my $line = <$nin>){
                                    chomp $line;
                                    next if $line =~ m/^\s.*$/;
                                    if ($line =~ m/^>/){
                                        last if $readseq;
                                        my $id = substr($line, 1);
                                        $id =~ s/\s.*$//;
                                        $readseq = 1 if $id eq $qseqid;
                                        next;
                                    }
                                    next unless $readseq;
                                    $line =~ s/\s//g;
                                    $seq .= $line;
                                }
                                close ($nin);
                                if ($seq){
                                    my $dist = $qend - $qstart + 1;
                                    my $aseq = substr($seq, $qstart - 1, $dist);
                                    if ($sstart > $send){
                                        $aseq = reverse $aseq;
                                        $aseq =~ tr/ACGTacgt/TGCAtgca/;
                                    }
                                    my $num = 0;
                                    $num = $novel_alleles{$allele}{"num"} if $novel_alleles{$allele}{"num"};
                                    if ($novel_alleles{$allele}{$aseq}){
                                        $num = $novel_alleles{$allele}{$aseq};
                                    } else {
                                        $num++;
                                        $novel_alleles{$allele}{$aseq} = $num;
                                        $novel_alleles{$allele}{"num"} = $num;
                                        my $nout;
                                        if ($nfiles{$allele}){
                                            $nout = $nfiles{$allele};
                                        } else {
                                            my $nfile = "$npref.$allele.fasta";
                                            open ($nout, ">$nfile");
                                            $nfiles{$allele} = $nout;
                                        }
                                        print $nout ">$allele\_N$num\n$aseq\n";
                                    }
                                    $allele_num = "N$num";
                                } else {
                                    die "ERROR: Problem trying to find sequence of novel allele: $line\n";
                                }
                            }
                            
                        } else {
                            #print STDERR "Potential truncated novel allele found: $line\n";
                            $out_allele = "NT";
                        }
                    }
                }
            }
            if ($allele_num){
                $out_allele = $allele_num;
            } else {
                if ($top_qid){
                    $out_allele = "N";
                } elsif (!$top_qid and $trunc_hit) {
                    $out_allele = "T";
                } else {
                    #print STDERR "No hits\n";
                }
            }
            push @allele_vals, $out_allele;
            #print STDERR "Top line: $top_line\n" if $top_line;
        } else {
            #print STDERR "No hits\n";
            push @allele_vals, $out_allele;
        }
    }
    
    if ($isgzip){
        unlink("$file");
    }
    
    # Determine sequence type
    my @tmpst = ("?") x $non_count;
    my $st = join("\t", @tmpst);
    my $pstring = join(",", @allele_vals);
    if ($profile{$pstring}){
        $st = $profile{$pstring};
    }
    
    # PA ST235 hack
    if ($st eq "?" and ($pstring eq "38&172,11,3,13,1,2,4" or $pstring eq "172&38,11,3,13,1,2,4")){
        $st = "235";
    }
    
    #if (%mlst_clade){
    #    $st .= "\t";
    #    if (my $mc = $mlst_clade{$pstring}){
    #        $st .= "$mc";
    #    }
    #}
    print "$fileid\t", join("\t", @allele_vals), "\t$st\n";
    print STDERR "\n";
}
unlink($aout_file);
foreach my $nout (keys %nfiles){
    close $nout;
}


#---------------------
sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($f) if -x $f;
        }
    }
    return();
}

sub slop{
    my $val1 = shift;
    my $val2 = shift;
    my $pm = shift;
    return(1) if abs($val1 - $val2) <= $pm;
    return(0);
}


