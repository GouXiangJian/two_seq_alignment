#!/usr/bin/env perl

#author : Xiangjian Gou
#date   : 2018/12/28

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw {gettimeofday};

#get the options
my ($fasta, $output, $match, $unmatch, $space, $help);
GetOptions(
              'fasta=s'   =>   \$fasta,
             'output=s'   =>   \$output,
              'match=i'   =>   \$match,
            'unmatch=i'   =>   \$unmatch,
              'space=i'   =>   \$space,
                'help+'   =>   \$help,
);

#describe program information
my $usage = <<__GUIDE__;

################################################################################
Function:
  Local alignment between two sequences by using Smith-Waterman algorithm.

Usage:
  perl Smith_Waterman.pl -f <fasta> -o <output> -m <match> -u <unmatch> -s <space>

Required Options: 
  -f|-fasta   <file> : fasta file contains 2 sequences

Optional Options:
  -o|-output  <file> : a file name to store the output (default: STDOUT)
  -m|-match   <int>  : a score when character match    (default:  1)
  -u|-unmatch <int>  : a score when character unmatch  (default: -1)
  -s|-space   <int>  : a score when character space    (default: -2)
  -h|-help           : show the help information
################################################################################

__GUIDE__

#show the help information
die $usage if $help;

#lack of required options
die $usage unless $fasta;

#set the default value
$match   =  1  unless defined $match;
$unmatch = -1  unless defined $unmatch;
$space   = -2  unless defined $space;

#set the output place
if (defined $output) {
    open my $out, '>', $output or die "can't generate $output:$!";
    select $out;
}else{
    select STDOUT;
}

#======================================================================
#start main program
#======================================================================

#step1 : start local alignment
my $id_seq = GetFastaInfo($fasta);
my @ids    = sort keys %$id_seq;
my $time1  = join '.', gettimeofday; #record the start time
my ($max_align_score, $align_seq_1, $align_seq_2) = SmithWaterman( $id_seq->{$ids[0]}, $id_seq->{$ids[1]} );
my $time2  = join '.', gettimeofday; #record the   end time

#step2 : add align line between the two sequences
my $align_line = '';
foreach my $sub ( 0 .. length($align_seq_1)-1 ) {
    $align_line .= substr($align_seq_1, $sub, 1) eq substr($align_seq_2, $sub, 1) ? '|' : ' ';
}

#step3 : get the number and ratio of aligned bases
my $line_sum = $align_line =~ s/(\|)/$1/g;
my $align_ratio = ($line_sum / length($align_seq_1)) * 100;

#step4 : output the result of alignment
print  "----------------------------------------------------------------------\n";
print  "the score of the best align  is :\t$max_align_score\n";
printf "the ratio of aligned  bases  is :\t%.2f%%\n", $align_ratio;
print  "the result of aligned seq 1  is :\t$align_seq_1\n";
print  "                                \t$align_line\n" ;
print  "the result of aligned seq 2  is :\t$align_seq_2\n";
printf "the run time of SW algorithm is :\t%.4fs\n", $time2 - $time1;
print  "----------------------------------------------------------------------\n";

#======================================================================
#end main program
#======================================================================

#read the fasta file
sub GetFastaInfo {
	my $file_name = shift;
    open my $in_fa,'<',$file_name or die "can't open $file_name:$!";
	my %id_seq;
	while (<$in_fa>) {
    	chomp;
    	next unless $_;
    	$id_seq{$1} .= uc unless /\A>(.*)/;
	}
	close $in_fa;
    die "Error: the fasta file may be incorrect!\n" if keys %id_seq != 2;
	return \%id_seq;
}

#local alignment by using Smith-Waterman algorithm
sub SmithWaterman {
    #get the options
    my ($seq1, $seq2) = @_;

    #get the dimension of matrix
    my $row_num = (length $seq1) + 2;
    my $col_num = (length $seq2) + 2;

    #create 2 matrices, one for score and one for arrow
    my (@matrix_score, @matrix_arrow);
    push @matrix_score, [ (0) x $col_num ] foreach 1 .. $row_num;
    push @matrix_arrow, [ (0) x $col_num ] foreach 1 .. $row_num;

    #step1----------------------------------------------------------------------
    #for score and arrow matrix, fill the seq of row 1 and col 1, and fill the score/arrow of row 2 and col 2
    #note : in arrow matrix, -1 mean left, 1 mean up, 0 mean upper left, and 2 mean empty
    foreach my $i (0 , 1) {
        foreach my $j (2 .. $#{$matrix_score[$i]}) {
            $matrix_score[$i][$j] = $i == 0 ? substr $seq2, $j - 2, 1 : 0;
            $matrix_arrow[$i][$j] = $i == 0 ? substr $seq2, $j - 2, 1 : 2;
        }
    }
    foreach my $j (0 , 1) {
        foreach my $i (2 .. $#matrix_score) {
            $matrix_score[$i][$j] = $j == 0 ? substr $seq1, $i - 2, 1 : 0;
            $matrix_arrow[$i][$j] = $j == 0 ? substr $seq1, $i - 2, 1 : 2 ;
        }
    }

    #step2----------------------------------------------------------------------
    #for score and arrow matrix, fill the other cells, respectively
    foreach my $i (2 .. $#matrix_score) {
        foreach my $j (2 .. $#{$matrix_score[$i]}) {
            my $road_1_score = $matrix_score[$i][$j-1] + $space; #reach each cell from the left
            my $road_2_score = $matrix_score[$i-1][$j] + $space; #reach each cell from the up
            my $road_3_score = $matrix_score[$i][0] eq $matrix_score[0][$j] ?
                               $matrix_score[$i-1][$j-1] + $match : $matrix_score[$i-1][$j-1] + $unmatch; #reach each cell from the upper left
            if ($road_1_score > $road_2_score) {
                $matrix_score[$i][$j] = $road_1_score > $road_3_score ? $road_1_score : $road_3_score;
            }
            else{
                $matrix_score[$i][$j] = $road_2_score > $road_3_score ? $road_2_score : $road_3_score;
            }
            #if there are multiple paths for the maximum score, road 3 will be used first, road 2 will be used second, and road 1 will be used last
            $matrix_arrow[$i][$j] = $matrix_score[$i][$j] == $road_3_score ? 0 : ($matrix_score[$i][$j] == $road_2_score ? 1 : -1);
            #if the cell of the score matrix is negative, replace it with 0, and the corresponding arrow matrix cell is empty
            ($matrix_score[$i][$j], $matrix_arrow[$i][$j]) = (0, 2) if $matrix_score[$i][$j] < 0;
        }
    }

    #step3----------------------------------------------------------------------
    #get a maximum alignment score
    my ($row, $col) = (1, 1);
    foreach my $i (2 .. $#matrix_score) {
        foreach my $j (2 .. $#{$matrix_score[$i]}) {
            ($row, $col) = ($i, $j) if $matrix_score[$i][$j] > $matrix_score[$row][$col];
        }
    }
    my $max_align_score = $matrix_score[$row][$col];
    #get two aligned sequences
    my ($align_seq_1, $align_seq_2);
    my ($i, $j) = ($row, $col);
    while ($matrix_score[$i][$j]) {
        if ($matrix_arrow[$i][$j] == -1) {
            $align_seq_1 .= '-';
            $align_seq_2 .= $matrix_arrow[0][$j];
            $j--;
        }
        elsif ($matrix_arrow[$i][$j] == 1) {
            $align_seq_1 .= $matrix_arrow[$i][0];
            $align_seq_2 .= '-';
            $i--;
        }
        else{
            $align_seq_1 .= $matrix_arrow[$i][0];
            $align_seq_2 .= $matrix_arrow[0][$j];
            $i--;
            $j--;
        }
    }
    $align_seq_1 = reverse $align_seq_1;
    $align_seq_2 = reverse $align_seq_2;

    #return the result
    return $max_align_score, $align_seq_1, $align_seq_2;
}
