#!/usr/bin/env perl

## This file is part of YOC.
##
## YOC is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## YOC is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with YOC.  If not, see <http://www.gnu.org/licenses/>

# author:Raluca Uricaru

# description: based on a chaining result file .chn (coordinates of fragments taken in the chain) and a YASS file, this script extracts the pieces of YASS alignments corresponding to the fragments taken in the chain and outputs them together with the fragment coordinates (in a YASS like format)

# input: <1st_fasta> <2nd_fasta> <result_folder> <chn_file> <yass_output_folder>

# output: <chain_yass_file>

my $seq_file1 = $ARGV[0];
my $seq_file2 = $ARGV[1];
my $results_dir = $ARGV[2];
my $output_dir = $ARGV[2];
my $overlapChaining = $ARGV[3];
my $yass_output_1 = $ARGV[4];
my $yass_output_2 = $ARGV[4];

my $ok_yass = 1;
my $InverseFile = "Inversed_genomes.txt";

#circular shiftings if they exist
my $leftShift = 0;
my $rightShift = 0;

print $overlapChaining."\n";
compose_output($overlapChaining, $results_dir);

sub compose_output {
    my $chain_file = shift;
    my $chain_folder = shift;
    my @tokens = split('_',$chain_file);
    
    my $seq1 = @tokens[0];
    my $seq2 = @tokens[2];
    my $eval = @tokens[4];
    my @small_tokens = split('\.',$eval);
    $eval = @small_tokens[0];
    if ( @small_tokens[1]=='1')
    {
        $eval = "0.1";
    }
    
    $leftShift = 0;
    $rightShift = 0;
    
    my $output_file = $seq1 . "_GR_" . $seq2 . "_GR_" . $eval . ".chn_yass";
    
    my $yass_file1 = $seq1 . "_GR_" . $seq2 . "_GR_" . $eval .".yass1";
    my $yass_file2 = $seq1 . "_GR_" . $seq2 . "_GR_" . $eval .".yass";
    #print $yass_file2 . "\n";
    
    #chain result file
    open(FILE, $chain_folder."/".$chain_file) || die "error";
    #print $overlapChaining . "/" . $chain_file . "\n";
    
    #chain result file with idp
    #if ( -e $output_dir . "/" . $output_file)
    #{
    #	close(FILE) || die "error";
    #	return;
    #}
    open(Result_FILE, ">" . $output_dir . "/" . $output_file);
    #open(Covs_FILE, ">>" . "covs/cov_output_YOC.txt");
    
    
    my $i = 0;
    my $nr = 0;
    my $nr1 = 0;
    my $alignment_text = "";
    my $idp = 0;
    my $alignment_len = 0;
    my $matches = 0;
    my $mismatches = 0;
    my $gaps_length = 0;
    my $total_alignment = 0;
    my $total_mp = 0;
    my $cov1 = 0;
    my $cov2 = 0;
    my @prev_segment;
    my $nb_fragments = 0;
    
    while(<FILE>)
    {
        if (index($_,"total time") >= 0 )
        {
            last;
        }
        if ( $i == 0 )
        {
            @tokens = split(/ /,$_);
            
            if ( @tokens>2 )
            {
                print "circular\n";
                $leftShift = @tokens[3];
                $rightShift = @tokens[4];
            }
            
        }
        if ( $i > 2 ) {
            @tokens = split(/[\[\],]/,$_);
            my @current_segment = ( @tokens[1] , @tokens[2], @tokens[4], @tokens[5]);
            
            if ( @tokens[1] == 0 && @tokens[2] == 0 && @tokens[4] == 0 && @tokens[5] == 0 )
            {
                print Result_FILE "# 0 0 0 0 0\n";
                print Result_FILE $_;
                if ( $ok_yass == 1 )
                {
                    print Result_FILE "\n\n";
                }
            }
            else
            {
                #print "fragment " . ($i-3) . "\n";
                
                $idp = 0;
                $alignment_len = 0;
                $matches = 0;
                $mismatches = 0;
                $gaps_length = 0;
                $alignment_text = "";
                $overlap1 = 0;
                $overlap2 = 0;
                $nb_fragments++;
                
                obtain_fragment_with_idp(\@current_segment, $yass_file2, \$nr, \$idp, \$alignment_len,
                \$matches, \$mismatches, \$gaps_length);
                print Result_FILE "# " . $alignment_len . " " . $matches . " " . $mismatches . " " . $gaps_length . " " . $idp ."\n";
                print Result_FILE "[". @current_segment[0] . "," . @current_segment[1] . "] [" . @current_segment[2] . "," . @current_segment[3] . "] " .$idp . "\n";
                
                
                
                if ( $ok_yass == 1 )
                {
                    obtain_fragment_with_idp_align(\@current_segment, $yass_file1, \$nr1, \$alignment_text);
                    print Result_FILE $alignment_text;
                }
                
                $total_mp += $matches;
                $total_alignment += $alignment_len;
                
                $cov1 += @current_segment[1] - @current_segment[0] + 1;
                if ( @current_segment[2] < @current_segment[3] )
                {
                    $cov2 += @current_segment[3]  - @current_segment[2]  + 1;
                }
                else
                {
                    $cov2 += @current_segment[2] - @current_segment[3] + 1;
                }
                
                if ( $i>3 )
                {
                    if (@prev_segment[0] <= @current_segment[1] )
                    {
                        $overlap1 = @current_segment[1] - @prev_segment[0] +1;
                        
                    }
                    
                    if (@current_segment[2] < @current_segment[3] )
                    {
                        $x = @current_segment[2];
                        $y = @current_segment[3];
                    }
                    else
                    {
                        $y = @current_segment[2];
                        $x = @current_segment[3];
                    }
                    
                    
                    if ( @prev_segment[2] < @prev_segment[3] )
                    {
                        $z = @prev_segment[2];
                        $t = @prev_segment[3];
                    }
                    else
                    {
                        $t = @prev_segment[2];
                        $z = @prev_segment[3];
                    }
                    
                    if ($z <= $y && !($x > $t && $y > $t))
                    {
                        $overlap2 = $y - $z +1;
                    }
                    else
                    {
                        if ( $y > $t && $x <= $t )
                        {
                            $overlap2 = $t - $x +1;
                        }
                    }
                    $cov1 -= $overlap1;
                    $cov2 -= $overlap2;
                    $total_mp -= int(($overlap1+$overlap2)/2);
                    
                }
                
                
                @prev_segment[0] = @current_segment[0];
                @prev_segment[1] = @current_segment[1];
                @prev_segment[2] = @current_segment[2];
                @prev_segment[3] = @current_segment[3];
            }
            
            
        }
        else
        {
            if ( $i == 1 )
            {
                #not used in this version
                #@seq_files = search_sequence_files($seq1."_GR",$seq2."_GR", $eval);
                $size1 = `awk -f genomeLength.awk $seq_file1`;
                $size2 = `awk -f genomeLength.awk $seq_file2`;
                print Result_FILE "# ".$seq1 . "_GR " . $seq2 . "_GR " .$eval ."\n";
                print Result_FILE "# ".$seq_file1 . " " . $size1 . "# " . $seq_file2. " " .$size2;
                print Result_FILE "# alignment_length matches mismatches gaps_length %id\n\n";
            }
            elsif ( $i == 2 )
            {
                print Result_FILE "# 0 0 0 0 0\n";
            }
            if ( $i != 1 )
            {
                print Result_FILE $_;
            }
            if ( $i == 2 && $ok_yass == 1 )
            {
                print Result_FILE "\n\n";
            }
        }
        
        $i++;
    }
    #debug
    #if ( $nr != $i-4 || ($ok_yass == 1 && $nr1 != $i-4) )
    #{
    #print "\n" . "ERROR1 " . $nr . " " . $nr1 . "  ". ($i-3) . " " .$hchain_file . "\n";
    #exit 1;
    #}
    
    $cov1 = $cov1-1;
    $cov2 = $cov2-1;
    $total_mp = $total_mp-1;
    #print Covs_FILE $output_file." ".int($size1)." ".int($size2)." ".$cov1." ".$cov2." ".$total_mp." ".$total_alignment." ".$nb_fragments."\n";
    
    close(FILE) || die "error close file";
    close(Result_FILE) || die "error close file";
    #close(Covs_FILE) || die "error";
}

sub obtain_fragment_with_idp {
    my $current_segment_ptr = shift;
    my @current_segment = @{$current_segment_ptr};
    my $yass_file2 = shift;
    my $nr_ptr = shift;
    my $idp_ptr = shift;
    my $alignment_len = shift;
    my $matches = shift;
    my $mismatches = shift;
    my $gaps_length = shift;
    
    open(FILE2, $yass_output_2 . "/" . $yass_file2) || die "error";
    
    my $gaps_long_seq;
    my $len1;
    my $len2;
    my $gaps_commun;
    
    
    my $i = 0;
    my $ok = 0;
    
    while(<FILE2>)
    {
        if ( $i > 0 ) {
            my @tokens = split(/\t/,$_);
            #print @tokens[6] . " " . @tokens[7] . " " . @tokens[8] . " " . @tokens[9];
            
            if ( $leftShift==0 && $rightShift==0 )
            {
                if (( @current_segment[0] ==  @tokens[6] && @current_segment[1] ==  @tokens[7] )
                    &&
                    ( @current_segment[2] ==  @tokens[8] && @current_segment[3] ==  @tokens[9]
                    ||
                    @current_segment[2] ==  @tokens[9] && @current_segment[3] ==  @tokens[8] ))
                {
                    ${$gaps_length} = @tokens[5];
                    $len1 = abs(@tokens[7] - @tokens[6])+1;
                    $len2 = abs(@tokens[9] - @tokens[8])+1;
                    ${$mismatches} = @tokens[4];
                    
                    ${$nr_ptr}++;
                    $ok = 1;
                    
                    @{$current_segment_ptr}[0] = @tokens[6];
                    @{$current_segment_ptr}[1] = @tokens[7];
                    @{$current_segment_ptr}[2] = @tokens[8];
                    @{$current_segment_ptr}[3] = @tokens[9];
                    #print @tokens[6] . " " . @tokens[7] . " " . @tokens[8] . " " . @tokens[9];
                    last;
                }
            }
            else
            {
                #circular case
                #print "Circular case";
                if (( @current_segment[0] ==  @tokens[6] && @current_segment[1] ==  @tokens[7] )
                    &&
                    ( @current_segment[2] ==  @tokens[8]+$rightShift && @current_segment[3] ==  @tokens[9]+$rightShift
                    ||
                    @current_segment[2] ==  @tokens[9]+$rightShift && @current_segment[3] ==  @tokens[8]+$rightShift
                    ||
                    @current_segment[2] ==  @tokens[8]-$leftShift && @current_segment[3] ==  @tokens[9]-$leftShift
                    ||
                    @current_segment[2] ==  @tokens[9]-$leftShift && @current_segment[3] ==  @tokens[8]-$leftShift
                    ))
                {
                    ${$gaps_length} = @tokens[5];
                    $len1 = abs(@tokens[7] - @tokens[6])+1;
                    $len2 = abs(@tokens[9] - @tokens[8])+1;
                    ${$mismatches} = @tokens[4];
                    
                    ${$nr_ptr}++;
                    $ok = 1;
                    
                    @{$current_segment_ptr}[0] = @tokens[6];
                    @{$current_segment_ptr}[1] = @tokens[7];
                    @{$current_segment_ptr}[2] = @tokens[8];
                    @{$current_segment_ptr}[3] = @tokens[9];
                    
                    last;
                }
            }
        }
        
        $i++;
    }
    
    if ( $ok == 1)
    {
        
        #in order to compute the correct %id, as Yass1,2 are wrong on this computation
        if ( $len1 > $len2 )
        {
            $gaps_long_seq = $len1 - $len2;
        }
        else
        {
            $gaps_long_seq = $len2 - $len1;
        }
        $gaps_commun = (${$gaps_length} - $gaps_long_seq)/2;
        if ( $len1 > $len2 )
        {
            ${$matches} = $len1 - ${$mismatches} - $gaps_commun - $gaps_long_seq;
            ${$alignment_len} = ${$matches} + ${$mismatches} + ${$gaps_length};
        }
        else
        {
            ${$matches} = $len2 - ${$mismatches} - $gaps_commun - $gaps_long_seq;
            ${$alignment_len} = ${$matches} + ${$mismatches} + ${$gaps_length};		
        }	
        
        #print $matches . " " . $alignment_length;
        ${$idp_ptr} = int((${$matches}/${$alignment_len})*10000)/100;
        #print "\n".${$idp_ptr}."\n";
    }
    else
    {
        #print "\n" . "ERROR2 " . $chain_file . "\n";
        #print @current_segment[0] . " " . @current_segment[1];
        #exit 1;
    }
    
    close(FILE2) || die "error close file";
}

sub obtain_fragment_with_idp_align {
    my $current_segment_ptr = shift;
    my @current_segment = @{$current_segment_ptr};
    my $yass_file1 = shift;
    my $nr1_ptr = shift;
    my $alignment_text = shift;
    
    
    open(FILE1, $yass_output_1 . "/" . $yass_file1) || die "error open file";
    
    my $i = 0;
    my $ok = 0;
    my $idx = 0;
    
    while(<FILE1>)
    {
        
        if ( substr($_, 0, 2) eq "*(" )
        {
            if ( $ok == 1 )
            {
                last;
            }
            
            my @tokens = split(/[\s-\(\)\*]/,$_);
            
            if ( (@current_segment[0] ==  @tokens[2] && @current_segment[1] ==  @tokens[3] ||
                @current_segment[0] ==  @tokens[3] && @current_segment[1] ==  @tokens[2]) &&
                (@current_segment[2] ==  @tokens[5] && @current_segment[3] ==  @tokens[6] ||
                @current_segment[2] ==  @tokens[6] && @current_segment[3] ==  @tokens[5]) )
            {
                #print join( "--", @tokens);
                ${$nr1_ptr}++;
                $ok = 1;
                $idx = 0;			
            }
        }
        elsif ( substr($_, 0, 1) ne "*" && $ok==1 && $_ ne "\n" )
        {
            $idx++;
            if ( $idx%5 == 0 )
            {
                ${$alignment_text} .= $_. "\n";
            }
            else
            {
                ${$alignment_text} .= $_;
            }
        }
        $i++;
    }
    #debug
    #if ( $ok == 0 )
    #{
    #print "\n" . "ERROR3 " .  @current_segment[0]. " " . @current_segment[1] . " " .$chain_file . "\n";
    #exit 1;
    #}
    close(FILE1) || die "error";
}


sub search_sequence_files {
    my $seq1 = shift;
    my $seq2 = shift;
    my $eval = shift;
    
    my @seq_files;
    
    my @seq_files_names;
    
    @seq_files_names[0] = "";
    @seq_files_names[1] = "";
    
    #print "We search for ".$seq1." & ".$seq2;
    
    opendir(DIR, $seq_dir);
    my @files = readdir(DIR);
    foreach $file ( @files ) {
        
        if ( index($file,"_GR") == length($file)-7 &&  $file =~ /^$seq1(.)*\.f.a$/ ) {
            @seq_files[0] = $abs_file . "/" . $file;
            @seq_files_names[0] = $file;
        }
        elsif ( index($file,"_GR") == length($file)-7 && $file =~ /^$seq2(.)*\.f.a$/ ) {
            @seq_files[1] = $abs_file . "/" . $file;
            @seq_files_names[1] = $file;
        }
        
        
        if ( @seq_files_names[0] ne "" && @seq_files_names[1] ne "" )
        {		
            my $text;
            $text = substr(@seq_files_names[0], 0, length(@seq_files_names[0])-4);
            $text = $text . "_" . substr(@seq_files_names[1], 0, length(@seq_files_names[1])-4);
            
            # look in Inversed_genomes.txt if there's a line with $text
            # take @seq_files[0] = @seq_files[0]"_"
            if ( system("grep $text $InverseFile") == 0 )
            {
                @seq_files[0] = @seq_files[0]."_";
            }
        }
        
        if ( @seq_files_names[0] ne "" && @seq_files_names[1] ne "" )
        {
            last;
        }
        
    }
    closedir(DIR);
    
    return @seq_files;
}


