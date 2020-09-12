#!/usr/bin/perl

# Algorithm for predicting amino acid changes

# USAGE:
# perl Algorithm_predicting_full_AA_change.pl VCF_input_file_in_tab_delimited_format.txt.append kgXref.txt hg19_CDSIntronWithSign.txt ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.txt.out.txt

my ($snpfile,$convertfile,$gene_outfile,$knowngenefile) = @ARGV;
#the Codon and AA conversion hash:
my %aa_hash=(ATG=>'M',TGG=>'W',TTT=>'F',TTC=>'F',TAT=>'Y',TAC=>'Y',TGT=>'C',TGC=>'C',
  CAT=>'H',CAC=>'H',CAA=>'Q',CAG=>'Q',AAT=>'N',AAC=>'N',AAA=>'K',AAG=>'K',GAT=>'D',
  GAC=>'D',GAA=>'E',GAG=>'E',ATT=>'I',ATC=>'I',ATA=>'I',CCT=>'P',CCC=>'P',CCA=>'P',
  CCG=>'P',ACT=>'T',ACC=>'T',ACA=>'T',ACG=>'T',GTT=>'V',GTC=>'V',GTA=>'V',GTG=>'V',
  GCT=>'A',GCC=>'A',GCA=>'A',GCG=>'A',GGT=>'G',GGC=>'G',GGA=>'G',GGG=>'G',TCT=>'S',
  TCC=>'S',TCA=>'S',TCG=>'S',AGT=>'S',AGC=>'S',TTA=>'L',TTG=>'L',CTT=>'L',CTC=>'L',
  CTA=>'L',CTG=>'L',CGT=>'R',CGC=>'R',CGA=>'R',CGG=>'R',AGA=>'R',AGG=>'R',TAG=>'*',
  TAA=>'*',TGA=>'*',);
#my $gene_outfile = "$genefile.out";

#Hash kgXref.txt UCSC_ID->Gene_Symbol information
my %geneHash=();
open INFILEA, "$convertfile";
while(my $lineA = <INFILEA>)
{
  chomp($lineA);
  my @lineA_array = split('\t', $lineA);
  $geneHash{$lineA_array[0]} = $lineA_array[4];
}
close(INFILEA);

#Hash kgXref.txt UCSC_ID->strand information
my %strandHash=();
open INFILEB, "$knowngenefile";
while(my $lineB = <INFILEB>)
{
  chomp($lineB);
  my @lineB_array = split('\t', $lineB);
  $strandHash{$lineB_array[0]} = $lineB_array[2];
}
close(INFILEB);

my $SNP_output = "$snpfile.out.txt";  #  "$snpfile.snp.out.txt";
#my $indel_output = "$snpfile.indel.out.txt";
open OUTFILE_SNP, ">$SNP_output";
#open OUTFILE_INDEL, ">$indel_output";

#print header info
print OUTFILE_SNP "Sample","\t","Chromosome","\t","Variant Position","\t","Gene Symbol","\t","UCSC ID","\t","Strand","\t","AA Position of Mutation (for CDSHIT)","\t","Variant Type","\t","Amino Acid Ref (Codon) -> AA SNP (Codon)","\t","Variant Class","\t","Ref AA chain","\t","Alt AA chain","\t","Hit Type","\t","Known dbSNP","\t","Ref nt","\t","Alt nt","\t","Quality","\t","Depth","\t","Allele Freq","\t","Read Categories","\t","Info","\n";

#read the SNP file
open INFILE0, "$snpfile";
while(my $line0 = <INFILE0>)
{
  chomp($line0);
  my $SNP_file_flag = 1;  #assume it is a SNP; 1 means true
  my $protein_flag = -1;
  my @line0_array = split(/\t/, $line0);
  my $snp_chromosome = $line0_array[0];
  my $snp_location = $line0_array[1];
  my $ref_char = $line0_array[3];  #assume the ref has only one case
  $ref_char =~ tr/[a-z]/[A-Z]/;  #force everything to upper case
  my @snp_chars;
  my $snp_char;
  my $case_counter;
  if($line0_array[4] =~ m/","/)
  {
    @snp_chars = split(/","/, $line0_array[4]);
    $case_counter = $#snp_chars + 1;
    for(my $k=0; $k<($#snp_chars + 1); $k++)
    {
      if(length($snp_chars[$k]) > 1)
      { #this is not a single SNP
        $SNP_file_flag = -1;
        last;
      }
    }
  }
  else
  {
    $snp_char = $line0_array[4];
    $snp_char =~ tr/[a-z]/[A-Z]/;
    $case_counter = 1;
    if(length($snp_char) > 1)
    { #this is not a single SNP
      $SNP_file_flag = -1;
    }
  }

  #get depth and ReadCategory info
  my @info_array = split(/\;/, $line0_array[7]);
  my @depth_array;
  my @alle_freq;
  my @read_category;
  my $hit_type;
  my $UCSC_ID;
  if($line0_array[7] =~ m/VDB/)  #for new samtools version 0.1.18
  {
    if(($#line0_array + 1) == 14) #three samples run from samtools
    {
      $hit_type = $line0_array[12];
      $UCSC_ID = $line0_array[13];
      print "Three samples & new samtools...\n";
    }
    elsif(($#line0_array + 1) == 13) #two samples run from samtools
    {
      $hit_type = $line0_array[11];
      $UCSC_ID = $line0_array[12];
      print "Two samples & new samtools...\n";
    }
    elsif(($#line0_array + 1) == 12) #one sample run from samtools
    {
      $hit_type = $line0_array[10];
      $UCSC_ID = $line0_array[11];
      print "One sample & new samtools...\n";
    }
  }
  else
  {
    $hit_type = $line0_array[10];
    $UCSC_ID = $line0_array[11];
    print "One sample & old samtools...\n";
  }

  my $UCSC_ID_flag = -1;
  my $strand;
  my $codon_change_string;
  if($hit_type eq "CDSHIT")
  {
    print "\n============================\n",$line0,"\t",$geneHash{$UCSC_ID},"\n";
    print "The possible ALT cases for this line are $case_counter\n";
    if((length($ref_char) == 1) && ($SNP_file_flag eq 1))  # if both ref and SNP calls are single "SNPs"
    {
      $protein_flag = 1;
      print OUTFILE_SNP substr($snpfile,0,index($snpfile,".")),"\t",substr($line0_array[0],3),"\t",$line0_array[1],"\t",$geneHash{$UCSC_ID},"\t",$UCSC_ID,"\t";
    }
    else #at least ref OR SNP calls are NOT a single "SNP"
    {
      print OUTFILE_SNP substr($snpfile,0,index($snpfile,".")),"\t",substr($line0_array[0],3),"\t",$line0_array[1],"\t",$geneHash{$UCSC_ID},"\t",$UCSC_ID,"\t";
    }
    #to check how many ALT cases need to be gone through
    my $looper;
    if($case_counter == 1)  #if only one ALT case
    {
      $looper = $case_counter;
      #my $strand;
      my $CDS_start;
      my $CDS_end;
      my $line_flag = -1;
      open(INFILE2, "$gene_outfile");
      while(my $line2 = <INFILE2>)
      {
        chomp($line2);
        if($line2 =~ m/">"/)
        {
          my @line2_array = split(/ /, $line2);
          my $UCSC_ID_check = substr($line2_array[0],16);
          if($UCSC_ID eq $UCSC_ID_check)  #if UCSC_ID is found in the CDSIntron file
          {
            my @line22_array = split(/:/, $line2_array[1]);
            my @line22_arrayChr = split('=', $line22_array[0]);
            my @line222_array = split('-', $line22_array[1]);
            $CDS_start = $line222_array[0];
            $CDS_end = $line222_array[1];
            if ($snp_chromosome eq $line22_arrayChr[1])  #also their chromosomes are the same
            {
              if(($snp_location < $CDS_start) || ($snp_location > $CDS_end))
              {
                print "SNP location is outside of CDS region!!! No checking occurs!\n";
                #exit(1);
              }
              $strand = substr($line2_array[4],-1,1);
              $line_flag = 1;
              $UCSC_ID_flag = 1;  #newly inserted condition
              next;
            }
          }
        }
        if($line_flag == 1)  #found CDS line
        {
          my $CDS_line = $line2;
          my $line2_noIntron = remove_intron($line2);
          my $CDS_line_noIntron;
          #print "\nThe length of CDS with intron for $UCSC_ID ($gene_name) is: ", length($CDS_line), "\n";
          my $coordinate;
          my $before_SNP_string;
          my $c_lower;
          if($strand eq "+")
          {
            #print "CDS start location is $CDS_start\n";
            $coordinate = $snp_location - $CDS_start;
            $before_SNP_string = substr($line2, 0, $coordinate);  #count how many lower case nucleotides (NON-CDS or intron nucleotide) before SNP location and adjust its coordinate
            $c_lower = $before_SNP_string =~ tr/a-z//;
            #print $snp_location, "\t", $CDS_start, "\t",$CDS_end, "\n";
            #print "Before replacement: ",substr($CDS_line,$coordinate,1),"\n";
            eval {substr($CDS_line,$coordinate,length($ref_char),$snp_char) || die "This is not the right CDS length: $!";};
            if($@)
            {
              #if($@ =~ m/substr outside of string/)
              #{
              print "There is an error for checking $UCSC_ID ($gene_name): $@";
              #}
            }
            #print "After replacement: ",substr($CDS_line,$coordinate,1),"\n";
            else
            { #process new CDS string
              $CDS_line_noIntron = remove_intron($CDS_line);
              #============== The last two parameters are just for print OUTFILE_SNP purpose
              protein_translation($line2_noIntron,$CDS_line_noIntron,$coordinate-$c_lower,$protein_flag,$strand,$case_counter,$looper);
            }
            $line_flag = -1;
            last;
          }
          else
          {
            #print "CDS start location is $CDS_end\n";
            $coordinate = $CDS_end - $snp_location;
            my $replaced_snp_char = revdnacomp($snp_char);
            $replaced_snp_char =~ tr/[a-z]/[A-Z]/;
            $before_SNP_string = substr($line2, 0, $coordinate);
            $c_lower = $before_SNP_string =~ tr/a-z//;
            #print $replaced_snp_char, "\n";
            #print $snp_location, "\t", $CDS_start, "\t",$CDS_end, "\n";
            #print "Before replacement: ",substr($CDS_line,$coordinate,1),"\n";
            ## need to go further back to the length of ref_char to replace it since it it on the reverse strand
            eval {substr($CDS_line,$coordinate-length($ref_char)+1,length($ref_char),$replaced_snp_char) || die "This is not the right CDS length: $!";};
            if($@)
            {
              print "There is an error for checking $UCSC_ID ($gene_name): $@";
            }
            #print "After replacement: ", substr($CDS_line, $coordinate, 1), "\n";
            else
            { #process new CDS string
            #  print $UCSC_ID, "\n";
            #  if ($UCSC_ID eq 'uc001wja.2')
            #  {
            #    print "New CDS line\n";
            #    print $CDS_line, "\n";
            #    print revdnacomp($ref_char), "\n";
            #    print $line2_noIntron, "\n";
            #    print "============$replaced_snp_char==============\nbbbbbbbbbbbbbb\n";
              $CDS_line_noIntron = remove_intron($CDS_line);
            #    print $CDS_line_noIntron, "\n";
            #  }
              #process CDS sequece with intron removed sequence to identify SNP new coordinate for use later
              #print $strand,"\t",$snp_location, "\t",$snp_location_new,"\t",length($CDS_line_noIntron),"\n";
              #print $line0,"\t",$geneHash{$UCSC_ID},"\n";
              protein_translation($line2_noIntron,$CDS_line_noIntron,$coordinate-length($ref_char)+1-$c_lower,$protein_flag,$strand,$case_counter,$looper);
            }
            $line_flag = -1;
            last;
          }
        }
      }
      close(INFILE2);
    }#if
    else  #there are more than one possible ALT cases or $#snp_chars > 0
    {
      $looper = $case_counter;
      for(my $case=0; $case<($#snp_chars + 1); $case++)  # for each possible ALT case
      {
        $snp_chars[$case] =~ tr/[a-z]/[A-Z]/;
        print "\n--------------------------------------\n";
        #my $strand;
        my $CDS_start;
        my $CDS_end;
        my $line_flag = -1;
        open(INFILE2, "$gene_outfile");
        while(my $line2 = <INFILE2>)
        {
          chomp($line2);
          if($line2 =~ m/">"/)
          {
            my @line2_array = split(/ /, $line2);
            my $UCSC_ID_check = substr($line2_array[0],16);
            if($UCSC_ID eq $UCSC_ID_check)  #if UCSC_ID is found in the CDSIntron file
            {
              my @line22_array = split(/:/, $line2_array[1]);
              my @line22_arrayChr = split('=', $line22_array[0]);
              my @line222_array = split('-', $line22_array[1]);
              $CDS_start = $line222_array[0];
              $CDS_end = $line222_array[1];
              if($snp_chromosome eq $line22_arrayChr[1])  #also their chromosomes are the same
              {
                if(($snp_location < $CDS_start) || ($snp_location > $CDS_end))
                {
                  print "SNP location is outside of CDS region!!! No checking occurs!\n";
                  #exit(1);
                }
                $strand = substr($line2_array[4],-1,1);
                $line_flag = 1;
                $UCSC_ID_flag = 1;  #newly inserted condition
                next;
              }
            }
          }
          if($line_flag == 1)  #found CDS line
          {
            my $CDS_line = $line2;
            my $line2_noIntron = remove_intron($line2);
            my $CDS_line_noIntron;
            #print "\nThe length of CDS with intron for $UCSC_ID ($gene_name) is: ", length($CDS_line), "\n";
            my $coordinate;
            my $before_SNP_string;
            my $c_lower;
            #replace original char with SNP
            if($strand eq "+")
            {
              #print "CDS start location is $CDS_start\n";
              $coordinate = $snp_location - $CDS_start;
              $before_SNP_string = substr($line2, 0, $coordinate);
              $c_lower = $before_SNP_string =~ tr/a-z//;
              #print $snp_location, "\t", $CDS_start, "\t",$CDS_end, "\n";
              #print "Before replacement: ", substr($CDS_line, $coordinate, 1), "\n";
              eval {substr($CDS_line,$coordinate,length($ref_char),$snp_chars[$case]) || die "This is not the right CDS length: $!";};
              if($@)
              {
                #if($@ =~ m/substr outside of string/)
                #{
                print "There is an error for checking $UCSC_ID ($gene_name): $@";
                #}
              }
              #print "After replacement: ",substr($CDS_line,$coordinate,1),"\n";
              else
              { #process new CDS string
                $CDS_line_noIntron = remove_intron($CDS_line);
                #process CDS sequece with intron removed sequence to identify SNP new coordinate for use later
                #print $strand, "\t", $snp_location, "\t", $snp_location_new, "\t", length($CDS_line_noIntron), "\n";
                #print $line0, "\t", $geneHash{$UCSC_ID}, "\n";
                $codon_change_string .= protein_translation($line2_noIntron,$CDS_line_noIntron,$coordinate-$c_lower,$protein_flag,$strand,$case_counter,$looper);
              }
              $line_flag = -1;
              last;
            }
            else
            {
              #print "CDS start location is $CDS_end\n";
              $coordinate = $CDS_end - $snp_location;
              my $replaced_snp_char = revdnacomp($snp_chars[$case]);
              $replaced_snp_char =~ tr/[a-z]/[A-Z]/;
              $before_SNP_string = substr($line2, 0, $coordinate);
              $c_lower = $before_SNP_string =~ tr/a-z//;
              #print $replaced_snp_char, "\n";
              #print $snp_location, "\t", $CDS_start, "\t",$CDS_end, "\n";
              #print "Before replacement: ", substr($CDS_line, $coordinate, 1), "\n";
              eval {substr($CDS_line,$coordinate-length($ref_char)+1,length($ref_char),$replaced_snp_char) || die "This is not the right CDS length: $!";};
              if($@)
              {
                print "There is an error for checking $UCSC_ID ($gene_name): $@";
              }
              #print "After replacement: ", substr($CDS_line, $coordinate, 1), "\n";
              else
              { #process new CDS string
                $CDS_line_noIntron = remove_intron($CDS_line);
                #process CDS sequece with intron removed sequence to identify SNP new coordinate for use later
                #print $strand, "\t", $snp_location, "\t", $snp_location_new, "\t", length($CDS_line_noIntron), "\n";
                #print $line0, "\t", $geneHash{$UCSC_ID}, "\n";
                $codon_change_string .= protein_translation($line2_noIntron,$CDS_line_noIntron,$coordinate-length($ref_char)+1-$c_lower,$protein_flag,$strand,$case_counter,$looper);
              }
              $line_flag = -1;
              last;
            }
          }
        }
        close(INFILE2);
        $looper--;
      }#for each ALT case
    }#else
    if((length($ref_char) == 1) && ($SNP_file_flag eq 1)) #single SNP case
    {
      @depth_array = split("=", $info_array[0]);
      if($line0_array[7] =~ m/VDB/)  #for new samtools version 0.1.18
      {
        @alle_freq = split("=", $info_array[2]);
        @read_category = split("=", $info_array[4]);
      }
      else #for old samtools version 0.1.12
      {
        @alle_freq = split("=", $info_array[1]);
        @read_category = split("=", $info_array[3]);
      }
      if(length($codon_change_string) == 0) #if there is no ALT case
      {
        print OUTFILE_SNP $hit_type,"\t",$line0_array[2],"\t",$line0_array[3],"\t",$line0_array[4],"\t",$line0_array[5],"\t",$depth_array[1],"\t",$alle_freq[1],"\t",$read_category[1],"\t",$line0_array[7],"\n";
      }
      else
      {
        print OUTFILE_SNP $codon_change_string,"\t",$hit_type,"\t",$line0_array[2],"\t",$line0_array[3],"\t",$line0_array[4],"\t",$line0_array[5],"\t",$depth_array[1],"\t",$alle_freq[1],"\t",$read_category[1],"\t",$line0_array[7],"\n";
      }
    }
    else #indel case
    {
      @depth_array = split("=", $info_array[1]);
      if($line0_array[7] =~ m/VDB/)  #for new samtools version 0.1.18
      {
        @alle_freq = split("=", $info_array[3]);
        @read_category = split("=", $info_array[5]);
      }
      else #for old samtools version 0.1.12
      {
        @alle_freq = split("=", $info_array[2]);
        @read_category = split("=", $info_array[4]);
      }
      print OUTFILE_SNP $strand,"\t","---","\t","INDEL","\t","---","\t","---","\t","---","\t","---","\t"; 
      print OUTFILE_SNP $hit_type,"\t",$line0_array[2],"\t",$line0_array[3],"\t",$line0_array[4],"\t",$line0_array[5],"\t",$depth_array[1],"\t",$alle_freq[1],"\t",$read_category[1],"\t",$line0_array[7],"\n";
    }
    if($UCSC_ID_flag == -1)
    {
      print "$UCSC_ID is not found!!!\n";
      #exit(1);
    }
  }#if CDSHIT
  #elsif(($hit_type eq "INTRONHIT") || ($hit_type eq "UPSTREAMHIT") || ($hit_type eq "DOWNSTREAMHIT"))
  #{
  #}
  else  #INTRON, UTR or UP-DOWNSTREAM cases
  {
    if((length($ref_char) == 1) && ($SNP_file_flag eq 1))
    {
      print OUTFILE_SNP substr($snpfile,0,index($snpfile,".")),"\t",substr($line0_array[0],3),"\t",$line0_array[1],"\t",$geneHash{$UCSC_ID},"\t",$UCSC_ID,"\t";
      @depth_array = split("=", $info_array[0]);
      if($line0_array[7] =~ m/VDB/)  #for new samtools version 0.1.18
      {
        @alle_freq = split("=", $info_array[2]);
        @read_category = split("=", $info_array[4]);
      }
      else #for old samtools version 0.1.12
      {
        @alle_freq = split("=", $info_array[1]);
        @read_category = split("=", $info_array[3]);
      }
      print OUTFILE_SNP $strandHash{$UCSC_ID},"\t","---","\t","SNP","\t","---","\t","---","\t","---","\t","---","\t";
      print OUTFILE_SNP $hit_type,"\t",$line0_array[2],"\t",$line0_array[3],"\t",$line0_array[4],"\t",$line0_array[5],"\t",$depth_array[1],"\t",$alle_freq[1],"\t",$read_category[1],"\t",$line0_array[7],"\n";
    }
    else #at least ref OR SNP calls are NOT a single "SNP"
    {
      print OUTFILE_SNP substr($snpfile,0,index($snpfile,".")),"\t",substr($line0_array[0],3),"\t",$line0_array[1],"\t",$geneHash{$UCSC_ID},"\t",$UCSC_ID,"\t";
      @depth_array = split("=", $info_array[1]);
      if($line0_array[7] =~ m/VDB/)  #for new samtools version 0.1.18
      {
        @alle_freq = split("=", $info_array[3]);
        @read_category = split("=", $info_array[5]);
      }
      else #for old samtools version 0.1.12
      {
        @alle_freq = split("=", $info_array[2]);
        @read_category = split("=", $info_array[4]);
      }
      print OUTFILE_SNP $strandHash{$UCSC_ID},"\t","---","\t","INDEL","\t","---","\t","---","\t","---","\t","---","\t";
      print OUTFILE_SNP $hit_type,"\t",$line0_array[2],"\t",$line0_array[3],"\t",$line0_array[4],"\t",$line0_array[5],"\t",$depth_array[1],"\t",$alle_freq[1],"\t",$read_category[1],"\t",$line0_array[7],"\n";
    }
  }
}
close(INFILE0);
close(OUTFILE_SNP);
#close(OUTFILE_INDEL);

#function which does the complement of nucleotide
#sub complement
#{
#  my ($str) = @_;
#  $str =~ tr/ACGT/TGCA/;
#  return $str;
#}

#### http://www.perlmonks.org/?node_id=197793
sub revdnacomp {
  # my $dna = @_;  # $dna gets the number of arguments in @_, since it's a scalar context!
  my $dna = shift;  # or  my $dna = shift @_; # ah, scalar context of scalar gives expected results.
  # my ($dna) = @_;  # would work, too
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

#function to remove intron using its lower case characterristics when download
sub remove_intron {
  my ($full_line) = @_;
  $full_line =~ s/a//g;
  $full_line =~ s/c//g;
  $full_line =~ s/g//g;
  $full_line =~ s/t//g;
  $full_line =~ s/n//g;
  return $full_line;
}

#function which does the AA translation
sub protein_translation {
  my($input_original_CDS_line,$input_CDS_line,$input_snp_location,$input_protein_flag,$input_strand,$input_case_counter,$input_looper)=@_;
  #print "Adjusted SNP location started at:", $input_snp_location, "\n";
  #print "original CDS length: ", length($input_original_CDS_line), "\n";
  #print "replaced CDS length: ", length($input_CDS_line), "\n";
  my $input_codon_change_string;
  my $truncate_orignal_CDS_line;
  my $truncate_orignal_AA_line;
  my $checked_codon;
  for(my $i=0; $i<length($input_original_CDS_line); $i+=3)  #translate CDS into AA sequence by moving 3 nt one time
  { #get the codon for original sequence
    my $original_codon = substr($input_original_CDS_line, $i, 3);
  #  if(($i-$input_snp_location) >= -3)  # if this is the starting mutated codon location three nc ahead
  #  {
      #print $i, "\n";
      $truncate_orignal_CDS_line .= $original_codon;
      $truncate_orignal_AA_line .= $aa_hash{$original_codon};
  #  }
    if((($input_snp_location-$i)>=0) && (($input_snp_location-$i)<=2) && ($input_protein_flag==1))  #if this is the codon position
    {
      if($input_case_counter == 1)
      {
        print OUTFILE_SNP $input_strand,"\t",$i/3+1,"\t","SNP","\t",$aa_hash{$original_codon},"(",$original_codon,")";
        $checked_codon = $aa_hash{$original_codon};
      }
      else
      { #if this is the first case of ALT
        if($input_looper == $input_case_counter)
        {
          print OUTFILE_SNP $input_strand,"\t",$i/3+1,"\t","SNP","\t",$aa_hash{$original_codon},"(",$original_codon,")";
          $checked_codon = $aa_hash{$original_codon};
        }
        else #if this is the other cases of ALT
        {
          print OUTFILE_SNP $aa_hash{$original_codon},"(",$original_codon,")";
          $checked_codon = $aa_hash{$original_codon};
        }
      }
    }
  }#for
  #print "\n";
  my $truncate_CDS_line;
  my $truncate_AA_line;
  for(my $j=0; $j<length($input_CDS_line); $j+=3)  #translate CDS into AA sequence by moving 3 nt one time
  { #get the codon for mutated sequence
    my $codon = substr($input_CDS_line, $j, 3);
  #  if(($j-$input_snp_location) >= -3)  #if this is the starting mutated codon location
  #  {
      #print $j, "\n";
      $truncate_CDS_line .= $codon;
      $truncate_AA_line .= $aa_hash{$codon};
  #  }
    if((($input_snp_location-$j)>=0) && (($input_snp_location-$j)<=2) && ($input_protein_flag==1))
    {
      $SNP_codon_mutant = $codon;
      $SNP_aa_mutant = $aa_hash{$codon};
      if($input_case_counter == 1)
      {
        print OUTFILE_SNP "->", $aa_hash{$codon}, "(", $codon, ")", "\t";
        if($checked_codon eq $aa_hash{$codon})
        {
          print OUTFILE_SNP "SYN", "\t";
        }
        else
        {
          if($aa_hash{$codon} eq "*")
          {
            print OUTFILE_SNP "NSN", "\t";
          }
          else
          {
            print OUTFILE_SNP "NSM", "\t";
          }
        }
      }
      else
      { #if this is the other case of ALT
        if($input_looper > 1)
        {
          print OUTFILE_SNP "->", $aa_hash{$codon}, "(", $codon, ")", ",";
          if($checked_codon eq $aa_hash{$codon})
          {
            $input_codon_change_string .= "SYN";
            $input_codon_change_string .= ",";
          }
          else
          {
            if($aa_hash{$codon} eq "*")
            {
              $input_codon_change_string .= "NSN";
              $input_codon_change_string .= ",";
            }
            else
            {
              $input_codon_change_string .= "NSM";
              $input_codon_change_string .= ",";
            }
          }
        }
        else #if this is the last case of ALT
        {
          print OUTFILE_SNP "->", $aa_hash{$codon}, "(", $codon, ")", "\t";
          if($checked_codon eq $aa_hash{$codon})
          {
            $input_codon_change_string .= "SYN";
          }
          else
          {
            if($aa_hash{$codon} eq "*")
            {
              $input_codon_change_string .= "NSN";
            }
            else
            {
              $input_codon_change_string .= "NSM";
            }
          }
          #print OUTFILE_SNP $input_codon_change_string, "\t";
        }
      }
    }
  }#for
  print $truncate_orignal_CDS_line, "\n";
  print $truncate_CDS_line, "\n";
  print $truncate_orignal_AA_line, "\n";
  print $truncate_AA_line, "\n";
  #newly added for printing full length of AA
  print OUTFILE_SNP $truncate_orignal_AA_line, "\t";
  print OUTFILE_SNP $truncate_AA_line, "\t";
  return $input_codon_change_string;
}#end sub
