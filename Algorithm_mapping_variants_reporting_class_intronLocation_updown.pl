#!/usr/bin/env perl

# Algorithm for mapping identified variants onto the genomic location and
# reporting the hit class

# USAGE:
# perl Algorithm_mapping_variants_reports_class_intronLocation.pl ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format.txt IntronExon_boundary_in_bp

use POSIX;

my $exonfile;
my $snpfile;
my $exonbuffer;

my $intronOption;

if(($#ARGV + 1)== 2 )
{
  ($exonfile,$snpfile) = @ARGV;
  print "The program assumes that you do NOT want to report how far the variant falls in the exon boundary.\n";
  $intronOption = -1;
}
elsif(($#ARGV + 1) == 3 )
{
  ($exonfile,$snpfile,$exonbuffer) = @ARGV;
  print "The program assumes that you DO want to report how far the variant falls in the exon boundary.\nOnly variants flanking its nearby exon with <= $exonbuffer bp is reported\n";
  $intronOption = 1;
}
else
{
  print "The input commands do not meet the requirement. Please see the README file and try again.\n";
  exit;
}

%chr_cds=();
%chr_intron=();
%chr_utr5=();
%chr_utr3=();
%chr_upstream=();
%chr_downstream=();
%gene_cds=();
%gene_intron=();
%gene_utr5=();
%gene_utr3=();
%gene_upstream=();
%gene_downstream=();

#For each feature, make a hash for all unique start as key and end as value
open(EXONFILECDSLINK, "$exonfile.cds_link_shrink");
while(<EXONFILECDSLINK>)
{
  chomp;
  @exonlinecds=split(/\t/,$_);
  $chr_cds{$exonlinecds[0]} = $exonlinecds[1];
}
close EXONFILECDSLINK;

open(EXONFILECDSGENE, "$exonfile.cds_gene");
while(<EXONFILECDSGENE>)
{
  chomp;
  @exonlinecdsgene=split(/\t/,$_);
  $gene_cds{$exonlinecdsgene[0]} = $exonlinecdsgene[1];
}
close EXONFILECDSGENE;

open(EXONFILEINTRONLINK, "$exonfile.intron_link_shrink");
while(<EXONFILEINTRONLINK>)
{
  chomp;
  @exonlineintron=split(/\t/,$_);
  $chr_intron{$exonlineintron[0]} = $exonlineintron[1];
}
close EXONFILEINTRONLINK;

open(EXONFILEINTRONGENE, "$exonfile.intron_gene");
while(<EXONFILEINTRONGENE>)
{
  chomp;
  @exonlineintrongene=split(/\t/,$_);
  $gene_intron{$exonlineintrongene[0]} = $exonlineintrongene[1];
}
close EXONFILEINTRONGENE;

open(EXONFILEUTR5LINK, "$exonfile.utr5_link_shrink");
while(<EXONFILEUTR5LINK>)
{
  chomp;
  @exonlineutr5=split(/\t/,$_);
  $chr_utr5{$exonlineutr5[0]} = $exonlineutr5[1];
}
close EXONFILEUTR5LINK;

open(EXONFILEUTR5GENE, "$exonfile.utr5_gene");
while(<EXONFILEUTR5GENE>)
{
  chomp;
  @exonlineutr5gene=split(/\t/,$_);
  $gene_utr5{$exonlineutr5gene[0]} = $exonlineutr5gene[1];
}
close EXONFILEUTR5GENE;

open(EXONFILEUTR3LINK, "$exonfile.utr3_link_shrink");
while(<EXONFILEUTR3LINK>)
{
  chomp;
  @exonlineutr3=split(/\t/,$_);
  $chr_utr3{$exonlineutr3[0]} = $exonlineutr3[1];
}
close EXONFILEUTR3LINK;

open(EXONFILEUTR3GENE, "$exonfile.utr3_gene");
while(<EXONFILEUTR3GENE>)
{
  chomp;
  @exonlineutr3gene=split(/\t/,$_);
  $gene_utr3{$exonlineutr3gene[0]} = $exonlineutr3gene[1];
}
close EXONFILEUTR3GENE;

open(EXONFILEUPSTREAMLINK, "$exonfile.upstream_link_shrink");
while(<EXONFILEUPSTREAMLINK>)
{
  chomp;
  @exonlineupstream=split(/\t/,$_);
  $chr_upstream{$exonlineupstream[0]} = $exonlineupstream[1];
}
close EXONFILEUPSTREAMLINK;

open(EXONFILEUPSTREAMGENE, "$exonfile.upstream_gene");
while(<EXONFILEUPSTREAMGENE>)
{
  chomp;
  @exonlineupstreamgene=split(/\t/,$_);
  $gene_upstream{$exonlineupstreamgene[0]} = $exonlineupstreamgene[1];
}
close EXONFILEUPSTREAMGENE;

open(EXONFILEDOWNSTREAMLINK, "$exonfile.downstream_link_shrink");
while(<EXONFILEDOWNSTREAMLINK>)
{
  chomp;
  @exonlinedownstream=split(/\t/,$_);
  $chr_downstream{$exonlinedownstream[0]} = $exonlinedownstream[1];
}
close EXONFILEDOWNSTREAMLINK;

open(EXONFILEDOWNSTREAMGENE, "$exonfile.downstream_gene");
while(<EXONFILEDOWNSTREAMGENE>)
{
  chomp;
  @exonlinedownstreamgene=split(/\t/,$_);
  $gene_downstream{$exonlinedownstreamgene[0]} = $exonlinedownstreamgene[1];
}
close EXONFILEDOWNSTREAMGENE;

#map vcf variant back onto genome location and report the hit class
open OUTFILE1, ">$snpfile.append";

open SNPFILE, "<$snpfile";
while(my $line = <SNPFILE>)
{
  my @tempArray_cds;
  my @tempArray_intron;
  my @tempArray_utr5;
  my @tempArray_utr3;
  my @tempArray_upstream;
  my @tempArray_downstream;

  chomp($line);
  if($line =~ m/^#/)
  {
    next;
  }
  else
  {
  @snpline = split(/\t/,$line);
  $snp_chr = $snpline[0];
  $snp_start = $snpline[1];
  open(EXONFILECDS, "$exonfile.cds");
  while (my $linecds = <EXONFILECDS>)
  {
    chomp($linecds);
    my @linecds_array = split(/\t/, $linecds);
    #given the same genome chromosome as variant's chrom
    if($snp_chr eq $linecds_array[0])
    { #only store all of start postions in the array (without chromosome information included) for binary search
      for ($j_cds=1; $j_cds<($#linecds_array + 1); $j_cds++)
      {
        push @tempArray_cds, $linecds_array[$j_cds];
      }
    }
  }
  close(EXONFILECDS);

      #cdsstart is 0-based, and cdsend is 1-based
  my @sortchromArray_cds = @tempArray_cds;
      $length_cds= scalar(@sortchromArray_cds);
      $cursor_cds=$length_cds/2;
      $localmin_cds=0;
      $localmax_cds=$length_cds;
      $traverse_cds=0;
  until (($localmax_cds-$localmin_cds)<=1) {
  #print "cursor=|$cursor|\n";
        if ($snp_start <= $sortchromArray_cds[$cursor_cds]) {
                        $localmax_cds=floor($cursor_cds);
                        &left_cds($cursor_cds, $localmin_cds);}
        if ($snp_start >= $sortchromArray_cds[$cursor_cds]) {
                        $localmin_cds=floor($cursor_cds);
                        &right_cds($cursor_cds, $localmax_cds);}

        $traverse_cds++;
                if ($traverse_cds >100) {
                        die "Excess Traverse: min=$localmin_cds\tmax=$localmax_cds\tcursor=$cursor_cds\tVPOS=$snp_start\tArrayCursor=$sortchromArray_cds[$cursor_cds]\t$sortchromArray_cds[$cursor_cds-1]\t$sortchromArray_cds[$cursor_cds+1]\n";
                       }
        } # end until

        if ( ($snp_start > $sortchromArray_cds[$cursor_cds]) && ($snp_start<=$chr_cds{"$snp_chr"._."$sortchromArray_cds[$cursor_cds]"})  && ($sortchromArray_cds[$cursor_cds] ne "NA") && ($chr_cds{"$snp_chr"._."$sortchromArray_cds[$cursor_cds]"} ne "NA")  )
        {
          print OUTFILE1 $line, "\t", "CDSHIT", "\t", $gene_cds{"$snp_chr"._."$sortchromArray_cds[$cursor_cds]"}, "\n";
          $cdsCount++;
        }

  open(EXONFILEINTRON, "$exonfile.intron");
  while (my $lineintron = <EXONFILEINTRON>)
  {
    chomp($lineintron);
    my @lineintron_array = split(/\t/, $lineintron);
    if($snp_chr eq $lineintron_array[0])
    {
      for ($j_intron=1; $j_intron<($#lineintron_array + 1); $j_intron++)
      {
        push @tempArray_intron, $lineintron_array[$j_intron];
      }
    }
  }
  close(EXONFILEINTRON);


      # inronstart is 1-based, add intronhit buffer number for parse_discrepancy_jan_inronSeparate_final.perl because intronend is 0-based
  my @sortchromArray_intron = @tempArray_intron;
      $length_intron= scalar(@sortchromArray_intron);
      $cursor_intron=$length_intron/2;
      $localmin_intron=0;
      $localmax_intron=$length_intron;
      $traverse_intron=0;
  until (($localmax_intron-$localmin_intron)<=1) {
        if ($snp_start <= $sortchromArray_intron[$cursor_intron]) {
                        $localmax_intron=floor($cursor_intron);
                        &left_intron($cursor_intron, $localmin_intron);}
        if ($snp_start >= $sortchromArray_intron[$cursor_intron]) {
                        $localmin_intron=floor($cursor_intron);
                        &right_intron($cursor_intron, $localmax_intron);}

        $traverse_intron++;
                if ($traverse_intron >100) {
                        die "Excess Traverse: min=$localmin_intron\tmax=$localmax_intron\tcursor=$cursor_intron\tVPOS=$snp_start\tArrayCursor=$sortchromArray_intron[$cursor_intron]\t$sortchromArray_intron[$cursor_intron-1]\t$sortchromArray_intron[$cursor_intron+1]\n";
                        }
        } # end until

  if($intronOption == 1)
  {
        	if ( ($snp_start >= $sortchromArray_intron[$cursor_intron]) && ($snp_start < $sortchromArray_intron[$cursor_intron] + $exonbuffer) || ($snp_start <= $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1)  && ($snp_start > $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1 - $exonbuffer) && ($sortchromArray_intron[$cursor_intron] ne "NA") && ($chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} ne "NA")  ) {
        		my $hitbuffer;
        		my $direction;
        		#12 -7 = 5 <6
        		if(($snp_start - $sortchromArray_intron[$cursor_intron]) < $exonbuffer )
        		{
          		  $hitbuffer = $snp_start - $sortchromArray_intron[$cursor_intron] + 1;
          		  $direction = "right";
        		}
        		if(($chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1 - $snp_start) < $exonbuffer )
        		{
          		  $hitbuffer = $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1 - $snp_start + 1;
          		  $direction = "left";
        		}
        		print  OUTFILE1 $line, "\t", "INTRONHIT.$hitbuffer.$direction", "\t", $gene_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"}, "\n";
          		$intronCount++;
                }
  }

  if($intronOption == -1)
  {
        	if ( ($snp_start >= $sortchromArray_intron[$cursor_intron]) && ($snp_start <= $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1)  && ($sortchromArray_intron[$cursor_intron] ne "NA") && ($chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} ne "NA")  ) {
        		print  OUTFILE1 $line, "\t", "INTRONHIT", "\t", $gene_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"}, "\n";
          		$intronCount++;

		}
  }


  open(EXONFILEUTR5, "$exonfile.utr5");
  while(my $lineutr5 = <EXONFILEUTR5>)
  {
    chomp($lineutr5);
    my @lineutr5_array = split(/\t/, $lineutr5);
    if($snp_chr eq $lineutr5_array[0])
    {
      for ($j_utr5=1; $j_utr5<($#lineutr5_array + 1); $j_utr5++)
      {
        push @tempArray_utr5, $lineutr5_array[$j_utr5];
      }
    }
  }
  close(EXONFILEUTR5);

  #utr5start is 0-based, and utr5end is 0-based
  my @sortchromArray_utr5 = @tempArray_utr5;
      $length_utr5= scalar(@sortchromArray_utr5);
      $cursor_utr5=$length_utr5/2;
      $localmin_utr5=0;
      $localmax_utr5=$length_utr5;
      $traverse_utr5=0;
  until (($localmax_utr5-$localmin_utr5)<=1) {
  #print "cursor=|$cursor|\n";
        if ($snp_start <= $sortchromArray_utr5[$cursor_utr5]) {
                        $localmax_utr5=floor($cursor_utr5);
                        &left_utr5($cursor_utr5, $localmin_utr5);}
        if ($snp_start >= $sortchromArray_utr5[$cursor_utr5]) {
                        $localmin_utr5=floor($cursor_utr5);
                        &right_utr5($cursor_utr5, $localmax_utr5);}

        $traverse_utr5++;
                if ($traverse_utr5 >100) {
                        die "Excess Traverse: min=$localmin_utr5\tmax=$localmax_utr5\tcursor=$cursor_utr5\tVPOS=$snp_start\tArrayCursor=$sortchromArray_utr5[$cursor_utr5]\t$sortchromArray_utr5[$cursor_utr5-1]\t$sortchromArray_utr5[$cursor_utr5+1]\n";
                        }
        } # end until

        if ( ($snp_start > $sortchromArray_utr5[$cursor_utr5]) && ($snp_start <= $chr_utr5{"$snp_chr"._."$sortchromArray_utr5[$cursor_utr5]"} + 1)  && ($sortchromArray_utr5[$cursor_utr5] ne "NA") && ($chr_utr5{"$snp_chr"._."$sortchromArray_utr5[$cursor_utr5]"} ne "NA")  ) {
			print  OUTFILE1 $line, "\t", "UTR5HIT", "\t", $gene_utr5{"$snp_chr"._."$sortchromArray_utr5[$cursor_utr5]"}, "\n";
          		$utr5Count++;
                }


  open(EXONFILEUPSTREAM, "$exonfile.upstream");
  while(my $lineupstream = <EXONFILEUPSTREAM>)
  {
    chomp($lineupstream);
    my @lineupstream_array = split(/\t/, $lineupstream);
    if($snp_chr eq $lineupstream_array[0])
    {
      for ($j_upstream=1; $j_upstream<($#lineupstream_array + 1); $j_upstream++)
      {
        push @tempArray_upstream, $lineupstream_array[$j_upstream];
      }
    }
  }
  close(EXONFILEUPSTREAM);

      #upstreamstart is 0-based, and upstreamend is 0-based
    my @sortchromArray_upstream = @tempArray_upstream;
        $length_upstream= scalar(@sortchromArray_upstream);
        $cursor_upstream=$length_upstream/2;
        $localmin_upstream=0;
        $localmax_upstream=$length_upstream;
        $traverse_upstream=0;
until (($localmax_upstream-$localmin_upstream)<=1) {
#print "cursor=|$cursor|\n";
        if ($snp_start <= $sortchromArray_upstream[$cursor_upstream]) {
                        $localmax_upstream=floor($cursor_upstream);
                        &left_upstream($cursor_upstream, $localmin_upstream);}
        if ($snp_start >= $sortchromArray_upstream[$cursor_upstream]) {
                        $localmin_upstream=floor($cursor_upstream);
                        &right_upstream($cursor_upstream, $localmax_upstream);}

        $traverse_upstream++;
                if ($traverse_upstream >100) {
                        die "Excess Traverse: min=$localmin_upstream\tmax=$localmax_upstream\tcursor=$cursor_upstream\tVPOS=$snp_start\tArrayCursor=$sortchromArray_upstream[$cursor_upstream]\t$sortchromArray_upstream[$cursor_upstream-1]\t$sortchromArray_upstream[$cursor_upstream+1]\n";
                        }
        } # end until

        if ( ($snp_start > $sortchromArray_upstream[$cursor_upstream]) && ($snp_start <= $chr_upstream{"$snp_chr"._."$sortchromArray_upstream[$cursor_upstream]"} + 1)  && ($sortchromArray_upstream[$cursor_upstream] ne "NA") && ($chr_upstream{"$snp_chr"._."$sortchromArray_upstream[$cursor_upstream]"} ne "NA")  ) {
			print  OUTFILE1 $line, "\t", "UPSTREAMHIT", "\t", $gene_upstream{"$snp_chr"._."$sortchromArray_upstream[$cursor_upstream]"}, "\n";
          		$upstreamCount++;
                }

  open(EXONFILEUTR3, "$exonfile.utr3");
  while (my $lineutr3 = <EXONFILEUTR3>)
  {
    chomp($lineutr3);
    my @lineutr3_array = split(/\t/, $lineutr3);
    if($snp_chr eq $lineutr3_array[0])
    {
      for ($j_utr3=1; $j_utr3<($#lineutr3_array + 1); $j_utr3++)
      {
        push @tempArray_utr3, $lineutr3_array[$j_utr3];
      }
    }
  }
  close(EXONFILEUTR3);


      #utr3start is 1-based, and utr3end is 1-based
    my @sortchromArray_utr3 = @tempArray_utr3;
        $length_utr3= scalar(@sortchromArray_utr3);
        $cursor_utr3=$length_utr3/2;
        $localmin_utr3=0;
        $localmax_utr3=$length_utr3;
        $traverse_utr3=0;
until (($localmax_utr3-$localmin_utr3)<=1) {
#print "cursor=|$cursor|\n";
        if ($snp_start <= $sortchromArray_utr3[$cursor_utr3]) {
                        $localmax_utr3=floor($cursor_utr3);
                        &left_utr3($cursor_utr3, $localmin_utr3);}
        if ($snp_start >= $sortchromArray_utr3[$cursor_utr3]) {
                        $localmin_utr3=floor($cursor_utr3);
                        &right_utr3($cursor_utr3, $localmax_utr3);}

        $traverse_utr3++;
                if ($traverse_utr3 >100) {
                        die "Excess Traverse: min=$localmin_utr3\tmax=$localmax_utr3\tcursor=$cursor_utr3\tVPOS=$snp_start\tArrayCursor=$sortchromArray_utr3[$cursor_utr3]\t$sortchromArray_utr3[$cursor_utr3-1]\t$sortchromArray_utr3[$cursor_utr3+1]\n";
                        }
        } # end until

        if ( ($snp_start >= $sortchromArray_utr3[$cursor_utr3]) && ($snp_start <= $chr_utr3{"$snp_chr"._."$sortchromArray_utr3[$cursor_utr3]"})  && ($sortchromArray_utr3[$cursor_utr3] ne "NA") && ($chr_utr3{"$snp_chr"._."$sortchromArray_utr3[$cursor_utr3]"} ne "NA")  ) {
       			print  OUTFILE1 $line, "\t", "UTR3HIT", "\t", $gene_utr3{"$snp_chr"._."$sortchromArray_utr3[$cursor_utr3]"}, "\n";
              		$utr3Count++;
               }

  open(EXONFILEDOWNSTREAM, "$exonfile.downstream");
  while (my $linedownstream = <EXONFILEDOWNSTREAM>)
  {
    chomp($linedownstream);
    my @linedownstream_array = split(/\t/, $linedownstream);
    if($snp_chr eq $linedownstream_array[0])
    {
      for ($j_downstream=1; $j_downstream<($#linedownstream_array + 1); $j_downstream++)
      {
        push @tempArray_downstream, $linedownstream_array[$j_downstream];
      }
    }
  }
  close(EXONFILEDOWNSTREAM);


      #downstreamstart is 1-based, and downstreamend is 1-based
    my @sortchromArray_downstream = @tempArray_downstream;
        $length_downstream= scalar(@sortchromArray_downstream);
        $cursor_downstream=$length_downstream/2;
        $localmin_downstream=0;
        $localmax_downstream=$length_downstream;
        $traverse_downstream=0;
until (($localmax_downstream-$localmin_downstream)<=1) {
#print "cursor=|$cursor|\n";
        if ($snp_start <= $sortchromArray_downstream[$cursor_downstream]) {
                        $localmax_downstream=floor($cursor_downstream);
                        &left_downstream($cursor_downstream, $localmin_downstream);}
        if ($snp_start >= $sortchromArray_downstream[$cursor_downstream]) {
                        $localmin_downstream=floor($cursor_downstream);
                        &right_downstream($cursor_downstream, $localmax_downstream);}

        $traverse_downstream++;
                if ($traverse_downstream >100) {
                        die "Excess Traverse: min=$localmin_downstream\tmax=$localmax_downstream\tcursor=$cursor_downstream\tVPOS=$snp_start\tArrayCursor=$sortchromArray_downstream[$cursor_downstream]\t$sortchromArray_downstream[$cursor_downstream-1]\t$sortchromArray_downstream[$cursor_downstream+1]\n";
                        }
        } # end until

        if ( ($snp_start >= $sortchromArray_downstream[$cursor_downstream]) && ($snp_start <= $chr_downstream{"$snp_chr"._."$sortchromArray_downstream[$cursor_downstream]"})  && ($sortchromArray_downstream[$cursor_downstream] ne "NA") && ($chr_downstream{"$snp_chr"._."$sortchromArray_downstream[$cursor_downstream]"} ne "NA")  ) {
			print  OUTFILE1 $line, "\t", "DOWNSTREAMHIT", "\t", $gene_downstream{"$snp_chr"._."$sortchromArray_downstream[$cursor_downstream]"}, "\n";
          		$downstreamCount++;
                }

print "Done for the SNP $snp_chr-$snp_start\n";
  }
}
close SNPFILE;
close(OUTFILE1);

sub left_cds {
    $cursor_cds = shift @_;
    $localmin_cds = shift @_;
    $cursor_cds = floor(($cursor_cds+$localmin_cds)/2);
    return $cursor_cds;
    }

sub right_cds {
    $cursor_cds = shift @_;
    $localmax_cds = shift @_;
    $cursor_cds = floor(($cursor_cds+$localmax_cds)/2);
    return $cursor_cds;
    }

sub left_intron {
    $cursor_intron = shift @_;
    $localmin_intron = shift @_;
    $cursor_intron = floor(($cursor_intron+$localmin_intron)/2);
    return $cursor_intron;
    }

sub right_intron {
    $cursor_intron = shift @_;
    $localmax_intron = shift @_;
    $cursor_intron = floor(($cursor_intron+$localmax_intron)/2);
    return $cursor_intron;
    }

sub left_utr5 {
    $cursor_utr5 = shift @_;
    $localmin_utr5 = shift @_;
    $cursor_utr5 = floor(($cursor_utr5+$localmin_utr5)/2);
    return $cursor_utr5;
    }

sub right_utr5 {
    $cursor_utr5 = shift @_;
    $localmax_utr5 = shift @_;
    $cursor_utr5 = floor(($cursor_utr5+$localmax_utr5)/2);
    return $cursor_utr5;
    }

sub left_utr3 {
    $cursor_utr3 = shift @_;
    $localmin_utr3 = shift @_;
    $cursor_utr3 = floor(($cursor_utr3+$localmin_utr3)/2);
    return $cursor_utr3;
    }

sub right_utr3 {
    $cursor_utr3 = shift @_;
    $localmax_utr3 = shift @_;
    $cursor_utr3 = floor(($cursor_utr3+$localmax_utr3)/2);
    return $cursor_utr3;
    }

sub left_upstream {
    $cursor_upstream = shift @_;
    $localmin_upstream = shift @_;
    $cursor_upstream = floor(($cursor_upstream+$localmin_upstream)/2);
    return $cursor_upstream;
    }

sub right_upstream {
    $cursor_upstream = shift @_;
    $localmax_upstream = shift @_;
    $cursor_upstream = floor(($cursor_upstream+$localmax_upstream)/2);
    return $cursor_upstream;
    }

sub left_downstream {
    $cursor_downstream = shift @_;
    $localmin_downstream = shift @_;
    $cursor_downstream = floor(($cursor_downstream+$localmin_downstream)/2);
    return $cursor_downstream;
    }

sub right_downstream {
    $cursor_downstream = shift @_;
    $localmax_downstream = shift @_;
    $cursor_downstream = floor(($cursor_downstream+$localmax_downstream)/2);
    return $cursor_downstream;
    }
