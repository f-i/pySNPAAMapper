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

if(($#ARGV + 1)== 2)
{
  ($exonfile,$snpfile) = @ARGV;
  print "The program assumes that you do NOT want to report how far the variant falls in the exon boundary.\n";
  $intronOption = -1;
}
elsif(($#ARGV + 1) == 3)
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
sub func0 {
  my ($f, $f0, $f1) = @_;
  open($f, "$f0");
  while(<$f>)
  {
    chomp;
    @exonline=split(/\t/,$_);
    ${$f1}{$exonline[0]} = $exonline[1];
  }
  close $f;
}

func0(EXONFILECDSLINK,$exonfile.cds_link_shrink,chr_cds);
func0(EXONFILECDSGENE,$exonfile.cds_gene,gene_cds);
func0(EXONFILEINTRONLINK,$exonfile.intron_link_shrink,chr_intron);
func0(EXONFILEINTRONGENE,$exonfile.intron_gene,gene_intron);
func0(EXONFILEUTR5LINK,$exonfile.utr5_link_shrink,chr_utr5);
func0(EXONFILEUTR5GENE,$exonfile.utr5_gene,gene_utr5);
func0(EXONFILEUTR3LINK,$exonfile.utr3_link_shrink,chr_utr3);
func0(EXONFILEUTR3GENE,$exonfile.utr3_gene,gene_utr3);
func0(EXONFILEUPSTREAMLINK,$exonfile.upstream_link_shrink,chr_upstream);
func0(EXONFILEUPSTREAMGENE,$exonfile.upstream_gene,gene_upstream);
func0(EXONFILEDOWNSTREAMLINK,$exonfile.downstream_link_shrink,chr_downstream);
func0(EXONFILEDOWNSTREAMGENE,$exonfile.downstream_gene,gene_downstream);

#map vcf variant back onto genome location and report the hit class
open OUTFILE1, ">$snpfile.append";

sub left {
  my ($f,$f0) = @_;
  ${$f} = shift @_;
  ${$f0} = shift @_;
  ${$f} = floor((${$f}+${$f0})/2);
  return ${$f};
  }

sub right {
  my ($f,$f0) = @_;
  ${$f} = shift @_;
  ${$f0} = shift @_;
  ${$f} = floor((${$f}+${$f0})/2);
  return ${$f};
  }

sub func1 {
  my ($f,$f0,$f1,$f2,$f3,$f4) = @_;
  open(EXONFILE, "$exonfile.$f");
  while(my $line = <EXONFILE>)
  {
    chomp($line);
    my @line_array = split(/\t/, $line);
    if($snp_chr eq $line_array[0])  #given the same genome chromosome as variant's chrom
    { #only store all of start postions in the array (without chromosome information included) for binary search
      for($j=1; $j<($#line_array + 1); $j++)
      {
        push @{$f1}, $line_array[$j];
      }
    }
  }
  close(EXONFILE);
  #my @{$f0};
  @{$f0} = @{$f1};
  $length=scalar(@{$f0});
  ${$f2}=$length/2;
  $localmin=0;
  $localmax=$length;
  $traverse=0;
  until(($localmax-$localmin)<=1)
  {
    #print "cursor=|$cursor|\n";
    if($snp_start <= ${$f0}[${$f2}])
    {
      $localmax=floor(${$f2});
      &left($f2, localmin);
    }
    if($snp_start >= ${$f0}[${$f2}])
    {
      $localmin=floor(${$f2});
      &right($f2, localmax);
    }
    $traverse++;
    if($traverse>100)
    {
      die "Excess Traverse: min=$localmin\tmax=$localmax\tcursor=${$f2}\tVPOS=$snp_start\tArrayCursor=${$f0}[${$f2}]\t${$f0}[${$f2}-1]\t${$f0}[${$f2}+1]\n";
    }
  } # end until
}

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

    # cdsstart is 0-based, and cdsend is 1-based
    func1(cds,sortchromArray_cds,tempArray_cds,cursor_cds);
    if(($snp_start > $sortchromArray_cds[$cursor_cds]) && ($snp_start<=$chr_cds{"$snp_chr"._."$sortchromArray_cds[$cursor_cds]"}) && ($sortchromArray_cds[$cursor_cds] ne "NA") && ($chr_cds{"$snp_chr"._."$sortchromArray_cds[$cursor_cds]"} ne "NA"))
    {
      print OUTFILE1 $line, "\t", "CDSHIT", "\t", $gene_cds{"$snp_chr"._."$sortchromArray_cds[$cursor_cds]"}, "\n";
      $cdsCount++;
    }

    # inronstart is 1-based, add intronhit buffer number for parse_discrepancy_jan_inronSeparate_final.perl because intronend is 0-based
    func1(intron,sortchromArray_intron,tempArray_intron,cursor_intron);
    if($intronOption == 1)
    {
      if(($snp_start >= $sortchromArray_intron[$cursor_intron]) && ($snp_start<$sortchromArray_intron[$cursor_intron] + $exonbuffer) || ($snp_start <= $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1) && ($snp_start > $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1 - $exonbuffer) && ($sortchromArray_intron[$cursor_intron] ne "NA") && ($chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} ne "NA"))
      {
        my $hitbuffer;
        my $direction;
        #12 -7 = 5 <6
        if(($snp_start - $sortchromArray_intron[$cursor_intron]) < $exonbuffer)
        {
          $hitbuffer = $snp_start - $sortchromArray_intron[$cursor_intron] + 1;
          $direction = "right";
        }
        if(($chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1 - $snp_start) < $exonbuffer)
        {
          $hitbuffer = $chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1 - $snp_start + 1;
          $direction = "left";
        }
        print OUTFILE1 $line, "\t", "INTRONHIT.$hitbuffer.$direction", "\t", $gene_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"}, "\n";
        $intronCount++;
      }
    }
    if($intronOption == -1)
    {
      if(($snp_start >= $sortchromArray_intron[$cursor_intron]) && ($snp_start<=$chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} + 1) && ($sortchromArray_intron[$cursor_intron] ne "NA") && ($chr_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"} ne "NA"))
      {
        print OUTFILE1 $line, "\t", "INTRONHIT", "\t", $gene_intron{"$snp_chr"._."$sortchromArray_intron[$cursor_intron]"}, "\n";
        $intronCount++;
      }
    }

    # utr5start is 0-based, and utr5end is 0-based
    func1(utr5,sortchromArray_utr5,tempArray_utr5,cursor_utr5);
    if(($snp_start > $sortchromArray_utr5[$cursor_utr5]) && ($snp_start<=$chr_utr5{"$snp_chr"._."$sortchromArray_utr5[$cursor_utr5]"} + 1) && ($sortchromArray_utr5[$cursor_utr5] ne "NA") && ($chr_utr5{"$snp_chr"._."$sortchromArray_utr5[$cursor_utr5]"} ne "NA"))
    {
      print OUTFILE1 $line, "\t", "UTR5HIT", "\t", $gene_utr5{"$snp_chr"._."$sortchromArray_utr5[$cursor_utr5]"}, "\n";
      $utr5Count++;
    }

    # upstreamstart is 0-based, and upstreamend is 0-based
    func1(upstream,sortchromArray_upstream,tempArray_upstream,cursor_upstream);
    if(($snp_start > $sortchromArray_upstream[$cursor_upstream]) && ($snp_start<=$chr_upstream{"$snp_chr"._."$sortchromArray_upstream[$cursor_upstream]"} + 1) && ($sortchromArray_upstream[$cursor_upstream] ne "NA") && ($chr_upstream{"$snp_chr"._."$sortchromArray_upstream[$cursor_upstream]"} ne "NA"))
    {
      print OUTFILE1 $line, "\t", "UPSTREAMHIT", "\t", $gene_upstream{"$snp_chr"._."$sortchromArray_upstream[$cursor_upstream]"}, "\n";
      $upstreamCount++;
    }

    # utr3start is 1-based, and utr3end is 1-based
    func1(utr3,sortchromArray_utr3,tempArray_utr3,cursor_utr3);
    if(($snp_start >= $sortchromArray_utr3[$cursor_utr3]) && ($snp_start<=$chr_utr3{"$snp_chr"._."$sortchromArray_utr3[$cursor_utr3]"}) && ($sortchromArray_utr3[$cursor_utr3] ne "NA") && ($chr_utr3{"$snp_chr"._."$sortchromArray_utr3[$cursor_utr3]"} ne "NA"))
    {
      print OUTFILE1 $line, "\t", "UTR3HIT", "\t", $gene_utr3{"$snp_chr"._."$sortchromArray_utr3[$cursor_utr3]"}, "\n";
      $utr3Count++;
    }

    # downstreamstart is 1-based, and downstreamend is 1-based
    func1(downstream,sortchromArray_downstream,tempArray_downstream,cursor_downstream);
    if(($snp_start >= $sortchromArray_downstream[$cursor_downstream]) && ($snp_start<=$chr_downstream{"$snp_chr"._."$sortchromArray_downstream[$cursor_downstream]"}) && ($sortchromArray_downstream[$cursor_downstream] ne "NA") && ($chr_downstream{"$snp_chr"._."$sortchromArray_downstream[$cursor_downstream]"} ne "NA"))
    {
      print OUTFILE1 $line, "\t", "DOWNSTREAMHIT", "\t", $gene_downstream{"$snp_chr"._."$sortchromArray_downstream[$cursor_downstream]"}, "\n";
      $downstreamCount++;
    }
    print "Done for the SNP $snp_chr-$snp_start\n";
  }
}
close SNPFILE;
close(OUTFILE1);
