#!/usr/bin/env perl

# Algorithm for preprocessing the gene structure to build annotation structure
# for each exon

# USAGE:
# perl Algorithm_preprocessing_exon_annotation_regulatoryRegion.pl ChrAll_knownGene.txt.exons

use POSIX;

my ($exonfile) = @ARGV;

%AllChromStart_cds=();
%AllChromStart_intron=();
%AllChromStart_utr5=();
%AllChromStart_utr3=();
%AllChromStart_upstream=();
%AllChromStart_downstream=();

my @sortchromArray_cds;
my @sortchromArray_intron;
my @sortchromArray_utr5;
my @sortchromArray_utr3;
my @sortchromArray_upstream;
my @sortchromArray_downstream;

my $out_cds="$exonfile".".cds";
my $out_intron="$exonfile".".intron";
my $out_utr5="$exonfile".".utr5";
my $out_utr3="$exonfile".".utr3";
my $out_upstream="$exonfile".".upstream";
my $out_downstream="$exonfile".".downstream";

my $out_cds_link="$exonfile".".cds_link";
my $out_intron_link="$exonfile".".intron_link";
my $out_utr5_link="$exonfile".".utr5_link";
my $out_utr3_link="$exonfile".".utr3_link";
my $out_upstream_link="$exonfile".".upstream_link";
my $out_downstream_link="$exonfile".".downstream_link";

open OUTFILECDS, ">$out_cds";
open OUTFILEINTRON, ">$out_intron";
open OUTFILEUTR5, ">$out_utr5";
open OUTFILEUTR3, ">$out_utr3";
open OUTFILEUPSTREAM, ">$out_upstream";
open OUTFILEDOWNSTREAM, ">$out_downstream";

open OUTFILECDSLINK, ">$out_cds_link";
open OUTFILEINTRONLINK, ">$out_intron_link";
open OUTFILEUTR5LINK, ">$out_utr5_link";
open OUTFILEUTR3LINK, ">$out_utr3_link";
open OUTFILEUPSTREAMLINK, ">$out_upstream_link";
open OUTFILEDOWNSTREAMLINK, ">$out_downstream_link";

sub func0 {
  my ($start, $end, $tempArray, $allchromstart, ) = @_;
  if($exonline[$start] ne "NA" || $exonline[$end] ne "NA")
  {
    push @{$tempArray}, $exonline[$start];  #only serves func1
    ${$allchromstart}{$exonline[2]} = [ @{$tempArray} ];  #only serves func1
    print {@_[-1]} "$exonline[2]"._."$exonline[$start]","\t",$exonline[$end],"\t",$exonline[1],"\n";
  }
}

open EXONFILE, '<', $exonfile;
chomp(my @lines = <EXONFILE>);
close EXONFILE;

@chrArray=`cut -f3 $exonfile | uniq`;
foreach $chromosome (@chrArray)
{
  chomp $chromosome;
  if(($chromosome =~ /\_/) || ($chromosome =~ 'chrom'))  #=~ is testing whether underscore is a substring of the value of $chromosome
  {
    next;
  } # this discards the strange chrom pieces that are mapped
  else
  { #initialize the array
    @tempArray_cds=();
    @tempArray_intron=();
    @tempArray_utr5=();
    @tempArray_utr3=();
    @tempArray_upstream=();
    @tempArray_downstream=();
    print "==========Loop chromosome $chromosome and Hash $exonfile start...\n";
    foreach (@lines)
    {
      chomp;
      if($_ =~ /name/)
      {
        next;
      }
      else
      {
        @exonline=split(/\t/,$_);
        if($chromosome eq $exonline[2])
        {
          func0(5,6,tempArray_cds,AllChromStart_cds,OUTFILECDSLINK);
          func0(10,11,tempArray_intron,AllChromStart_intron,OUTFILEINTRONLINK);
          func0(12,13,tempArray_utr5,AllChromStart_utr5,OUTFILEUTR5LINK);
          func0(14,15,tempArray_utr3,AllChromStart_utr3,OUTFILEUTR3LINK);
          func0(16,17,tempArray_upstream,AllChromStart_upstream,OUTFILEUPSTREAMLINK);
          func0(18,19,tempArray_downstream,AllChromStart_downstream,OUTFILEDOWNSTREAMLINK);
        }#firstif
      }
    }
    print "==========Loop chromosome $chromosome and Hash $exonfile done!!\n";
  }#else
} #for chromsome

sub func1 {
  my ($acstart, $f) = @_;
  foreach $key (sort keys %{$acstart})
  {
    @sortchromArray = sort {$a<=>$b} @{${$acstart}{$key}};
    print $f $key, "\t";
    for($i=0; $i<($#sortchromArray + 1); $i++)  # $#array is the last index of @array
    {
      print $f $sortchromArray[$i], "\t";
    }
    print $f "\n";
  }
}
func1(AllChromStart_cds,OUTFILECDS);
func1(AllChromStart_intron,OUTFILEINTRON);
func1(AllChromStart_utr5,OUTFILEUTR5);
func1(AllChromStart_utr3,OUTFILEUTR3);
func1(AllChromStart_upstream,OUTFILEUPSTREAM);
func1(AllChromStart_downstream,OUTFILEDOWNSTREAM);

close(OUTFILECDS);
close(OUTFILEINTRON);
close(OUTFILEUTR5);
close(OUTFILEUTR3);
close(OUTFILEUPSTREAM);
close(OUTFILEDOWNSTREAM);

close(OUTFILECDSLINK);
close(OUTFILEINTRONLINK);
close(OUTFILEUTR5LINK);
close(OUTFILEUTR3LINK);
close(OUTFILEUPSTREAMLINK);
close(OUTFILEDOWNSTREAMLINK);

my $out_cds_link_shrink="$exonfile".".cds_link_shrink";
my $out_intron_link_shrink="$exonfile".".intron_link_shrink";
my $out_utr5_link_shrink="$exonfile".".utr5_link_shrink";
my $out_utr3_link_shrink="$exonfile".".utr3_link_shrink";
my $out_upstream_link_shrink="$exonfile".".upstream_link_shrink";
my $out_downstream_link_shrink="$exonfile".".downstream_link_shrink";

my $out_cds_gene="$exonfile".".cds_gene";
my $out_intron_gene="$exonfile".".intron_gene";
my $out_utr5_gene="$exonfile".".utr5_gene";
my $out_utr3_gene="$exonfile".".utr3_gene";
my $out_upstream_gene="$exonfile".".upstream_gene";
my $out_downstream_gene="$exonfile".".downstream_gene";

open OUTFILECDSLINKSHRINK, ">$out_cds_link_shrink";
open OUTFILEINTRONLINKSHRINK, ">$out_intron_link_shrink";
open OUTFILEUTR5LINKSHRINK, ">$out_utr5_link_shrink";
open OUTFILEUTR3LINKSHRINK, ">$out_utr3_link_shrink";
open OUTFILEUPSTREAMLINKSHRINK, ">$out_upstream_link_shrink";
open OUTFILEDOWNSTREAMLINKSHRINK, ">$out_downstream_link_shrink";

open OUTFILECDSGENE, ">$out_cds_gene";
open OUTFILEINTRONGENE, ">$out_intron_gene";
open OUTFILEUTR5GENE, ">$out_utr5_gene";
open OUTFILEUTR3GENE, ">$out_utr3_gene";
open OUTFILEUPSTREAMGENE, ">$out_upstream_gene";
open OUTFILEDOWNSTREAMGENE, ">$out_downstream_gene";

sub func2 {
  my ($f1, $f2, $f3) = @_;
  open EXONFILELINK, '<', $f1;
  chomp(my @lins = <EXONFILELINK>);
  close EXONFILELINK;
  my @array_link;
  my $i=0;
  foreach (@lins)
  {
    if($array_link[$i] != 1)
    {
      chomp;
      @exonline=split(/\t/,$_);
      my $chrstart = $exonline[0];
      my $maxStop = $exonline[1];
      my $maxGene = $exonline[2];
      my $j=0;
      foreach (@lins)
      {
        chomp;
        @exonlinelink=split(/\t/,$_);
        if(($chrstart eq $exonlinelink[0])&&($exonlinelink[1]>$maxStop))
        {
          $maxStop = $exonlinelink[1];
          $maxGene = $exonlinelink[2];
          $array_link[$j] = 1;
        }
        else
        {
          if($chrstart eq $exonlinelink[0])
          {
            $array_link[$j] = 1;
          }
        }
        $j++;
      }
      print $f2 $chrstart, "\t", $maxStop, "\n";
      print $f3 $chrstart, "\t", $maxGene, "\n";
    }
    $i++;
  }
}

func2("$exonfile.cds_link", OUTFILECDSLINKSHRINK, OUTFILECDSGENE);
func2("$exonfile.intron_link", OUTFILEINTRONLINKSHRINK, OUTFILEINTRONGENE);
func2("$exonfile.utr5_link", OUTFILEUTR5LINKSHRINK, OUTFILEUTR5GENE);
func2("$exonfile.utr3_link", OUTFILEUTR3LINKSHRINK, OUTFILEUTR3GENE);
func2("$exonfile.upstream_link", OUTFILEUPSTREAMLINKSHRINK, OUTFILEUPSTREAMGENE);
func2("$exonfile.downstream_link", OUTFILEDOWNSTREAMLINKSHRINK, OUTFILEDOWNSTREAMGENE);

close(OUTFILECDSLINKSHRINK);
close(OUTFILEINTRONLINKSHRINK);
close(OUTFILEUTR5LINKSHRINK);
close(OUTFILEUTR3LINKSHRINK);
close(OUTFILEUPSTREAMLINKSHRINK);
close(OUTFILEDOWNSTREAMLINKSHRINK);

close(OUTFILECDSGENE);
close(OUTFILEINTRONGENE);
close(OUTFILEUTR5GENE);
close(OUTFILEUTR3GENE);
close(OUTFILEUPSTREAMGENE);
close(OUTFILEDOWNSTREAMGENE);
