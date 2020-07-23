#!/usr/bin/env perl

#Algorithm for preprocessing the gene structure to build annotation structure for each exon
#USAGE: perl Algorithm_preprocessing_exon_annotation_regulatoryRegion.pl ChrAll_knownGene.txt.exons

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

#start the program
@chrArray=`cut -f3 $exonfile | uniq | head -n2`;
foreach $chromosome (@chrArray)
{
  chomp $chromosome;
  if(($chromosome =~ /\_/) || ($chromosome =~ 'chrom'))
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
    print "======================Loop chromosome $chromosome and Hash $exonfile start...\n";
    open EXONFILE, "<$exonfile";
    while(<EXONFILE>)
    {
      my $cdsstart;
      my $cdsend;
      my $intronstart;
      my $intronend;
      my $utr5start;
      my $utr5end;
      my $utr3start;
      my $utr3end;
      my $upstreamstart;
      my $upstreamend;
      my $downstreamstart;
      my $downstreamend;
      chomp;
      if($_ =~ /bin/)
      {
        next;
      }
      else
      {
        @exonline=split(/\t/,$_);
        $exon_gene=$exonline[1];
        $chr=$exonline[2];
        if($chromosome eq $chr)
        {
          $cdsstart=$exonline[5];
          $cdsend=$exonline[6];
          $intronstart=$exonline[10];
          $intronend=$exonline[11];
          $utr5start=$exonline[12];
          $utr5end=$exonline[13];
          $utr3start=$exonline[14];
          $utr3end=$exonline[15];
          $upstreamstart=$exonline[16];
          $upstreamend=$exonline[17];
          $downstreamstart=$exonline[18];
          $downstreamend=$exonline[19];
          if($cdsstart ne "NA" || $cdsend ne "NA")
          {
            push @tempArray_cds, $cdsstart;
            $AllChromStart_cds{$chr} = [ @tempArray_cds ];
            $chr_cds{"$chr"._."$cdsstart"}=$cdsend;
            $gene_cds{"$chr"._."$cdsstart"}=$exon_gene;
            print OUTFILECDSLINK "$chr"._."$cdsstart", "\t", $chr_cds{"$chr"._."$cdsstart"}, "\t";
            print OUTFILECDSLINK $gene_cds{"$chr"._."$cdsstart"}, "\n";
          }
          if($intronstart ne "NA" || $intronend ne "NA")
          {
            push @tempArray_intron, $intronstart;
            $AllChromStart_intron{$chr} = [ @tempArray_intron ];
            $chr_intron{"$chr"._."$intronstart"}=$intronend;
            $gene_intron{"$chr"._."$intronstart"}=$exon_gene;
            print OUTFILEINTRONLINK "$chr"._."$intronstart", "\t", $chr_intron{"$chr"._."$intronstart"}, "\t";
            print OUTFILEINTRONLINK $gene_intron{"$chr"._."$intronstart"}, "\n";
          }
          if($utr5start ne "NA" || $utr5end ne "NA")
          {
            push @tempArray_utr5, $utr5start;
            $AllChromStart_utr5{$chr} = [ @tempArray_utr5 ];
            $chr_utr5{"$chr"._."$utr5start"}=$utr5end;
            $gene_utr5{"$chr"._."$utr5start"}=$exon_gene;
            print OUTFILEUTR5LINK "$chr"._."$utr5start", "\t", $chr_utr5{"$chr"._."$utr5start"}, "\t";
            print OUTFILEUTR5LINK $gene_utr5{"$chr"._."$utr5start"}, "\n";
          }
          if($utr3start ne "NA" || $utr3end ne "NA")
          {
            push @tempArray_utr3, $utr3start;
            $AllChromStart_utr3{$chr} = [ @tempArray_utr3 ];
            $chr_utr3{"$chr"._."$utr3start"}=$utr3end;
            $gene_utr3{"$chr"._."$utr3start"}=$exon_gene;
            print OUTFILEUTR3LINK "$chr"._."$utr3start", "\t", $chr_utr3{"$chr"._."$utr3start"}, "\t";
            print OUTFILEUTR3LINK $gene_utr3{"$chr"._."$utr3start"}, "\n";
          }
          if($upstreamstart ne "NA" || $upstreamend ne "NA")
          {
            push @tempArray_upstream, $upstreamstart;
            $AllChromStart_upstream{$chr} = [ @tempArray_upstream ];
            $chr_upstream{"$chr"._."$upstreamstart"}=$upstreamend;
            $gene_upstream{"$chr"._."$upstreamstart"}=$exon_gene;
            print OUTFILEUPSTREAMLINK "$chr"._."$upstreamstart", "\t", $chr_upstream{"$chr"._."$upstreamstart"}, "\t";
            print OUTFILEUPSTREAMLINK $gene_upstream{"$chr"._."$upstreamstart"}, "\n";
          }
          if($downstreamstart ne "NA" || $downstreamend ne "NA")
          {
            push @tempArray_downstream, $downstreamstart;
            $AllChromStart_downstream{$chr} = [ @tempArray_downstream ];
            $chr_downstream{"$chr"._."$downstreamstart"}=$downstreamend;
            $gene_downstream{"$chr"._."$downstreamstart"}=$exon_gene;
            print OUTFILEDOWNSTREAMLINK "$chr"._."$downstreamstart", "\t", $chr_downstream{"$chr"._."$downstreamstart"}, "\t";
            print OUTFILEDOWNSTREAMLINK $gene_downstream{"$chr"._."$downstreamstart"}, "\n";
          }
        }#firstif
      }
    }
    close EXONFILE;
    print "======================Loop chromosome $chromosome and Hash $exonfile done!!\n";
  }#else
} #for chromsome

foreach $key (sort keys %AllChromStart_cds)
{
  @sortchromArray_cds = sort {$a<=>$b} @{$AllChromStart_cds{$key}};
  print OUTFILECDS $key, "\t";
  for($i=0; $i<($#sortchromArray_cds + 1); $i++)
  {
    print OUTFILECDS $sortchromArray_cds[$i], "\t";
  }
  print OUTFILECDS "\n";
}
foreach $key (sort keys %AllChromStart_intron)
{
  @sortchromArray_intron = sort {$a<=>$b} @{$AllChromStart_intron{$key}};
  print OUTFILEINTRON $key, "\t";
  for($i=0; $i<($#sortchromArray_intron + 1); $i++)
  {
    print OUTFILEINTRON $sortchromArray_intron[$i], "\t";
  }
  print OUTFILEINTRON "\n";
}
foreach $key (sort keys %AllChromStart_utr5)
{
  @sortchromArray_utr5 = sort {$a<=>$b} @{$AllChromStart_utr5{$key}};
  print OUTFILEUTR5 $key, "\t";
  for($i=0; $i<($#sortchromArray_utr5 + 1); $i++)
  {
    print OUTFILEUTR5 $sortchromArray_utr5[$i], "\t";
  }
  print OUTFILEUTR5 "\n";
}
foreach $key (sort keys %AllChromStart_utr3)
{
  @sortchromArray_utr3 = sort {$a<=>$b} @{$AllChromStart_utr3{$key}};
  print OUTFILEUTR3 $key, "\t";
  for($i=0; $i<($#sortchromArray_utr3 + 1); $i++)
  {
    print OUTFILEUTR3 $sortchromArray_utr3[$i], "\t";
  }
  print OUTFILEUTR3 "\n";
}
foreach $key (sort keys %AllChromStart_upstream)
{
  @sortchromArray_upstream = sort {$a<=>$b} @{$AllChromStart_upstream{$key}};
  print OUTFILEUPSTREAM $key, "\t";
  for($i=0; $i<($#sortchromArray_upstream + 1); $i++)
  {
    print OUTFILEUPSTREAM $sortchromArray_upstream[$i], "\t";
  }
  print OUTFILEUPSTREAM "\n";
}
foreach $key (sort keys %AllChromStart_downstream)
{
  @sortchromArray_downstream = sort {$a<=>$b} @{$AllChromStart_downstream{$key}};
  print OUTFILEDOWNSTREAM $key, "\t";
  for($i=0; $i<($#sortchromArray_downstream + 1); $i++)
  {
    print OUTFILEDOWNSTREAM $sortchromArray_downstream[$i], "\t";
  }
  print OUTFILEDOWNSTREAM "\n";
}
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

`cp "$exonfile.cds_link" "$exonfile.cds_link_copy"`;
`cp "$exonfile.intron_link" "$exonfile.intron_link_copy"`;
`cp "$exonfile.utr5_link" "$exonfile.utr5_link_copy"`;
`cp "$exonfile.utr3_link" "$exonfile.utr3_link_copy"`;
`cp "$exonfile.upstream_link" "$exonfile.upstream_link_copy"`;
`cp "$exonfile.downstream_link" "$exonfile.downstream_link_copy"`;

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

my @array_cds_link;
my $i_cds_link=0;
open(EXONFILECDSLINK, "$exonfile.cds_link");
while(<EXONFILECDSLINK>)
{
  if($array_cds_link[$i_cds_link] != 1)
  {
    chomp;
    @exonlinecds=split(/\t/,$_);
    my $chrstart_cds = $exonlinecds[0];
    my $maxStop_cds = $exonlinecds[1];
    my $maxGene_cds = $exonlinecds[2];
    open(EXONFILECDSLINKCOPY, "$exonfile.cds_link_copy");
    my $j_cds_link=0;
    while(<EXONFILECDSLINKCOPY>)
    {
      chomp;
      @exonlinecdslink=split(/\t/,$_);
      if(($chrstart_cds eq $exonlinecdslink[0])&&($exonlinecdslink[1]>$maxStop_cds))
      {
        $maxStop_cds = $exonlinecdslink[1];
        $maxGene_cds = $exonlinecdslink[2];
        $array_cds_link[$j_cds_link] = 1;
      }
      else
      {
        if($chrstart_cds eq $exonlinecdslink[0])
        {
          $array_cds_link[$j_cds_link] = 1;
        }
      }
      $j_cds_link++;
    }
    close(EXONFILECDSLINKCOPY);
    print OUTFILECDSLINKSHRINK $chrstart_cds, "\t", $maxStop_cds, "\n";
    print OUTFILECDSGENE $chrstart_cds, "\t", $maxGene_cds, "\n";
  }
  $i_cds_link++;
}
close EXONFILECDSLINK;

my @array_intron_link;
my $i_intron_link=0;
open(EXONFILEINTRONLINK, "$exonfile.intron_link");
while(<EXONFILEINTRONLINK>)
{
  if($array_intron_link[$i_intron_link] != 1)
  {
    chomp;
    @exonlineintron=split(/\t/,$_);
    my $chrstart_intron = $exonlineintron[0];
    my $maxStop_intron = $exonlineintron[1];
    my $maxGene_intron = $exonlineintron[2];
    open(EXONFILEINTRONLINKCOPY, "$exonfile.intron_link_copy");
    my $j_intron_link=0;
    while(<EXONFILEINTRONLINKCOPY>)
    {
      chomp;
      @exonlineintronlink=split(/\t/,$_);
      if(($chrstart_intron eq $exonlineintronlink[0])&&($exonlineintronlink[1]>$maxStop_intron))
      {
        $maxStop_intron = $exonlineintronlink[1];
        $maxGene_intron = $exonlineintronlink[2];
        $array_intron_link[$j_intron_link] = 1;
      }
      else
      {
        if($chrstart_intron eq $exonlineintronlink[0])
        {
          $array_intron_link[$j_intron_link] = 1;
        }
      }
      $j_intron_link++;
    }
    close(EXONFILEINTRONLINKCOPY);
    print OUTFILEINTRONLINKSHRINK $chrstart_intron, "\t", $maxStop_intron, "\n";
    print OUTFILEINTRONGENE $chrstart_intron, "\t", $maxGene_intron, "\n";
  }
  $i_intron_link++;
}
close EXONFILEINTRONLINK;

my @array_utr5_link;
my $i_utr5_link=0;
open(EXONFILEUTR5LINK, "$exonfile.utr5_link");
while(<EXONFILEUTR5LINK>)
{
  if($array_utr5_link[$i_utr5_link] != 1)
  {
    chomp;
    @exonlineutr5=split(/\t/,$_);
    my $chrstart_utr5 = $exonlineutr5[0];
    my $maxStop_utr5 = $exonlineutr5[1];
    my $maxGene_utr5 = $exonlineutr5[2];
    open(EXONFILEUTR5LINKCOPY, "$exonfile.utr5_link_copy");
    my $j_utr5_link=0;
    while(<EXONFILEUTR5LINKCOPY>)
    {
      chomp;
      @exonlineutr5link=split(/\t/,$_);
      if(($chrstart_utr5 eq $exonlineutr5link[0])&&($exonlineutr5link[1]>$maxStop_utr5))
      {
        $maxStop_utr5 = $exonlineutr5link[1];
        $maxGene_utr5 = $exonlineutr5link[2];
        $array_utr5_link[$j_utr5_link] = 1;
      }
      else
      {
        if($chrstart_utr5 eq $exonlineutr5link[0])
        {
          $array_utr5_link[$j_utr5_link] = 1;
        }
      }
      $j_utr5_link++;
    }
    close(EXONFILEUTR5LINKCOPY);
    print OUTFILEUTR5LINKSHRINK $chrstart_utr5, "\t", $maxStop_utr5, "\n";
    print OUTFILEUTR5GENE $chrstart_utr5, "\t", $maxGene_utr5, "\n";
  }
  $i_utr5_link++;
}
close EXONFILEUTR5LINK;

my @array_utr3_link;
my $i_utr3_link=0;
open(EXONFILEUTR3LINK, "$exonfile.utr3_link");
while(<EXONFILEUTR3LINK>)
{
  if($array_utr3_link[$i_utr3_link] != 1)
  {
    chomp;
    @exonlineutr3=split(/\t/,$_);
    my $chrstart_utr3 = $exonlineutr3[0];
    my $maxStop_utr3 = $exonlineutr3[1];
    my $maxGene_utr3 = $exonlineutr3[2];
    open(EXONFILEUTR3LINKCOPY, "$exonfile.utr3_link_copy");
    my $j_utr3_link=0;
    while(<EXONFILEUTR3LINKCOPY>)
    {
      chomp;
      @exonlineutr3link=split(/\t/,$_);
      if(($chrstart_utr3 eq $exonlineutr3link[0])&&($exonlineutr3link[1]>$maxStop_utr3))
      {
        $maxStop_utr3 = $exonlineutr3link[1];
        $maxGene_utr3 = $exonlineutr3link[2];
        $array_utr3_link[$j_utr3_link] = 1;
      }
      else
      {
        if($chrstart_utr3 eq $exonlineutr3link[0])
        {
          $array_utr3_link[$j_utr3_link] = 1;
        }
      }
      $j_utr3_link++;
    }
    close(EXONFILEUTR3LINKCOPY);
    print OUTFILEUTR3LINKSHRINK $chrstart_utr3, "\t", $maxStop_utr3, "\n";
    print OUTFILEUTR3GENE $chrstart_utr3, "\t", $maxGene_utr3, "\n";
  }
  $i_utr3_link++;
}
close EXONFILEUTR3LINK;

my @array_upstream_link;
my $i_upstream_link=0;
open(EXONFILEUPSTREAMLINK, "$exonfile.upstream_link");
while(<EXONFILEUPSTREAMLINK>)
{
  if($array_upstream_link[$i_upstream_link] != 1)
  {
    chomp;
    @exonlineupstream=split(/\t/,$_);
    my $chrstart_upstream = $exonlineupstream[0];
    my $maxStop_upstream = $exonlineupstream[1];
    my $maxGene_upstream = $exonlineupstream[2];
    open(EXONFILEUPSTREAMLINKCOPY, "$exonfile.upstream_link_copy");
    my $j_upstream_link=0;
    while(<EXONFILEUPSTREAMLINKCOPY>)
    {
      chomp;
      @exonlineupstreamlink=split(/\t/,$_);
      if(($chrstart_upstream eq $exonlineupstreamlink[0])&&($exonlineupstreamlink[1]>$maxStop_upstream))
      {
        $maxStop_upstream = $exonlineupstreamlink[1];
        $maxGene_upstream = $exonlineupstreamlink[2];
        $array_upstream_link[$j_upstream_link] = 1;
      }
      else
      {
        if($chrstart_upstream eq $exonlineupstreamlink[0])
        {
          $array_upstream_link[$j_upstream_link] = 1;
        }
      }
      $j_upstream_link++;
    }
    close(EXONFILEUPSTREAMLINKCOPY);
    print OUTFILEUPSTREAMLINKSHRINK $chrstart_upstream, "\t", $maxStop_upstream, "\n";
    print OUTFILEUPSTREAMGENE $chrstart_upstream, "\t", $maxGene_upstream, "\n";
  }
  $i_upstream_link++;
}
close EXONFILEUPSTREAMLINK;

my @array_downstream_link;
my $i_downstream_link=0;
open(EXONFILEDOWNSTREAMLINK, "$exonfile.downstream_link");
while(<EXONFILEDOWNSTREAMLINK>)
{
  if($array_downstream_link[$i_downstream_link] != 1)
  {
    chomp;
    @exonlinedownstream=split(/\t/,$_);
    my $chrstart_downstream = $exonlinedownstream[0];
    my $maxStop_downstream = $exonlinedownstream[1];
    my $maxGene_downstream = $exonlinedownstream[2];
    open(EXONFILEDOWNSTREAMLINKCOPY, "$exonfile.downstream_link_copy");
    my $j_downstream_link=0;
    while(<EXONFILEDOWNSTREAMLINKCOPY>)
    {
      chomp;
      @exonlinedownstreamlink=split(/\t/,$_);
      if(($chrstart_downstream eq $exonlinedownstreamlink[0])&&($exonlinedownstreamlink[1]>$maxStop_downstream))
      {
        $maxStop_downstream = $exonlinedownstreamlink[1];
        $maxGene_downstream = $exonlinedownstreamlink[2];
        $array_downstream_link[$j_downstream_link] = 1;
      }
      else
      {
        if($chrstart_downstream eq $exonlinedownstreamlink[0])
        {
          $array_downstream_link[$j_downstream_link] = 1;
        }
      }
      $j_downstream_link++;
    }
    close(EXONFILEDOWNSTREAMLINKCOPY);
    print OUTFILEDOWNSTREAMLINKSHRINK $chrstart_downstream, "\t", $maxStop_downstream, "\n";
    print OUTFILEDOWNSTREAMGENE $chrstart_downstream, "\t", $maxGene_downstream, "\n";
  }
  $i_downstream_link++;
}
close EXONFILEDOWNSTREAMLINK;

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
