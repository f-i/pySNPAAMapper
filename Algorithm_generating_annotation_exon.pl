#!/usr/bin/perl

# Algorithm for generating annotation information structure for each exon (perl script to parse exons from ChrAll_knownGene.txt file)
# USAGE: perl Algorithm_generating_annotation_exon.perl ChrAll_knownGene.txt

if(! @ARGV)
{
  die  "Correct Syntax is: Algorithm_generating_annotation_exon.perl <fileneame> \n\nPlease supply a file to parse.  The output file will be the input file name with \".exons\" appended\n\n";
}

$infile = shift @ARGV;

$outfile = $infile . ".exons";

#number of bases on exon boundary
$exonboundary_offset=1;
$up_flank=2000;
$down_flank=500;

open INFILE, "<$infile";
open OUTFILE, ">$outfile";
print OUTFILE "bin\tname\tchrom\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tExonNumber\texStart\texEnd\tintronStart\tintronEnd\tUTR5Start\tUTR5End\tUTR3Start\tUTR3End\tupstreamStart\tupstreamEnd\tdownstreamStart\tdownstreamEnd\n";

while(<INFILE>)  #use the Perl while loop to read a file line by line to the end of the file
{
  chomp;  # avoid \n on last field
  if($_ =~ /^\#/)  #boolean condition; \# if u do want # character, you'll have to quote it with a backslash
  {
    next;  #The next command starts the next iteration of the loop
  }

   ($name, $chrom, $strand, $txStart, $txEnd, $cdsStart,
  $cdsEnd, $exonCount, $exonStarts, $exonEnds, $proteinID, $alignID)=split(/\t/,$_,12);

  @exstarts=split(/,/,$exonStarts);
  @exends=split(/,/,$exonEnds);

  my $bin = "NOINFO";  #for knownGene table parse only!!

  my $CDSStart = "NA";
  my $CDSEnd = "NA";
  my $intronStart = "NA";
  my $intronEnd = "NA";
  my $UTR5Start = "NA";
  my $UTR5End = "NA";
  my $UTR3Start = "NA";
  my $UTR3End = "NA";
  my $upstream = "NA";
  my $downstream = "NA";


  $count=0;
  while($count < $exonCount)
  {
    #if this is the first exon
    if($txStart == $exstarts[$count])
    {
      $intronStart = "NA";
      $intronEnd = "NA";
    }
    else #assign intron to the next Exon
    {
      $intronStart = $exends[$count - 1] + $exonboundary_offset;
      $intronEnd = $exstarts[$count] - $exonboundary_offset;
    }
    #if cdsStart is located in this exon
    if($strand eq "+")
    {
      $ex_num = $count + 1;
      if($exstarts[$count] > $cdsStart) #1
      {
        if($exends[$count] > $cdsEnd) #1.1
        {
          if($exstarts[$count] > $cdsEnd) #1.1.1
          {
            $CDSStart = "NA";
            $CDSEnd = "NA";
            $UTR3Start = $exstarts[$count];
          }
          else#1.1.2
          {
            $CDSStart = $exstarts[$count];
            $CDSEnd = $cdsEnd;
            $UTR3Start = $cdsEnd + 1;
          }
          $UTR3End = $exends[$count];
          if($exends[$count] == $txEnd) #if this is the last exon
          {
            $downstreamStart = $txEnd + 1;
            $downstreamEnd = $txEnd + $down_flank;
          }
          else
          {
            $downstreamStart = "NA";
            $downstreamEnd = "NA";
          }
        }
        else#1.2
        {
          $CDSStart = $exstarts[$count];
          $CDSEnd = $exends[$count];
          $UTR3Start = "NA";
          $UTR3End = "NA";
          $downstreamStart = "NA";
          $downstreamEnd = "NA";
        }
        $UTR5Start = "NA";
        $UTR5End = "NA";
        $upstreamStart = "NA";
        $upstreamEnd = "NA";
      }
      else#2
      {
        if($exends[$count] < $cdsStart) #2.1
        {
          $CDSStart = "NA";
          $CDSEnd = "NA";
          $UTR3Start = "NA";
          $UTR3End = "NA";
          $UTR5End = $exends[$count];
          $downstreamStart = "NA";
          $downstreamEnd = "NA";
        }
        else #2.2
        {
          $CDSStart = $cdsStart;
          $UTR5End = $cdsStart - 1;
          if($exends[$count] < $cdsEnd) #2.2.1
          {
            $CDSEnd = $exends[$count];
            $UTR3Start = "NA";
            $UTR3End = "NA";
            $downstreamStart = "NA";
            $downstreamEnd = "NA";
          }
          else #2.2.2
          {
            $CDSEnd = $cdsEnd;
            $UTR3Start = $cdsEnd + 1;
            $UTR3End = $exends[$count];
            if($exends[$count] == $txEnd) #if this is the last exon
            {
              $downstreamStart = $txEnd + 1;
              $downstreamEnd = $txEnd + $down_flank;
            }
            else
            {
              $downstreamStart = "NA";
              $downstreamEnd = "NA";
            }
          }
        }
        $UTR5Start = $exstarts[$count];
        if($exstarts[$count] == $txStart)#if this is the first exon
        {
          $upstreamStart = $txStart - $up_flank;
          $upstreamEnd = $txStart - 1;
        }
        else
        {
          $upstreamStart = "NA";
          $upstreamEnd = "NA";
        }
      }
    }
    else #negative strand
    {
      $ex_num = $exonCount - $count;
      if($exstarts[$count] > $cdsStart)#1
      {
        if($exends[$count] > $cdsEnd)#1.1
        {
          if($exstarts[$count] > $cdsEnd) #1.1.1
          {
            $CDSStart = "NA";
            $CDSEnd = "NA";
            $UTR5Start = $exstarts[$count];
          }
          else#1.1.2
          {
            $CDSStart = $exstarts[$count];
            $CDSEnd = $cdsEnd;
            $UTR5Start = $cdsEnd + 1;
          }
          $UTR5End = $exends[$count];
          if($exends[$count] == $txEnd) #if this is the last exon
          {
            $upstreamStart = $txEnd + 1;
            $upstreamEnd = $txEnd + $up_flank;
          }
          else
          {
            $upstreamStart = "NA";
            $upstreamEnd = "NA";
          }
        }
        else#1.2
        {
          $CDSStart = $exstarts[$count];
          $CDSEnd = $exends[$count];
          $UTR5Start = "NA";
          $UTR5End = "NA";
          $upstreamStart = "NA";
          $upstreamEnd = "NA";
        }
        $UTR3Start = "NA";
        $UTR3End = "NA";
        $downstreamStart = "NA";
        $downstreamEnd = "NA";
      }
      else#2
      {
#if($name eq "uc009ygf.2")
#{
#  print $exstarts[$count], "\t", $txStart, "\n";
#}
        if($exends[$count] < $cdsStart)#2.1
        {
          $CDSStart = "NA";
          $CDSEnd = "NA";
          $UTR3End = $exends[$count];
          $UTR5Start = "NA";
          $UTR5End = "NA";
          $upstreamStart = "NA";
          $upstreamEnd = "NA";
        }
        else#2.2
        {
          if($exends[$count] < $cdsEnd)#2.2.1
          {
            $CDSEnd = $exends[$count];
            $UTR5Start = "NA";
            $UTR5End = "NA";
            $upstreamStart = "NA";
            $upstreamEnd = "NA";
          }
          else#2.2.2
          {
            $CDSEnd = $cdsEnd;
            $UTR5Start = $cdsEnd + 1;
            $UTR5End = $exends[$count];
            if($exends[$count] == $txEnd) #if this is the last exon
            {
              $upstreamStart = $txEnd + 1;
              $upstreamEnd = $txEnd + $up_flank;
            }
            else
            {
              $upstreamStart = "NA";
              $upstreamEnd = "NA";
            }
          }
          $CDSStart = $cdsStart;
          $UTR3End = $cdsStart - 1;
        }
        $UTR3Start = $exstarts[$count];
        if($exstarts[$count] == $txStart) #if this is the first exon
        {
          $downstreamStart = $txStart - $down_flank;
          $downstreamEnd = $txStart - 1;
        }
        else
        {
          $downstreamStart = "NA";
          $downstreamEnd = "NA";
        }
      }
    }
    print OUTFILE "$bin\t$name\t$chrom\t$txStart\t$txEnd\t$CDSStart\t$CDSEnd\t$ex_num\t$exstarts[$count]\t$exends[$count]\t$intronStart\t$intronEnd\t$UTR5Start\t$UTR5End\t$UTR3Start\t$UTR3End\t$upstreamStart\t$upstreamEnd\t$downstreamStart\t$downstreamEnd\n";
    $count++;
  }##inside while
}###outside while
