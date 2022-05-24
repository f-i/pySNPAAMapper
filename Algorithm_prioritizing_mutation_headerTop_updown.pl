#!/usr/bin/perl

#Algorithm for prioritizing mutation effects
#USAGE: perl Algorithm_prioritizing_mutation.pl VCF_input_file_in_tab_delimited_format.txt.append.out.txt

my ($inputfile1) = @ARGV;

my $comma = ",";
my $i = 1;
my $rowHash = ();
my $scoreHash = ();
my $infoHash = ();

open OUTFILE01, ">$inputfile1.prioritized_out";


open INFILE01, "$inputfile1";
while (my $line01 = <INFILE01>) {
  chomp($line01);
  my @line01_array = split('\t', $line01);
  if ($line01_array[0] eq "Sample")
  {
    print OUTFILE01 $line01, "\n";
  }
  else
  {
  $rowHash{$i} = $line01;
  #get rid of right reciprocal pairs
  if($line01_array[7] eq "SNP")
  {
    if($line01_array[12] eq "CDSHIT")
    {
      if($line01_array[9] =~ m/$comma/)
      {
        my @snp_chars = split(/$comma/, $line01_array[9]);
        #$case_counter = $#snp_chars + 1;
        for(my $k=0; $k<($#snp_chars + 1); $k++)
        {
          if($snp_chars[$k] eq "NSM")
          {
            $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 6;
          }
          if($snp_chars[$k] eq "NSN")
          {
            $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 5;
          }
          if($snp_chars[$k] eq "SYN")
          {
            $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 4;
          }
        }
        $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
      }
      else
      {
        if($line01_array[9] eq "NSM")
        {
          $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 6;
        }
        if($line01_array[9] eq "NSN")
        {
          $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 5;
        }
        if($line01_array[9] eq "SYN")
        {
          $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 4;
        }
        $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
      }
    }
    elsif($line01_array[12] eq "UPSTREAMHIT")
    {
      #use combination key: chr_SNPLocation_geneSymbol_rowNumber
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 3;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "DOWNSTREAMHIT")
    {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 2;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "UTR5HIT")
    {
      #use combination key: chr_SNPLocation_geneSymbol_rowNumber
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 1;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "UTR3HIT")
    {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 0;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    else
     {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = -1;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
  }#if
  else #INDEL case
  {
    if($line01_array[12] eq "CDSHIT")
    {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 3.5;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "UPSTREAMHIT")
    {
      #use combination key: chr_SNPLocation_geneSymbol_rowNumber
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 3;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "DOWNSTREAMHIT")
    {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 2;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "UTR5HIT")
    {
      #use combination key: chr_SNPLocation_geneSymbol_rowNumber
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 1;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    elsif($line01_array[12] eq "UTR3HIT")
    {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = 0;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
    else
    {
      $scoreHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = -1;
      $infoHash{$line01_array[1]."+".$line01_array[2]."+".$line01_array[3]."+".$i} = [@line01_array];
    }
  }
  $i++;
  }
}
close(INFILE01);

my %outputHash = ();

for my $key21 (keys(%infoHash))
{
  my @key21_array = split(/\+/, $key21);
  #$outputHash{$key21_array[0]."+".$key21_array[1]."+".$key21_array[2]} = 0;;
  $outputHash{$key21_array[0]."+".$key21_array[1]."+".$key21_array[2]} = $key21_array[3];
  for my $key22 (keys(%infoHash))
  {
    my @key22_array = split(/\+/, $key22);
    if (($key21_array[0] eq $key22_array[0]) && ($key21_array[1] eq $key22_array[1]) && ($key21_array[2] eq $key22_array[2]))
    {
      if($scoreHash{$key22_array[0]."+".$key22_array[1]."+".$key22_array[2]."+".$key22_array[3]} > $scoreHash{$key21_array[0]."+".$key21_array[1]."+".$key21_array[2]."+".$key21_array[3]})
      {
        $outputHash{$key21_array[0]."+".$key21_array[1]."+".$key21_array[2]} = $key22_array[3]; #save the row number with the highest score
      }
    } #if
  }
}

for my $keyOut (keys(%outputHash))
{
  my @keyOut_array = split(/\+/, $keyOut);
  print OUTFILE01 $rowHash{$outputHash{$keyOut}}, "\n";
}

close(OUTFILE01);



