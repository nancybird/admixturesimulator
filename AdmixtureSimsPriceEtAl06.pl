## SIMULATES ADMIXED INDS ACCORDING TO PRICE ET AL 2006 (HAPMIX) PAPER PROTOCOL

## usage:  perl AdmixtureSimsPriceEtAl06.pl

###########################################
## INPUT:

$simnamePRE="example_sim"; #name of sim prefix
@lambdavec=(15); #number of generations ago
@probcopyvec=(0.5); ##percentage mixture
$nreps=30;          ## number of simulated inds for each gen/prop combo
@sourcevec=("French", "Malawi"); #populations you are mixing
@nhapssourcevec=(40,68);       ## use first XX haps from each source's ordered inds to simulate, rest as donors in painting (!EACH ENTRY NEEDS TO BE AN EVEN NUMBER!)
$idfile=""; ##Chromopainter id file
$genofilePRE=""; ##Chromopainter inp file before chromosome
$genofilePOST=""; ##Chromopainter inp file after chromosome

                               ## THIS PROGRAM'S OUTPUT FILES:
$outfileMID=".Chr";
$outfilePOST=".chromopainter.haps";
$outfilePOSTid=".idfile.txt";
$outfilePOSTtruth=".TRUTH.samples.out";
$pathwayTRUTH = "";
$pathwaySIM = "";

                               ## YOUR CHROMOPAINTER INPUT:
$recomratesPRE=""; ##Chromopainter recomrates file
$recomratesPOST=""; ##Chromopainter recomrates file
@chromovec=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22");
$ploidy=2;

$nooverlap='no';    ## TRY TO MAKE SURE THAT FOR ANY GIVEN REGION EACH DONOR IND IS COPIED AT MOST ONCE
$analyseIND='no';   ## ANALYSE EACH IND SEPARATELY?

############################################
## PROGRAM:

$numchrom=@chromovec;
$numsources=@sourcevec;
$numlambda=@lambdavec;
$numcopyprob=@probcopyvec;
$nhaps=$numcopyprob*$numlambda*$nreps*2;

        ## (I) GET POP LABELS:
open(IN,"$idfile");
@uniquepops=();
@nindpops=();
@idvec=();
@tokeepvec=();
@popvec=();
@popvecNUM=();
$numpops=0;
while(<IN>)
{
    $line=$_;
    @linearray=split(/\s+/,$line);
    push(@idvec,$linearray[0]);
    push(@popvec,$linearray[1]);
    if ($linearray[2]==0)
    {
        push(@tokeepvec,0);
        push(@popvecNUM,-9);
    }
    if ($linearray[2]==1)
    {
        push(@tokeepvec,1);
        $popfind=0;
        for ($k=0; $k < $numpops; $k+=1)
        {
            if ($uniquepops[$k] eq $linearray[1])
            {
                $popfind=1;
                $nindpops[$k]=$nindpops[$k]+1;
                push(@popvecNUM,$k);
                last;
            }
        }
        if ($popfind==0)
        {
            push(@uniquepops,$linearray[1]);
            push(@nindpops,1);
            push(@popvecNUM,$numpops);
            $numpops=$numpops+1;
        }
    }
}
close(IN);
$nindTOT=@idvec;
print "$nindTOT $numpops\n";
for ($k=0; $k < $numpops; $k+=1)
{
    print "$k $uniquepops[$k] $nindpops[$k]\n";
    for ($i=0; $i < $numsources; $i+=1)
    {
        if ($uniquepops[$k] eq $sourcevec[$i])
        {
            if (($ploidy*$nindpops[$k])<$nhapssourcevec[$i])
            {
                print "Want to simulate using $nhapssourcevec[$i] haps from $sourcevec[$i], but only have ($ploidy*$nindpops[$k]) total haps! Exiting...\n";
                die;
            }
            last;
        }
    }
}

        ## (II) GET INDS TO USE FOR SIMULATING:
@tokeepvecNEW=();
@numsourceindsSELECT=();
@allsourceSELECT=();
for ($a=0; $a < $numsources; $a+=1)
{
    push(@numsourceindsSELECT,0);
}
for ($a=0; $a < $numsources; $a+=1)
{
    for ($i=0; $i < $nindTOT; $i+=1)
    {
        if ($tokeepvec[$i]==1)
        {
            if ($popvec[$i] eq $sourcevec[$a])
            {
                if ($numsourceindsSELECT[$a]<$nhapssourcevec[$a])
                {
                    for ($m=0; $m < $ploidy; $m+=1)
                    {
                        push(@allsourceSELECT,$ploidy*$i+$m);
                    }
                    $numsourceindsSELECT[$a]=$numsourceindsSELECT[$a]+$ploidy;
                    if ($numsourceindsSELECT[$a]>=$nhapssourcevec[$a])
                    {
                        last;
                    }
                }
            }
        }
    }
}
$numselect=@allsourceSELECT/$ploidy;
print "$numselect\n";
for ($a=0; $a < $numsources; $a+=1)
{
    print "$sourcevec[$a] $numsourceindsSELECT[$a]\n";
}
@tokeepvecNEW=();
for ($i=0; $i < $nindTOT; $i+=1)
{
    $selectfind=0;
    for ($j=0; $j < $numselect; $j+=1)
    {
        if (($i*$ploidy)==$allsourceSELECT[($j*$ploidy)])
        {
            $selectfind=1;
            last;
        }
    }
    if ($selectfind==0 && $tokeepvec[$i]==1)
    {
        push(@tokeepvecNEW,1);
    }
    if ($selectfind==1 || $tokeepvec[$i]==0)
    {
        push(@tokeepvecNEW,0);
    }
}
print "@allsourceSELECT\n";
for ($i=0; $i < $numselect; $i+=1)
{
    print "$i $numselect $allsourceSELECT[($i*$ploidy)] $popvec[($allsourceSELECT[($i*$ploidy)]/$ploidy)]\n";
}

        ## (III) PRINT OUT NEW ID FILE FOR SIMS:
open(OUT,">${pathwaySIM}${simnamePRE}${outfilePOSTid}");
for ($i=0; $i < $nindTOT; $i+=1)
{
    print OUT "$idvec[$i] $popvec[$i] $tokeepvecNEW[$i]\n";
}
for ($g=0; $g < $numlambda; $g+=1)
{
    for ($p=0; $p < $numcopyprob; $p+=1)
    {
        for ($i=0; $i < $nreps; $i+=1)
        {
            $j=$i+1;
            if ($analyseIND eq 'yes')
            {
                print OUT "${simnamePRE}_$lambdavec[$g]gen_$probcopyvec[$p]percent_rep${j} ${simnamePRE}_$lambdavec[$g]gen_$probcopyvec[$p]percent_rep${j} 1\n";
            }
            if ($analyseIND ne 'yes')
            {
                print OUT "${simnamePRE}_$lambdavec[$g]gen_$probcopyvec[$p]percent_rep${j} ${simnamePRE}_$lambdavec[$g]gen_$probcopyvec[$p]percent 1\n";
            }
        }
    }
}
close(OUT);
$totalhaps=($ploidy*$nindTOT)+$nhaps;

        ## (IV) SIMULATE AND PRINT OUT HAPLOTYPES:
for ($b=0; $b < $numchrom; $b+=1)
{
    $chromnum = $chromovec[$b];
    print "chromo $chromnum\n";

    $outfile = "${pathwaySIM}${simnamePRE}${outfileMID}${chromnum}${outfilePOST}";
    open(OUT,">$outfile");
    $outfile2 = "${pathwayTRUTH}${simnamePRE}${outfileMID}${chromnum}${outfilePOSTtruth}";
    open(OUT2,">$outfile2");

         ## (i) GET POSITIONS & RECOMB MAP:
    $infile = "${recomratesPRE}${chromnum}${recomratesPOST}";
    open(IN,"$infile");
    $line=<IN>;   ## header
    @posvec=();
    @recomrateperbp=();
    while(<IN>)
    {
        $line=$_;
        @linearray=split(/\s+/,$line);
        push(@posvec,$linearray[0]);
        push(@recomrateperbp,$linearray[1]);
    }
    close(IN);

          ## PRINT OUT ALL HAPLOTYPES IN ORIGINAL FILE:
    print OUT "$totalhaps\n";
    open(IN,"${genofilePRE}${chromnum}${genofilePOST}");
    #print "${genofilePRE}${chromnum}${genofilePOST}\n";
    $line=<IN>;
    $line=<IN>;
    @linearray=split(/\s+/,$line);
    $nsites=$linearray[0];
    $line=<IN>;
    @posvec=split(/\s+/,$line);
    print OUT "$nsites\n";
    print OUT "@posvec\n";
    for ($i=0; $i < $nindTOT; $i+=1)
    {
        for ($m=0; $m < $ploidy; $m+=1)
        {
            $line=<IN>;
            print OUT "$line";
        }
    }
    close(IN);

          ## (ii) READ IN DONOR POP DATA:
    $sourceHAPcount=0;
    $totalsourceHAPcount=0;
    @genomatDONOR=();
    for ($a=0; $a < $numsources; $a+=1)
    {
        $totalsourceHAPcount=$totalsourceHAPcount+$nhapssourcevec[$a];
        open(IN,"${genofilePRE}${chromnum}${genofilePOST}");
        #print "$a $numsources ${genofilePRE}${chromnum}${genofilePOST}\n";
        $line=<IN>;
        $line=<IN>;
        @linearray=split(/\s+/,$line);
        $nsites=$linearray[0];
        $line=<IN>;
        @posvec=split(/\s+/,$line);
        shift(@posvec);
        for ($i=0; $i < $nindTOT; $i+=1)
        {
            for ($m=0; $m < $ploidy; $m+=1)
            {
                $line=<IN>;
                if (($ploidy*$i+$m) == $allsourceSELECT[$sourceHAPcount])
                {
                    #print "$sourcevec[$a] $i $m $popvec[$i] $sourceHAPcount $allsourceSELECT[$sourceHAPcount]\n";
                    @linearray=split(//,$line);
                    for ($j=0; $j < $nsites; $j+=1)
                    {
                        $genomatDONOR[$sourceHAPcount][$j]=$linearray[$j];
                    }
                    $sourceHAPcount=$sourceHAPcount+1;
                }
            }
            if ($sourceHAPcount>=$totalsourceHAPcount)
            {
                last;
            }
        }
        close(IN);
    }
    print "total number of haplotypes to simulate = $sourceHAPcount\n";
#    for ($i=0; $i < $sourceHAPcount; $i+=1)
#    {
#       for ($j=0; $j < $nsites; $j+=1)
#       {
#           print "$genomatDONOR[$i][$j]";
#       }
#       print "\n";
#   }
#    die;

    print OUT2 "EM_iter = NA, nsamples = 1, N_e fixed, populations = (@sourcevec), haplotypes = (@nhapssourcevec), time of admixture = (@lambdavec), copy probs = (@probcopyvec), numchroms = $nhaps, TRUE RESULTS, ALL HAPS ADMIXED (RANDOMLY PICKED NEW HAP FOR EACH CHUNK TRYING TO NOT PICK OTHERS PREVIOUSLY PICKED), dataset = $outfile\n";
    print OUT2 "CHROMO ${chromnum}\n";

           ## (iii) RANDOMLY GENERATE HAPLOTYPES AND PRINT OUT:
    @indDONOR=();
    $outofdonorcount=0;
    $hapcount=0;
    for ($h=0; $h < $numlambda; $h+=1)
    {
        for ($p=0; $p < $numcopyprob; $p+=1)
        {
            @probcopy=();
            push(@probcopy,$probcopyvec[$p],1.0-$probcopyvec[$p]);
            for ($i=0; $i < $nreps; $i+=1)
            {
                for ($m=0; $m < $ploidy; $m+=1)
                {
                    $hapcount=$hapcount+1;
                    print "simulating hap $hapcount of $nhaps ($i $m) .....\n";

                    $geneticlengthHAP = 0;
                    $totalgeneticdist = 0;
                    $start = 0;
                    $j = $start;
                    $hapind = $m+1;
                    print OUT2 "HAP ${simnamePRE}_$lambdavec[$h]gen_$probcopyvec[$p]percent_rep${j}_hap${hapind}\n1";
                    while ($j < $nsites)
                    {
                                 ## (I) PICK POP TO COPY:
                        $random_number_pop = rand(1);
                        $sumprobs = 0.0;
                        $ndonorsadd=0;
                        for ($k=0; $k < $numsources; $k+=1)
                        {
                            $sumprobs = $sumprobs + $probcopy[$k];
                            if ($sumprobs >= $random_number_pop)
                            {
                                last;
                            }
                            $ndonorsadd=$ndonorsadd+$nhapssourcevec[$k];
                        }

                                ## (II) PICK SIZE OF CHUNK:
                        $meanrate = $lambdavec[$h];
                        $random_number2 = rand(1);
                        $geneticdistCHUNK = 100*(log(1-$random_number2))/(-1.0*$meanrate);
                        $geneticlengthHAP = $geneticlengthHAP + $geneticdistCHUNK;

                        if ($nooverlap ne 'yes')
                        {
                               ## (III)-a PICK HAP TO COPY:
                            $random_number_hap = rand(1);
                            $sumprobs = 0.0;
                            for ($g=0; $g < $nhapssourcevec[$k]; $g+=1)
                            {
                                $sumprobs = $sumprobs + 1.0/$nhapssourcevec[$k];
                                if ($sumprobs >= $random_number_hap)
                                {
                                    last;
                                }
                            }
                            $hapval=0;
                            for ($l=0; $l < $k; $l+=1)
                            {
                                $hapval = $hapval+$nhapssourcevec[$l];
                            }
                            $hapval=$hapval+$g+1;
                        }
                        if ($nooverlap eq 'yes')
                        {
                             ## (III)-b FIND WHICH DONOR HAPS HAVE ALREADY BEEN COPIED AT THIS SNP:
                            $totalgeneticdistNEW=$totalgeneticdist;
                            @numtaken=();
                            @selectedhap=();
                            for ($a=0; $a < $numsources; $a+=1)
                            {
                                push(@numtaken,0);
                            }
                            for ($a=0; $a < $totaldonors; $a+=1)
                            {
                                push(@selectedhap,0);
                            }
                            for ($h = $start; $h < $nsites; $h+=1)
                            {
                                if ($totalgeneticdistNEW <= $geneticlengthHAP)
                                {
                                    for ($m=0; $m < $i; $m+=1)
                                    {
                                        if ($selectedhap[$indDONOR[$h][$j]]==0)
                                        {
                                            $selectedhap[$indDONOR[$h][$j]]=1;
                                            $sumdonors=0;
                                            for ($a=0; $a < $numsources; $a+=1)
                                            {
                                                if ($indDONOR[$m][$j] < ($sumdonors+$nhapssourcevec[$a]))
                                                {
                                                    $numtaken[$a]=$numtaken[$a]+1;
                                                    last;
                                                }
                                                $sumdonors=$sumdonors+$nhapssourcevec[$a];
                                            }
                                        }
                                    }
                                }
                                if ($totalgeneticdistNEW > $geneticlengthHAP)
                                {
                                    last;
                                }
                                if ($j < ($nsites - 1))
                                {
                                    $totalgeneticdistNEW = $totalgeneticdistNEW + $recomrateperbp[$j]*($posvec[($j+1)]-$posvec[$j])*100;
                                }
                            }
                            $outofdonorind=0;
                            for ($a=0; $a < $numsources; $a+=1)
                            {
                                if ($numtaken[$a] >= $nhapssourcevec[$a])
                                {
                                    print "OUT OF DONORS ( $numtaken[$a] $nhapssourcevec[$a] ) !!! -- hap $i snp $j pop $k \n";
                                    #print "@selectedhap\n";
                                    $outofdonorind=1;
                                    $outofdonorcount=$outofdonorcount+1;
                                    #die;
                                }
                            }

                                  ## (IV) PICK HAP TO COPY:
                            $exitind=0;
                            while($exitind==0)
                            {
                                $random_number_hap = rand(1);
                                $sumprobs = 0.0;
                                for ($g=0; $g < $nhapssourcevec[$k]; $g+=1)
                                {
                                    $sumprobs = $sumprobs + 1.0/$nhapssourcevec[$k];
                                    if ($sumprobs >= $random_number_hap)
                                    {
                                        last;
                                    }
                                }
                                if (($selectedhap[($ndonorsadd+$g)]==0) || ($outofdonorind==1))
                                {
                                    $exitind=1;
                                }
                            }
                            $hapval=$ndonorsadd+$g+1;
                        }
                        if ($i==0)
                        {
                            #print "$i $j $start $nsites $hapval $k $sourcevec[$k] $geneticlengthHAP $geneticdistCHUNK $totalgeneticdist\n";
                        }
                                    ## (V) PRINT OUT GENETIC MATERIAL:
                        for ($j = $start; $j < $nsites; $j+=1)
                        {
                            if ($totalgeneticdist <= $geneticlengthHAP)
                            {
                                if ($i==0 && $j < 1000)
                                {
                                    #print "$j $hapval ";
                                }
                                print OUT "$genomatDONOR[($hapval-1)][$j]";
                                #print OUT2 " $hapval";
                                print OUT2 " $allsourceSELECT[($hapval-1)]";
                                $indDONOR[$i][$j]=$hapval-1;
                            }
                            if ($totalgeneticdist > $geneticlengthHAP)
                            {
                                $start = $j;
                                last;
                            }
                            if ($j < ($nsites - 1))
                            {
                                $totalgeneticdist = $totalgeneticdist + $recomrateperbp[$j]*($posvec[($j+1)]-$posvec[$j])*100;
                                ##$totalgeneticdist = $totalgeneticdist + $recomrateperbp[$j]*($posvec[($j+1)]-$posvec[$j]);
                            }
                        }
                    }
                    print OUT "\n";
                    print OUT2 "\n";
                }
            }
        }
    }
    close(OUT);
    close(OUT2);
    if ($nooverlap eq 'yes')
    {
        print "unhappy: $outofdonorcount times.\n";
    }
    $command="gzip -f $outfile2";
    system ($command);
}
