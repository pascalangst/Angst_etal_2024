# get shared heterozygotes and presence of alternative alleles

# dimensions of VCF
nSites=$(bcftools query -f '[\t%SAMPLE=%GT]\n'  VCF.vcf | wc -l)

# for each sample
samples=$(bcftools query -l  VCF.vcf)

for sample in $samples; do 
	
  echo $sample > $sample.samplename
	vcftools --vcf VCF.vcf --keep $sample.samplename --recode --recode-INFO-all --out VCF.$sample

	nSamples=$(bcftools query -l VCF.$sample.recode.vcf | wc -l)

  # get missing info (count of how many samples are ./. or .|. per site), heterozygote sites (count how many samples have reads for both reference and alternative alleles per site), and number of samples (just add a column with the number of samples, same for each row)
	paste <( bcftools query -f '[\t%SAMPLE=%GT]\n' VCF.$sample.recode.vcf |\
    	awk 'BEGIN {print "nMiss"} {print gsub(/\.\|\.|\.\/\./, "")}') \
    	\
  		<( bcftools query -f '[\t%AD]\n' VCF.$sample.recode.vcf | sed -e 's/\.,\.//g' | sed -e 's/[0-9]*,0/0/g' | sed -E 's/([^0-9])0,[0-9]*/\12/g' | sed -e 's/[0-9]*,[0-9]*/1/g' |\
     	awk 'BEGIN {print "nHet"} {print gsub(/1/, "")}') \
     	\
     	<(for ((i=1; i<=$nSites; i++)); do echo -e "$nSamples"; done) | perl -s -p -e 's/nHet\t$number/nHet\tnSamples/' -- -number="$nSamples" > VCF.shared_het.$sample

done 
