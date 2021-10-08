vcfPath=$1
outVCF=$2
echo /nfs/seqscratch10/am5153/bin/TEST/SURVIVOR-1.0.7/Debug/SURVIVOR merge ${vcfPath} 1000 1 1 0 0 0 ${outVCF}.vcf
/nfs/seqscratch10/am5153/bin/TEST/SURVIVOR-1.0.7/Debug/SURVIVOR merge ${vcfPath} 1000 1 1 0 0 0 ${outVCF}.vcf && vcf-sort ${outVCF}.vcf > ${outVCF}.tmp.vcf && mv ${outVCF}.tmp.vcf ${outVCF}.vcf && bgzip -c ${outVCF}.vcf > ${outVCF}.vcf.gz && tabix -p vcf ${outVCF}.vcf.gz 