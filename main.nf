#!/usr/bin/env nextflow

VERSION="0.30"

log.info "===================================================================="
log.info "GATK4 Best Practice Nextflow Pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run Biocentric/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz"
  log.info " "
  log.info "Mandatory arguments:"
  log.info "    --fastq1        FILE               Fastq(.gz) file for read1"
  log.info "    --fastq2        FILE               Fastq(.gz) file for read2"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --outdir        DIR                Output directory(default: ./Results)"
  log.info "    --samplename    STRING             Sample name(dafault: fastq1 basename)"
  log.info "    --rg            STRING             Read group tag(dafault: fastq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}


fastq1 = file("$params.fastq1")
fastq2 = file("$params.fastq2")
params.outdir = "./Results"
params.samplename = fastq1.baseName
params.rg = fastq1.baseName

process get_reference {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'
	
	output:
	file "Homo_sapiens_assembly38.fasta" into reference
	file "Homo_sapiens_assembly38.dict" into reference_dict
	file "Homo_sapiens_assembly38.fasta.fai" into reference_fai
	
	"""
	gunzip -dc /data/Homo_sapiens_assembly38.fasta.gz > Homo_sapiens_assembly38.fasta
  gunzip -dc /data/Homo_sapiens_assembly38.dict.gz > Homo_sapiens_assembly38.dict
  gunzip -dc /data/Homo_sapiens_assembly38.fasta.fai.gz > Homo_sapiens_assembly38.fasta.fai
	"""
}

process get_dbSNP {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'

	output:
	file "Homo_sapiens_assembly38.dbsnp138.vcf" into dbsnp
	file "Homo_sapiens_assembly38.dbsnp138.vcf.idx" into dbsnp_idx

	"""
	gunzip -dc /data/Homo_sapiens_assembly38.dbsnp138.vcf.gz > Homo_sapiens_assembly38.dbsnp138.vcf
  gunzip -dc /data/Homo_sapiens_assembly38.dbsnp138.vcf.idx.gz > Homo_sapiens_assembly38.dbsnp138.vcf.idx
	"""
}

process get_golden_indel {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'

	output:
	file "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" into golden_indel
	file "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" into golden_indel_idx

	"""
	cp /data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ./Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  gunzip -dc /data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi.gz > Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
	"""
}

process get_hapmap {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'

	output:
	file "hapmap_3.3.hg38.vcf.gz" into hapmap
	file "hapmap_3.3.hg38.vcf.gz.tbi" into hapmap_idx

	"""
	cp /data/hapmap_3.3.hg38.vcf.gz ./hapmap_3.3.hg38.vcf.gz
  gunzip -dc /data/hapmap_3.3.hg38.vcf.gz.tbi.gz > hapmap_3.3.hg38.vcf.gz.tbi
	"""
}

process get_omni {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'

	output:
	file "1000G_omni2.5.hg38.vcf.gz" into omni
	file "1000G_omni2.5.hg38.vcf.gz.tbi" into omni_idx

	"""
	cp /data/1000G_omni2.5.hg38.vcf.gz ./1000G_omni2.5.hg38.vcf.gz
  gunzip -dc /data/1000G_omni2.5.hg38.vcf.gz.tbi > 1000G_omni2.5.hg38.vcf.gz.tbi
	"""
}

process get_phase1_SNPs {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'

	output:
	file "1000G_phase1.snps.high_confidence.hg38.vcf.gz" into phase1_snps
	file "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" into phase1_snps_idx

	"""
	cp /data/1000G_phase1.snps.high_confidence.hg38.vcf.gz ./1000G_phase1.snps.high_confidence.hg38.vcf.gz
  gunzip -dc /data/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi.gz > 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
	"""
}

process get_BWA_index {
	publishDir "${params.outdir}/reference"
	container 'biocentric/hg38-gatk4'

	output:
	set "Homo_sapiens_assembly38.fasta.64.alt", "Homo_sapiens_assembly38.fasta.64.amb", "Homo_sapiens_assembly38.fasta.64.ann", "Homo_sapiens_assembly38.fasta.64.bwt", "Homo_sapiens_assembly38.fasta.64.pac", "Homo_sapiens_assembly38.fasta.64.sa" into bwa_index

	"""
	gunzip -dc /data/Homo_sapiens_assembly38.fasta.64.alt.gz > Homo_sapiens_assembly38.fasta.64.alt
  gunzip -dc /data/Homo_sapiens_assembly38.fasta.64.amb.gz > Homo_sapiens_assembly38.fasta.64.amb
  gunzip -dc /data/Homo_sapiens_assembly38.fasta.64.ann.gz > Homo_sapiens_assembly38.fasta.64.ann
  gunzip -dc /data/Homo_sapiens_assembly38.fasta.64.bwt.gz > Homo_sapiens_assembly38.fasta.64.bwt
  gunzip -dc /data/Homo_sapiens_assembly38.fasta.64.pac.gz > Homo_sapiens_assembly38.fasta.64.pac
  gunzip -dc /data/Homo_sapiens_assembly38.fasta.64.sa.gz > Homo_sapiens_assembly38.fasta.64.sa
	"""
}

process BWA {
	publishDir "${params.outdir}/MappedRead"
	container 'dukegcb/bwa-samtools'

	input:
	file reference
	file bwa_index
	file fastq1
	file fastq2

	output:
	file 'aln-pe.bam' into samfile
	
	"""
	bwa mem -M -t 19 -R '@RG\\tID:${params.rg}\\tSM:${params.samplename}\\tPL:Illumina' $reference $fastq1 $fastq2 | samtools view -bS - > aln-pe.bam
	"""
		
}

process BWA_sort {
	publishDir "${params.outdir}/MappedRead"
	container 'mgibio/samtools:1.9'
	
	input:
	file samfile

	output:
	file 'aln-pe-sorted.bam' into bam_sort

	"""
	samtools sort -T temp -@ 10 -m 2G -o aln-pe-sorted.bam -O BAM $samfile
	"""

}

process MarkDuplicates {
	publishDir "${params.outdir}/MappedRead"
	container 'broadinstitute/gatk'
	
	input:
	file bam_sort

	output:
	file 'aln-pe_MarkDup.bam' into bam_markdup

	"""
	gatk MarkDuplicates -I $bam_sort -M metrics.txt --TMP_DIR=${params.outdir}/tmp -O aln-pe_MarkDup.bam	
	"""

}

process BaseRecalibrator {
	publishDir "${params.outdir}/BaseRecalibrator"
	container 'broadinstitute/gatk:latest'
	
	input:
	file reference
	file reference_fai
	file reference_dict
	file bam_markdup
	file dbsnp
	file dbsnp_idx
	file golden_indel
	file golden_indel_idx

	output:
	file 'recal_data.table' into BaseRecalibrator_table

	"""
	gatk BaseRecalibrator \
	-I $bam_markdup \
	--known-sites $dbsnp \
	--known-sites $golden_indel \
	-O recal_data.table \
	-R $reference
	"""
}

process ApplyBQSR {
	publishDir "${params.outdir}/BaseRecalibrator"
	container 'broadinstitute/gatk:latest'
	
	input:
	file BaseRecalibrator_table
	file bam_markdup

	output:
	file 'aln-pe_bqsr.bam' into bam_bqsr
	
	script:
	"""
	gatk ApplyBQSR -I $bam_markdup -bqsr $BaseRecalibrator_table -O aln-pe_bqsr.bam
	"""
}

process HaplotypeCaller {
	publishDir "${params.outdir}/HaplotypeCaller"
	container 'broadinstitute/gatk:latest'
	
	input:
	file reference
	file reference_fai
	file reference_dict
	file bam_bqsr

	output:
	file 'haplotypecaller.g.vcf' into haplotypecaller_gvcf
	
	script:
	"""
	gatk HaplotypeCaller -I $bam_bqsr -O haplotypecaller.g.vcf --emit-ref-confidence GVCF -R $reference
	"""
}

process HaplotypeGVCF {
	publishDir "${params.outdir}/HaplotypeCaller"
	container 'broadinstitute/gatk:latest'
	
	input:
	file reference
	file reference_fai
	file reference_dict
	file bam_bqsr

	output:
		
	script:
	"""
	gatk HaplotypeCaller -I $bam_bqsr -O haplotypecaller.gvcf --emit-ref-confidence BP_RESOLUTION -R $reference
	"""
}

process GenotypeGVCFs {
	publishDir "${params.outdir}/HaplotypeCaller"
	container 'broadinstitute/gatk:latest'
	
	input:
	file reference
	file reference_fai
	file reference_dict
	file haplotypecaller_gvcf

	output:
	file 'haplotypecaller.vcf' into haplotypecaller_vcf
	
	script:
	"""
	gatk GenotypeGVCFs --variant haplotypecaller.g.vcf -R $reference -O haplotypecaller.vcf
	"""
}



process VariantRecalibrator_SNPs {
	publishDir "${params.outdir}/VariantRecalibrator"
	container 'broadinstitute/gatk:latest'
	
	input:
	file reference
	file reference_fai
	file reference_dict
	file haplotypecaller_vcf
	file hapmap
	file hapmap_idx
	file omni
	file omni_idx
	file phase1_snps
	file phase1_snps_idx
	file dbsnp
	file dbsnp_idx

	output:
	file 'recalibrate_SNP.recal' into variantrecalibrator_recal
	file 'recalibrate_SNP.recal.idx' into variantrecalibrator_recal_idx
	file 'recalibrate_SNP.tranches' into variantrecalibrator_tranches

	script:
	"""
	gatk VariantRecalibrator \
	-V $haplotypecaller_vcf \
 	-R $reference \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./$hapmap \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 ./$omni \
    	--resource:1000G,known=false,training=true,truth=false,prior=10.0 ./$phase1_snps \
    	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./$dbsnp \
	-an DP \
    	-an QD \
	-an FS \
    	-an SOR \
    	-an MQ \
    	-an MQRankSum \
    	-an ReadPosRankSum \
    	-mode SNP \
    	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	--max-gaussians 8 \
    	-O recalibrate_SNP.recal \
    	--tranches-file recalibrate_SNP.tranches \
	"""
}



process ApplyVQSR_SNPs {
	publishDir "${params.outdir}/VariantRecalibrator"
	container 'broadinstitute/gatk:latest'
	
	input:
	file haplotypecaller_vcf
	file variantrecalibrator_recal
	file variantrecalibrator_recal_idx
	file variantrecalibrator_tranches

	output:
	file 'recalibrated_snps_raw_indels.vcf' into recalibrated_snps_raw_indels
	
	script:
	"""
	gatk ApplyVQSR \
	-V $haplotypecaller_vcf \
	--recal-file $variantrecalibrator_recal \
	--tranches-file $variantrecalibrator_tranches \
	-mode SNP \
	-ts-filter-level 99.0 \
	-O recalibrated_snps_raw_indels.vcf 
	"""
}



process VariantRecalibrator_INDELs {
	publishDir "${params.outdir}/VariantRecalibrator"
	container 'broadinstitute/gatk:latest'
	
	input:
	file reference
	file reference_fai
	file reference_dict
	file recalibrated_snps_raw_indels
	file dbsnp
	file dbsnp_idx
	file golden_indel
	file golden_indel_idx

	output:
	file 'recalibrate_INDEL.recal' into variantrecalibrator_indel_recal
	file 'recalibrate_INDEL.recal.idx' into variantrecalibrator_indel_recal_idx
	file 'recalibrate_INDEL.tranches' into variantrecalibrator_indel_tranches
	
	script:
	"""
	gatk VariantRecalibrator \
	-V $recalibrated_snps_raw_indels \
 	-R $reference \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 ./$golden_indel \
    	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./$dbsnp \
	-an QD \
    	-an DP \
    	-an FS \
	-an SOR \
    	-an MQRankSum \
    	-an ReadPosRankSum \
    	-mode INDEL \
    	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	--max-gaussians 4 \
    	-O recalibrate_INDEL.recal \
    	--tranches-file recalibrate_INDEL.tranches \
	"""
}

process ApplyVQSR_INDELs {
	publishDir "${params.outdir}/VariantRecalibrator"
	container 'broadinstitute/gatk:latest'
	
	input:
	file recalibrated_snps_raw_indels
	file variantrecalibrator_indel_recal
	file variantrecalibrator_indel_recal_idx
	file variantrecalibrator_indel_tranches

	output:
	file 'recalibrated_variants.vcf' into recalibrated_variants_vcf
	
	script:
	"""
	gatk ApplyVQSR \
	-V $recalibrated_snps_raw_indels \
	--recal-file $variantrecalibrator_indel_recal \
	--tranches-file $variantrecalibrator_indel_tranches \
	-mode INDEL \
	-ts-filter-level 99.0 \
	-O recalibrated_variants.vcf
	"""
}

process copy {
	publishDir "${params.outdir}", mode: 'copy'

	input:
	file recalibrated_variants_vcf

	output:
	file "${params.samplename}.vcf"
	
	"""
	mv $recalibrated_variants_vcf ${params.samplename}.vcf
	"""

}

