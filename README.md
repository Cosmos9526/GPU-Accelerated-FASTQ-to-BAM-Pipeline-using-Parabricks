# GPU-Accelerated FASTQ to VCF Pipeline using Parabricks

This README provides a step-by-step guide to execute a **GPU-accelerated GATK Best Practices pipeline** using NVIDIA Parabricks. The workflow includes both alignment (FASTQ to BAM) and variant calling (BAM to VCF), leveraging the high-performance `cosmos9526/fq2bam_ref` Docker image.

---

## 📦 Pipeline Overview

### Step 1: FASTQ to BAM (`fq2bam`)

Performs:

* Alignment with BWA-MEM
* Sorting by coordinate
* Marking duplicates
* Base Quality Score Recalibration (BQSR)

### Step 2: BAM to VCF (`haplotypecaller`)

Performs:

* Local reassembly of active regions
* SNP and indel calling
* Genotype likelihood calculation
* Variant calling in VCF format

All steps are fully GPU-accelerated.

---

## 🚀 Run Commands

### 1️⃣ FASTQ to BAM

```bash
docker run --rm -it \
  --runtime=nvidia \
  --gpus '"device=1,2"' \
  --privileged \
  -v /home/dgx/0MMilad/input_example:/workdir \
  -v /home/dgx/0MMilad/output:/outputdir \
  cosmos9526/fq2bam_ref:tagname \
  pbrun fq2bam \
    --ref /refrence_genome/hg38.fa \
    --in-fq /workdir/normal_R1.fastq.gz /workdir/normal_R2.fastq.gz \
    --knownSites /refrence_genome/dbsnp_138.hg38.vcf.gz \
    --knownSites /refrence_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites /refrence_genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --out-bam /outputdir/output.bam \
    --out-recal-file /outputdir/recal.txt \
    --num-gpus 2 \
    --low-memory
```

### 2️⃣ BAM to VCF

```bash
docker run --rm -it \
  --runtime=nvidia \
  --gpus '"device=1,2"' \
  --privileged \
  -v /home/dgx/0MMilad/input_example:/workdir \
  -v /home/dgx/0MMilad/output:/outputdir \
  cosmos9526/fq2bam_ref:tagname \
  pbrun haplotypecaller \
    --ref /refrence_genome/hg38.fa \
    --in-bam /outputdir/output.bam \
    --out-variants /outputdir/variants.vcf \
    --num-gpus 2
```

---

## 📁 Input Requirements

| File                                | Description                  |
| ----------------------------------- | ---------------------------- |
| `normal_R1.fastq.gz`, `R2.fastq.gz` | Paired-end FASTQ files       |
| `hg38.fa`                           | Reference genome FASTA       |
| `*.vcf.gz`                          | Known variant sites for BQSR |

> Ensure all reference and known sites are properly indexed and available in the mounted volume.

---

## 📤 Output Files

| File             | Description                     |
| ---------------- | ------------------------------- |
| `output.bam`     | Aligned BAM file after BQSR     |
| `output.bam.bai` | BAM index file (auto-generated) |
| `recal.txt`      | BQSR recalibration report       |
| `variants.vcf`   | Called variants in VCF format   |

---

## 🧪 Validation

### View BAM:

```bash
samtools view output.bam | head
```

### Stats:

```bash
samtools flagstat output.bam
```

### View VCF:

```bash
less /home/dgx/0MMilad/output/variants.vcf
```

### Count variants:

```bash
grep -v "^#" /home/dgx/0MMilad/output/variants.vcf | wc -l
```

---

## 🧠 Tips & Enhancements

* Use `--log-file /outputdir/fq2bam.log` to save logs.
* Tune memory with `--bwa-options "-K 500000"`.
* Confirm BAM index (`output.bam.bai`) is present before variant calling.

---

## 🧬 Software Versions

* **Parabricks**: `4.3.2-1`
* **Docker Image**: `cosmos9526/fq2bam_ref:tagname`
* **CUDA**: `>= 11.0`

---

## 👨‍🔬 Author

Milad Bagheri — GPU-powered genomics workflow designer

---

## 📚 References

* NVIDIA Parabricks: [https://docs.nvidia.com/clara/#parabricks](https://docs.nvidia.com/clara/#parabricks)
* GATK HaplotypeCaller: [https://gatk.broadinstitute.org](https://gatk.broadinstitute.org)
