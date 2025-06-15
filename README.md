# GPU-Accelerated FASTQ to BAM Pipeline using Parabricks

This README provides instructions to execute a GPU-accelerated variant of the GATK Best Practices alignment pipeline using NVIDIA Parabricks' `fq2bam` tool.

## 📦 Description

This pipeline takes paired-end FASTQ files and performs:

* Alignment with BWA-MEM
* Sorting by coordinate
* Marking duplicates
* Base Quality Score Recalibration (BQSR)

All operations are GPU-accelerated using the `cosmos9526/fq2bam_ref` Docker image.

---

## 🚀 Run Command

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

---

## 🧠 Notes

* `--low-memory` reduces GPU memory usage to \~16GB per GPU.
* `--num-gpus 2` ensures both specified GPUs are used.
* Input and reference files must exist in the mounted volumes.

---

## 📁 Expected Output

After successful execution, the following files will be present in `/home/dgx/0MMilad/output`:

* `output.bam` — Aligned and processed BAM file
* `output.bam.bai` — BAM index file (auto-generated)
* `recal.txt` — Base recalibration report

---

## ✅ Post-run Validation

You can validate the output using `samtools`:

### View BAM content:

```bash
samtools view output.bam | head
```

### Check BAM statistics:

```bash
samtools flagstat output.bam
```

---

## 🧪 Optional Enhancements

* Add `--log-file /outputdir/fq2bam.log` to store logs.
* Reduce `--bwa-options "-K 500000"` to limit memory per block.

---

## 🧬 Reference Genome

Ensure `/refrence_genome/hg38.fa` and its accompanying index files exist in the container.
You can generate `.fai` index with:

```bash
samtools faidx hg38.fa
```

---

## 🧩 Version Info

This example uses:

* Parabricks version: `4.3.2-1`
* Docker image: `cosmos9526/fq2bam_ref:tagname`
* CUDA version: >= 11.0

---

## 👨‍🔬 Author

Milad Bagheri — GPU-accelerated genomics workflow setup.

---

## 📚 Documentation

* Parabricks Docs: [https://docs.nvidia.com/clara/#parabricks](https://docs.nvidia.com/clara/#parabricks)
