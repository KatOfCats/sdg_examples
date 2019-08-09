## Triploid genome analysis (F. cylindrus)

### Contig assembly and content analysis

We create a contigs_raw.gfa output using w2rap-contigger.

```shell
w2rap-contigger/bin/w2rap-contigger -t 16 -r "readfile_1.fastq, readfile_2.fastq" -o F.cylindrus_assembly -p F.cylindrus
```

This is split into a .gfa and accompanying .fasta using custom python scripts. Verification of read content in the assembly is performed by KAT k-mer spectra.

```shell
kat comp -o F.cylindrus_contigs_vs_reads -t8 -H 2000000000 -I 1000000000 "readfile_1.fastq readfile_2.fastq" F.cylindrus_contigs_raw.fasta
```

Total genomic content is estimated from the spectra matrix file.
$$
\frac{total\; shared\; k-mers}{k-mer\; coverage}
$$
K-mer coverage is obtained from the mode of the fundamental frequency distribution of the spectra. 

Bacterial contamination is isolated by k-mer frequency in a KAT sect. This contamination can be identified by BLAST, filtered for 95% coverage and 95% identity. The main F. cylindrus distribution can also be confirmed by BLAST. 

```shell
kat sect -t8 -H 2000000000 F.cylindrus_contigs_raw.fasta readfile_1.fastq readfile_2.fastq

cat kat-sect-stats.tsv |awk '{if ($2 < 25 || $2 > 220) print;}' > bacterial_isolate.fasta
```



### Creating the WorkSpace

We first create a datastrore of the short reads, long reads and k-mer counter. 

```shell
sdg-datastore make -t paired -n illumina -o short-reads.datastore -1 readfile_1.fastq -2 readfile_2.fastq

sdg-datastore make -t long -n nanopore -o long-reads.datastore -L longread_file.fastq

sdg-kmercounts make -g F.cylindrus_contigs_raw.gfa -f readfile_1.fastq -f readfile_2.fastq -o k-mer-count -n kmer-count
```

We then use SDG to create a workspace from the W2rap-contigger .gfa, paired-end reads, long reads and k-mer count. 

```python
import sys
sys.path.append('path/to/sdg/build')
import pysdg as SDG
ws = SDG.WorkSpace()
ws.sdg.load_from_gfa("F.cylindrus_contigs_raw.gfa")
ws.add_paired_reads_datastore("short-reads.datastore.prseq", "PE1")
ws.add_long_reads_datastore("long-reads.datastore.loseq","LR1")
ws.add_kmer_counter("pe.count","kmer-count")
kcount = ws.get_kmer_counter("kmer-count")
```

And check the datastores and k-mer counter are available before saving the datastore to disk.

```
ws.list_paired_reads_datastores()
ws.list_long_reads_datastores()
ws.list_kmer_counters()
kcount.list_names()
ws.dump_to_disk('sr-lr-kc.sdgws')
```



### Analysing bubble content to check for ploidy signatures

Check for bubbles and how their k-mers are distirbuted. This will give an indication of the number of haplotypes in syntenic regions of the graph. 

```python
total number of bubbles code
turn this into a figure containing a number of graphs for both parallel sides
```

Graph structure can be viewed using BANDAGE. We can trace three potential paths through this fragmented region.

![DBG-double-bubble](/Users/kathodgkinson/Documents/Documents/Publication/DBG-double-bubble.png)

### Long read mapping



```
lords = ws.get_long_reads_datastore('LR1')
```

