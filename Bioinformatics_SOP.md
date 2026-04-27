# GENOSAR Bioinformatics SOP


#### 1.  Study naming for isolates and sequences

Separate sections with dashes "-". 

1. Prefix: "GSAR"
2. Country code 
    * NG: Nigeria
    * PE: Peru
    * PK: Pakistan
    * ZA: South Africa
3. Species
    * Primary study:
    
    Abbv | Species
    --- | ---
    EC | <i>Escherichia coli</i>
    KP | <i>Klebsiella pneumoniae</i>

    * Additional:
    
    Abbv | Species
    --- | ---
    PA | <i>Pseudomonas aeruginosa</i>
    AB | <i>Acinetobacter baumannii</i> complex
    KA | <i>Klebsiella aerogenes</i>
    KO | <i>Klebsiella oxytoca</i>
    CL | <i>Enterobacter cloacae</i>
    TM | <i>Stenotrophomonas maltophilia</i>
    
    
4. Five-digit number (serial or random non-repeating)
5. <i>(Optional)</i> Extra identifier (maximum length 4, numbers [0-9] and/or letters [A-Z])
<br>

<strong>Example:</strong> GSAR-ZA-KP-01234
<strong>Example:</strong> GSAR-PK-EC-86420
<strong>Example with optional identifier:</strong> GSAR-PE-KP-13579-2B

<br>

#### 2.	MinKNOW basecalling and demultiplexing:

*  Run name: GSAR-XX-YYYY-MM-DD 
    *Replace "XX" with your site code and "YYYY-MM-DD" with the sequencing date. Add any optional modifiers after another dash character.*
*	Run limit: 72 hours
*	Pore scan frequency: 1.5 hrs
*	Basecalling: Super-accurate basecalling 400bps
*	Modified basecalling: Off
*	Trim barcodes: On
*	Min Q Score: 10
*	BAM file output: Off

<br>

#### 3.	Save report.html. Upload to REDCap
<br>

#### 4.  Create barcode list to REDCap as Excel or tab-separated text file. 

Must have at least two columns with first two columns being, in order, "barcode" and "id". Additional columns can be included, if wanted. 


<strong> Example: </strong>

*barcodes.txt*

barcode | id
--- | ---
barcode01 | GSAR-ZA-KP-01238
barcode02 | GSAR-ZA-KP-18293
barcode03 | GSAR-ZA-EC-28274

#### 5. Merge pass_filter reads by barcode id

In the run output directory create a directory named `fastq` and combine the demultiplexed read files in `fastq_pass` into separate files named for the isolate IDs. For example:
```Shell
cd <GSAR-XX-YYYY-MM-DD>/no_sample_id/<run_id>

mkdir -p fastq

while read bc id
do
    cat fastq_pass/${bc}/*.fastq.gz > fastq/${id}.fastq.gz
done < barcodes.txt
```

>[!NOTE] All subsequent steps will be repeated for each sequenced isolate
> Perform these steps in a data analysis directory separate from the MinKnow output directory

#### 6.  Set the analysis variables for your system. 

```Shell
threads=8
busco_db=`realpath path/to/enterobacteriaceae_odb12`
sourmash_db=`realpath /path/to/gammaproteobacteria.lca.json.gz`
mlst=`realpath /path/to/mlst_profiler.pl`
```

#### 7.  For each isolate sequence, set variables for the isolate name and path to the sequence reads.

```Shell
iso="GSAR-ZA-KP-01238" ## Example isolate name
reads=`realpath /path/to/GSAR-XX-YYYY-MM-DD.fastq.gz` ## Example path to reads
mkdir -p $iso
cd $iso
```
<br>

#### 8. Activate the first conda environment

```Shell
conda activate genosar_pt1_env
```
<br>


#### 9.	Use Nanoplot to assess read counts and quality

```Shell
nanoplot \
    -t $threads \
    -o nanoplot \
    --fastq_rich $reads 
```

Upload the `nanopore/NanoStats.txt` file.

<br>

#### 10.	Assemble the genome with autocycler

a. Subsamble nanopore reads

```Shell
autocycler subsample \
    --reads $reads \
    --out_dir 01_read_subsets \
    --genome_size 5500000 \
    --count 8 2> >(tee subsample.err)
``` 

b. Generate assemblies

```Shell
mkdir -p 02_assemblies
for subset in {1..8}
do
    subset=$(printf "%02d" $subset)
    for assembler in flye miniasm raven
    do
    echo "Running ${assembler} on read subset ${subset}"
    autocycler helper \
        ${assembler} \
        --reads 01_read_subsets/sample_${subset}.fastq \
        --out_prefix 02_assemblies/assembly_${subset}_${assembler} \
        --genome_size 5500000 \
        --read_type ont_r10 \
        --threads ${threads} \
        > >(tee 02_assemblies/assembly_${subset}_${assembler}.out) \
        2> >(tee 02_assemblies/assembly_${subset}_${assembler}.err)
    done
done
```

c. Reconcile assemblies

```Shell
autocycler compress \
    --assemblies_dir 02_assemblies \
    --autocycler_dir 03_autocycler \
    --threads ${threads}

autocycler cluster --autocycler_dir 03_autocycler

for c in 03_autocycler/clustering/qc_pass/cluster_*
do
    c_name=`basename $c`
    echo "$c_name"
    autocycler trim -c $c
    autocycler resolve -c $c
    BandageNG image $c/5_final.gfa $c/5_final.png
done

autocycler combine \
    --autocycler_dir 03_autocycler \
    --in_gfas 03_autocycler/clustering/qc_pass/cluster_*/5_final.gfa

autocycler table > metrics.tsv
autocycler table -a \. -n $iso >> metrics.tsv

BandageNG image \
    03_autocycler/consensus_assembly.gfa \
    03_autocycler/consensus_assembly.png
```

d. View the `03_autocycler/consensus_assembly.png` file to evaluate for assembly issues.

<br>

#### 11.	Medaka polishing and contig rotation

```Shell
medaka_consensus \
  -i $reads \
  -d 03_autocycler/consensus_assembly.fasta \
  -o 04_medaka \
  -t $threads \
  -f --bacteria  

mkdir -p 05_analysis
cd 05_analysis
cp ../04_medaka/consensus.fasta .
dnaapler all \
  -i consensus.fasta \
  -o dnaapler \
  -t $threads \
  -a nearest \
  -f
cp dnaapler/dnaapler_reoriented.fasta ${iso}.fasta
```

<br>

#### 12. Activate the second conda environment

```Shell
conda deactivate
conda activate genosar_pt2_env
```

<br>

#### 13. Quality control with BUSCO

Use BUSCO to assess completeness and contamination levels.

 > :bulb: To avoid having to download the database file (147 Mb) every time you run this command, create a `busco_database` directory. In the database directory run the following command: `busco --download enterobacteriaceae_odb12`. The database will be downloaded into `busco_downloads/lineages/enterobacteriaceae_odb12/`. Give the path to the downloaded database (set in **Step 6** above) to the `-l` option in the command below.

```Shell
busco -i ${iso}.fasta \
    -m genome \
    -l ${busco_db} \
    --offline \
    -o busco \
    -c $threads
```

<br>

#### 14.  Species identification

Using sourmash seems to be one of the quickest and most straighforward ways to do rough species identification programatically.

We'll use a [custom Gammaproteobacteria database](sourmash_db.md) generated from NCBI reference sequences and override the default threshold of 5 to improve the specificity of the species identification.

```Shell
sourmash sketch dna \
    --name $iso \
    -o ${iso}.fasta.sig \
    ${iso}.fasta
sourmash lca classify \
    --db ${sourmash_db} \
    --threshold 10 \
    --query ${iso}.fasta.sig > sourmash.txt
```

#### 15. Sequence typing

Will provide profiles and allele sequences from Pasteur (K pneumo) / PubMLST (E coli)
For the `-o` option, use "KP" for K pneumonaie and "EC" for E coli.

```Shell
${mlst} \
    -o KP \ 
    -f ${iso}.fasta > mlst.txt
```

#### 16. AMRFinder

Options for `-O` are "Klebsiella_pneumoniae" or "Escherichia"

```
amrfinder \
    -n ${iso}.fasta \
    --name $iso \
    -O Klebsiella_pneumoniae \
    -o amrfinder.txt 
```