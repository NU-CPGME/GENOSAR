# Custom Sourmash LCA database

## 1. Download reference genome sequences

Requires NCBI datasets command line tools. 
Conda package: `conda-forge::ncbi-datasets-cli`
Version used: 18.23.0

Gammaproteobacteria class bacterial genomes (taxid:1236), reference genomes only

```Shell
datasets download genome taxon 1236 --reference --dehydrated
unzip ncbi_dataset.zip -d gammaproteobacteria
datasets rehydrate --directory gammaproteobacteria/
```

## 2. Gather all genome sequences into a single directory

Symlink and rename the sequence files 

```Shell
mkdir gammaproteobacteria_renamed
cd gammaproteobacteria_renamed

i=0
for file in ../gammaproteobacteria/ncbi_dataset/data/*/*.fna
do
  ln -s $file $i.fna
  i=$(($i+1))
done
```

## 3. Use sourmash to create sketches of the reference genome sequences

Conda package: `bioconda::sourmash`
Version used: 4.9.4

```Shell
sourmash sketch dna --name-from-first *.fna
```

## 4. Download the taxonomy list using NCBI datasets

```Shell
cd ..
datasets summary taxonomy taxon 1236 --children --as-json-lines | \
    dataformat tsv taxonomy --template tax-summary > gammaproteobacteria.taxonomy.txt
```

## 5. Reformat the taxonomy list

Use the custom script in the `scripts` directory on Github to convert the taxonomy list and filter it to include only the downloaded genomes with unambiguous lineages.

It's a bit hacky and will omit several genomes, but these are mostly unnamed and rare species likely highly divergent from anything that will end up sequenced in this study. Ommitted genome sequences will be listed in the output as the program runs. 

```Shell
perl ~/Scripts/taxonomy_to_sourmash.pl \
    gammaproteobacteria.taxonomy.txt \
    gammaproteobacteria_renamed > gammaproteobacteria_lineages.csv
```

## 6. Create the sourmash LCA database

```Shell
sourmash lca index \
    gammaproteobacteria_lineages.csv \
    gammaproteobacteria.lca.json \
    gammaproteobacteria_renamed/*.sig \
    --split-identifiers
gzip gammaproteobacteria.lca.json
```
