# GenoMerge

Genomic dataset merge and format conversion tool. Uses PLINK 1.9 and PLINK 2.0.

## Dependencies

- Python >= 3.9 (standard library only)
- PLINK 1.9: https://www.cog-genomics.org/plink/1.9/
- PLINK 2.0: https://www.cog-genomics.org/plink/2.0/

Place `plink` and `plink2` on PATH or in the working directory.

Check tool paths:
```bash
genomerge --check-tools
```

## Installation

```bash
pip install .
```
Provides `genomerge` command.

## Input Formats

- PLINK binary (`.bed`, `.bim`, `.fam`)
- EIGENSTRAT (`.geno`, `.snp`, `.ind`)
- Raw DTC arrays: 23andMe, AncestryDNA, MyHeritage, FTDNA, LivingDNA (`.txt`, `.csv`, `.tsv`, `.zip`, `.gz`)
- VCF (`.vcf`, `.vcf.gz`)

## Output Formats

- EIGENSTRAT (default)
- PLINK binary (`--plink`)

## Usage

### Merge
```bash
# Merge datasets (outputs EIGENSTRAT)
genomerge --ref Data_Ref --sub Data1 --out merged

# Output PLINK binary
genomerge --ref Data_Ref --sub Data1 --plink --out merged

# Multiple cohorts
genomerge --ref Data_Ref --sub Data1 Data2 --out merged

# Batch file
genomerge --ref Data_Ref --batch list.txt --out merged

# Raw array or VCF input
genomerge --ref Data_Ref --raw raw_array.txt --pop Pop1 --id ID1 --out merged
genomerge --ref Data_Ref --vcf sample.vcf.gz --pop Pop1 --id ID1 --out merged

# Autosome filter (chr 1-22)
genomerge --ref Data_Ref --sub Data1 --autosome --out merged

# KING kinship filter
genomerge --ref Data_Ref --sub Data1 --king --out merged

# Modern LD pruning (50 10 0.2)
genomerge --ref Data_Ref --sub Data1 --ld-1 --out merged

# Ancient DNA (aDNA) LD pruning (200 25 0.4)
genomerge --ref Data_Ref --sub Data1 --ld-2 --out merged

# Exclude high-LD inversion regions (HLA, 8p23, 17q21, LCT)
genomerge --ref Data_Ref --sub Data1 --exclude-high-ld --out merged

# QC filters (MAF, missingness)
genomerge --ref Data_Ref --sub Data1 --maf 0.01 --geno 0.10 --out merged

# Extract / exclude specific SNPs from .snplist file
genomerge --ref Data_Ref --sub Data1 --extract snps_to_keep.snplist --out merged
genomerge --ref Data_Ref --sub Data1 --exclude bad_snps.snplist --out merged

# Dry run check
genomerge --ref Data_Ref --sub Data1 --check
```

### Convert (No Reference)
```bash
# PLINK to EIGENSTRAT
genomerge --convert Data1 --eig --out Data1_eig

# EIGENSTRAT to PLINK
genomerge --convert Data1 --plink --out Data1_plink

# Batch convert
genomerge --convert Data1 Data2 --plink --out Data1_p Data2_p

# Convert with aDNA LD pruning & autosome
genomerge --convert Data1 --plink --autosome --ld-2 --out Data1_clean
```

## Options

| Flag | Description |
| :--- | :--- |
| `--check-tools` | Check PLINK and PLINK2 paths |
| `--ref` | Reference dataset prefix |
| `--sub` | Subject dataset prefix(es) |
| `--convert` | Dataset(s) to convert without reference |
| `--out` | Output prefix(es) |
| `--eig` | Output EIGENSTRAT format |
| `--plink` | Output PLINK binary format |
| `--raw` | Raw DTC array file(s) |
| `--vcf` | VCF file(s) |
| `--batch` | File listing datasets, one per line |
| `--autosome` | Filter to chromosomes 1-22 |
| `--chr` | Filter to specified chromosomes |
| `--king` | KING kinship cutoff filter |
| `--ld-1` | Modern LD pruning preset (50 10 0.2) |
| `--ld-2` | Ancient DNA (aDNA) LD pruning preset (200 25 0.4) |
| `--ld` | Custom LD pruning (window step r2) |
| `--exclude-high-ld` | Exclude high-LD regions (HLA, 8p23, 17q21, LCT) |
| `--maf` | Filter by minor allele frequency |
| `--mac` | Filter by minor allele count |
| `--geno` | Filter variants by missing rate |
| `--mind` | Filter samples by missing rate |
| `--extract` | Extract variants from PLINK .snplist file |
| `--exclude` | Exclude variants from PLINK .snplist file |
| `--keep-pop` | Populations to keep |
| `--keep-ind` | Individuals to keep |
| `--check` | Dry run validation |
| `-v` | Verbose output |
