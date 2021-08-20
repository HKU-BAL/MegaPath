### LSAM format for MegaPath Internal Use

1. Annotated FASTQ
2. LSAM
3. LSAM.id
4. Summary

```
Annotated FASTQ == (fastq2lsam) ==> LSAM
m8 == (m8_to_lsam) ==> LSAM
LSAM, LSAM == (r2c_to_r2g) ==> LSAM
LSAM == (extractFromLSAM) ==> FASTQ
LSAM == (taxonomyLookup) ==> LSAM.id
LSAM.id == (rassign) ==> LSAM.id
LSAM.id == (genCountTable) ==> cnt[,ann]
```

TODO:
- fastq2lsam: support optional tags [x]
- taxonomyLookup: recalculate maximum score [x] 
- reassign [x], translate: support LSAM.id format