YOC
===

YOC is a whole genome alignment tool designed for collinear genomes.
It has been developed by Raluca Uricaru (LIRMM/LaBRI).

Usage
-----

The main script to run is `runYOC.csh`, which takes:

1. the first genome to compare in fasta format
2. the second genome to compare in fasta format
3. the E-value for YASS similarity search
4. the overlap ratio for `OverlapChainer`

Use case example:

```
./runYOC.csh fasta_files/CP000412_GR.fsa fasta_files/CR954253_GR.fsa 10 10
```

The results are generated in two folders:

`yass_output`
: YASS intermediary results in `.yass` and `.mat` format

`results`
: chaining result (`.chn`) and coverage statistics


Installation
------------

1. Compile the translator : `make yass_blast2chainer`
2. Compile the java code (from the root dir): `javac -sourcepath OverlapChainer/src -d OverlapChainer/bin OverlapChainer/src/overlapChainer/*.java OverlapChainer/src/tests/*.java`
3. Get Yass v1.14: <http://bioinfo.lifl.fr/yass/download.php> and link it (or copy it) into the root folder under the name: `yass_1.14` 
5. Enjoy ;)
