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
: YASS intermediary results in `.yass`, `.yass1`, and `.mat` formats

`results`
: simple chaining result with coordinates only (`.chn`), chaining result with alignments (`.chn_yass`), and coverage statistics


Installation
------------

1. Compile the translator : `make yass_blast2chainer`
2. Compile the java code (from the root dir): `javac -sourcepath OverlapChainer/src -d OverlapChainer/bin OverlapChainer/src/overlapChainer/*.java OverlapChainer/src/tests/*.java`
3. Get Yass v1.14: <http://bioinfo.lifl.fr/yass/download.php> and link it (or copy it) into the root folder under the name: `yass_1.14` 
5. Enjoy ;)


Requirements
------------

To run YOC,
- any version of tcsh, awk, Perl and a C compiler
- Java at least 1.6 (should work with minor adaptations downto 1.4)

are required.

YOC works on Apple Mac OS X and Linux systems.


Output
------

YOC produces several result files, as follows:

- .yass file (generated by YASS with -d 2 option), which gives information for every local similarity fragment (e.g., coordinates, length, identity percentage, alignment length, mismatches, gaps openings, E-value, ...). See YASS manual page for more details.

- .yass1 file (generated by YASS with -d 1 option), which gives the alignment for each fragment. See YASS manual page for more details.

- .mat file (obtained from the .yass file) gives the number of exact matches in each fragment and its coordinates.

- .chn file (generated by OverlapChainer) gives those fragments that were taken in the chain, among the ones contained in the .mat file (the .chn format is similar to .mat format).

- .chn_yass file (obtained from the .chn and the .yass1 files) gives the same fragments as in the .chn file, together with their alignments (the .chn_yass format is similar to .yass1 format).