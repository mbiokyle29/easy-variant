<h1>EasyVariant.pl</h1>

This 'script' is a miniature, basic variant caller, it processes a sam file containing matched alignments
Each SAM read is formed into a Moose object, containing Sequence objects and Cigar objects
The SAM is then manipulated, and the Cigar string is parsed to determine what exactly the alignment determined
The number of A,T,C, and G is tracked for each base position, along with X for deletions and I for insertions
An insertion at what would be the 'new' base 10 is counted at the old base 9 (the base preceding the insertion)

The optional arguments allow for a decent amount of fine tuning, including:
 -choosing specific ranges to only include
 -ignoreing a set of ranges
 -designating a set of ranges as repeat regions (For coverage statistics)
 -getting not only the variant calls, but the ATCGXI tally for each call

Notes|Warnings|Disclaimer

 	Using the --start --end flags, with the other range functions may cause unintended behvaior
	best to only use one (start/end) or (ignore and repeat)

	Be sure to set use lib to the correct dir for the Moose modules

	Range strings must be entered as int-int,int-int !!!
	
	The Genome package is backed by mysql, and might take some configuration
