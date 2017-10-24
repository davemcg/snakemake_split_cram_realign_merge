# snakemake_split_cram_realign_merge

Code outlining issue with getting Snakemake to split a cram file by read groups (each cram file potentially has unique set of read groups), re-aligning each individual lane cram, and then merging the lane crams into one file. 

So going from one to many to one. 

Discussion and trouble-shooting here:
https://stackoverflow.com/questions/46714560/snakemake-how-do-i-use-a-function-that-takes-in-a-wildcard-and-returns-a-value/46916902#46916902

