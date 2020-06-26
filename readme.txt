Requirements:
bam_readcount
samtools

python
pandas
numpy

R 
RVAideMemoire

Also requires a directory named transcripts with "transcript"".fa within.
And a transcript.txt file with "transcript" names listed.


The test.sh script is a sort of wrapper script that specifies the threads and 
calls the actual script adenosine-m6A-KO1-P1.sh script for the actual analysis.

If all is complete there should be a directory containing intermediate steps for each
step of the pipeline concluding with the gTest directory.

Run code:
./test.sh -r transcripts.txt -c safe/Ad5_transcript.M3KO1-alt.sorted.bam -t safe/Ad5_transcript.M3P1-alt.sorted.bam


Things to add.
1. Incorporate find_candidates.py and finally genomic_locations.py scripts into pipeline.
2. Add arguments for threads
3. Fix bugs with odds ratio test
4. Fix visualization issues (per candidate and whole transcript)
5. Optimize code and add comments