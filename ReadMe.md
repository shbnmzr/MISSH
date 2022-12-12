# MISSH : Efficient Hashing of Multiple Spaced Seeds with Application

Alignment-free analysis of sequences has enabled high-throughput processing of sequencing data in many bioinformatics pipelines. 
Hashing k-mers is a common function across many alignment-free applications and it is widely used for indexing, querying and rapid similarity search. Recently, spaced seeds, a special type of pattern that accounts for errors or mutations, are routinely used instead of k-mers. Spaced seeds allow to improve the sensitivity, with respect to k-mers, in many applications, however the hashing of spaced seeds increases substantially the computational time. Moreover, if multiple spaced seeds are used the accuracy can further increases at the cost of running time. 
In this paper we address the problem of efficient multiple spaced seed hashing. The proposed algorithms exploit the similarity of adjacent spaced seed hash values in an input sequence in order to efficiently compute the next hashes. 
We report the results on several tests which show that our methods
significantly outperform the previously proposed algorithms, with a speedup that can reach 20x. 
We also apply these efficient spaced seeds hashing algorithms to an application in the field of metagenomic, the classification of reads performed by Clark-S, and we shown that a significant speedup can be obtained, thus resolving the slowdown introduced by the use of multiple spaced seeds.

----------------------------------------------------------------------------------------------------------------------------------------------------



The purpose of this program is to benchmark different approaches that solve the problem of computing the hashing of DNA sequences by using spaced seeds.
The original code was written by Samuele Girotto for the paper that presented the FSH approach (paper: http://www.dei.unipd.it/~ciompin/main/fsh.html repo: https://bitbucket.org/samu661/fsh/src/master/ ), was later modified by Enrico Petrucci by adding the implementation of the ISSH, first for the version single seed(paper and repo: http://www.dei.unipd.it/~ciompin/main/issh.html ), and then -
together with Eleonora Miani - by adding the version supporting multiple spaced seeds at the same time.

The datasets used for testing are available at: https://bitbucket.org/samu661/metaprob/src/master/


To compile the program open the terminal inside the main folder and use:

```sh
$ make all
```

To run the program call it from the build folder.

Execution parameters:

- -si followed by the relative path of the file containing the single_end sequence.
```sh    
$ ./MISSH -si ../TestInputFile/reads_800.fa
```
    
- -pi followed by the two relative paths of the two files containing the paired_end sequence.
 
```sh   
$ ./MISSH -pi ../TestInputFile/paired.fna.1 ../TestInputFile/paired.fna.2
```
    
- -dirO followed by the relative path where the program will save the processing times. Default: ../output/
    
```sh   
$ ./MISSH -si ../TestInputFile/reads_800.fa -dirO ../output/test1/
```
    
- -q followed by the relative path of the file containing the spaced seeds which will be used.
  If not used default spaced seeds are:
  - 1111011101110010111001011011111 -> CLARK-S paper
  - 1111101011100101101110011011111 -> CLARK-S paper
  - 1111101001110101101100111011111 -> CLARK-S paper
  - 1111010111010011001110111110111 -> rasbhari minimizing overlap complexity
  - 1110111011101111010010110011111 -> rasbhari minimizing overlap complexity
  - 1111101001011100111110101101111 -> rasbhari minimizing overlap complexity
  - 1111011110011010111110101011011 -> rasbhari maximizing sensitivity
  - 1110101011101100110100111111111 -> rasbhari maximizing sensitivity
  - 1111110101101011100111011001111 -> rasbhari maximizing sensitivity 
    
```sh   
$ ./MISSH -si ../TestInputFile/reads_800.fa -q ../Seeds/Seed_test.fna
```
- -test followed by single or multi specify the kind of test that we want to perform, by default both tests are performed.
    
```sh   
$ ./MISSH -si ../TestInputFile/reads_800.fa -q ../Seeds/Seed_test.fna -test multi
```
- -threads followed by a number set the number of thread for the openmp library (in order to use multiple threads with the "naive" method decommenting a line of code is needed) it is used in the methods that are implemented for supporting the use of multiple cores (MISSHmulticolumnParallel) 

```sh   
$ ./MISSH -si ../TestInputFile/reads_800.fa -q ../Seeds/Seed_test.fna -test multi -threads 4
```

All the examples use test files which are present in this program.

Please note that the hashing results are not saved, but are discarded once computed in order to require only a small amount of memory. This is reasonable because we are interested in the time required for the computation and not in obtaining the resulting hashing.

Citation
---------
Efficient Hashing of Multiple Spaced Seeds with Application
Eleonora Miani, Enrico Petrucci, Cinzia Pizzi and Matteo Comin
Accepted at BIOINFORMATICS 2023 - 14th International Conference on Bioinformatics Models, Methods and Algorithms



