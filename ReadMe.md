# MISSH

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
$ ./ISSH -si ../TestInputFile/reads_800.fa
```
    
- -pi followed by the two relative paths of the two files containing the paired_end sequence.
 
```sh   
$ ./ISSH -pi ../TestInputFile/paired.fna.1 ../TestInputFile/paired.fna.2
```
    
- -dirO followed by the relative path where the program will save the processing times. Default: ../output/
    
```sh   
$ ./ISSH -si ../TestInputFile/reads_800.fa -dirO ../output/test1/
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
$ ./ISSH -si ../TestInputFile/reads_800.fa -q ../Seeds/Seed_test.fna
```
- -test followed by single or multi specify the kind of test that we want to perform, by default both tests are performed.
    
```sh   
$ ./ISSH -si ../TestInputFile/reads_800.fa -q ../Seeds/Seed_test.fna -test multi
```
- -threads followed by a number set the number of thread for the openmp library (in order to use multiple threads with the "naive" method decommenting a line of code is needed) it is used in the methods that are implemented for supporting the use of multiple cores (ISSHmulticolumnParallel) 

```sh   
$ ./ISSH -si ../TestInputFile/reads_800.fa -q ../Seeds/Seed_test.fna -test multi -threads 4
```

All the examples use test files which are present in this program.

Please note that the hashing results are not saved, but are discarded once computed in order to require only a small amount of memory. This is reasonable because we are interested in the time required for the computation and not in obtaining the resulting hashing.

Citation
---------
Enrico Petrucci, Laurent Noé, Cinzia Pizzi, and Matteo Comin,
"Iterative Spaced Seed Hashing: Closing the Gap Between Spaced Seed Hashing and k-mer Hashing"
Journal of Computational Biology 2020 27:2, 223-233 
http://doi.org/10.1089/cmb.2019.0298

etrucci E., Noé L., Pizzi C., Comin M. (2019) "Iterative Spaced Seed Hashing: Closing the Gap Between Spaced Seed Hashing and k-mer Hashing"
In Bioinformatics Research and Applications. ISBRA 2019. Lecture Notes in Computer Science, vol 11490, pp 208-219.
https://doi.org/10.1007/978-3-030-20242-2_18

