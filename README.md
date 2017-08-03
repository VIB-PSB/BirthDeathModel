# BirthDeathModel

## About

Given as input a species tree (with additional WGM nodes) and a gene family size profile across the species concerned, the model computes the likelihood of the observed gene family size profile conditioned on a duplication/loss rate λ and a gene family size r at the root of the species tree. 
For any given gene family, this likelihood was maximized as a function of λ for any fixed r between 1 and 20 using a cutting-plane optimization method

## Installation

Download the code as a zip file or `git clone` the repository from [GitHub][].

[GitHub]: https://github.com/VIB-PSB/BirthDeathModel

You'll need [Maven][] to build the executable jar file:

```
mvn clean install
```

If the build is successful, a jar file with embedded dependencies can now be found in the `target` folder.


[Maven]: http://maven.apache.org/

## Usage

### Input ##
The main program takes 5 parameters:

1. A text file with the Newick format of the tree in it,
ending in the root node , followed by ";".

2. A text file representing WGM events. Each WGM must be
present in one line. In case successive WGMs happened on one
branch the events should be present from the oldest to the
youngest. For example WGD from node A to B with branch length x is
presented as: WGD, A, B, x where x is the distance of the WGD node
to the parent node.

3. Gene family counts. this file must contain at every line,
a gene family ID, followed by gene counts in different species in
the same order that the species appear in the Newick format tree.

4. An integer representing the number of gene family in the
gene counts file, starting from 0.

5. An integer specifying the gene counts at the root node (root size).

### Output ##
The program returns the following line for every run:

```
Gene_Family_ID \t root_size \t optimal_lambda \t log_likelihood \t p-value 
```