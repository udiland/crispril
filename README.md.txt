
# CRISPys: target multiple genes CRISPR sites
CRISPys is used to find a set of CRISPR sgRNA guides that will capture multiple genes in a group of user defined genes family.

## Features
CRISPys can be used to find set of sgRNAs that will target as many genes in the family. 
It can also be used to consider the homology of the genes and find the best guides for each homology group, It is done by constracting a gene tree and by finding sgRNA guides for each gene subgroup (sometimes refferd to as 'internal node') 
This feature can be controlled by choosing "alg A" to search in all genes, or "alg E" to consider homology.

To restrict the search for targets to capture in only part of the gene use "where_in_gene" argument (integer 0-1).
The scoring function that ranks the cleavage efficiency of sgRNA is controlled with the "df_targets" argument ("cfd_funct" will use the CFD score).
The threshold of the scoring value to consider a target to cleave a site can be control with the "Omege" argument (integer 0-1). 
If choosing the "alg E" option, the number of guides to consider per internal node can be control with "internal_node_candidates" (integer).
The number of polymorphic sites (the number of difrent bases between the guides) to consider when designing sgRNA can be controlled with "PS_number" (integer).


## Dependencies
See `app/requirements.txt`
Also, some programs are installed while building the docker container.

## Run using docker 
1. Build docker image:
~~~
docker build -t crispys ./ -f Dockerfile
docker run -p8000:80 crispys
~~~
2. Upload a file of genes you want to target in Fasta format (each Fasta header is a gene name, exons should follow one another with the same header).
~~~
curl -F "file=@/path/to/fasta/file" localhost:8000/upload
~~~
4. Run the program with the described parameters as shown in the following example:  
~~~
curl --location --request POST 'localhost:8000/crispys' \
--header 'Content-Type: application/json' --data-raw '{
    "fasta_file": "name_of_fasta_file",	
    "alg": "E",
    "where_in_gene": 0.6,
    "Omega": 0.8,
    "df_targets": "cfd_funct",
    "internal_node_candidates": 200,
    "PS_number":12
}'
~~~
## Get the results (they will be saved as 'output.csv')
~~~
curl -o output.csv localhost:8000/download
~~~

## Links
CRISPys paper:
[Hyams *et. al*. *L Mol Biol* (2018)](https://www.sciencedirect.com/science/article/abs/pii/S0022283618301682?via%3Dihub)
CRISPys online:
http://multicrispr.tau.ac.il/
