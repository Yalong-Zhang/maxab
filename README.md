# Integral Densest Subgraph Search on Directed Graphs

# Datasets

Because some datasets used in the paper are too large to be uploaded to GitHub, we have summarized the download links for the dataset in the table below.

Datasets used for performance studies.

| Dataset | Link |
| --- | --- |
| IM | https://networkrepository.com/IMDB.php |
| HE | https://networkrepository.com/cit-HepPh.php |
| WI | https://networkrepository.com/web-wikipedia-link-fr.php |
| IN | https://networkrepository.com/web-indochina-2004-all.php |
| OR | https://networkrepository.com/aff-orkut-user2groups.php |
| UK | https://networkrepository.com/web-uk-2005-all.php |
| IT | https://networkrepository.com/web-it-2004-all.php |
| TW | https://networkrepository.com/soc-twitter-mpi-sws.php |
| Epinions | https://snap.stanford.edu/data/soc-Epinions1.html |

# Preprocess

The dataset file needs to be preprocessed as the format below, where for every directed edge (u_i, v_i), u_i is the staring node and v_i is the ending node.

```
<|V|> <|E|>
<u_1> <v_1>
<u_2> <v_2>
<u_3> <v_3>
...
```

# Usage of Exact Algorithms

Create a folder named "output" to store the output subgraph.

```
g++ main.cpp -o main -std=c++11 -O3
./main <dataset_address> <argument>
```

When argument is -a, algorithm GetIDS is used. When argument is -d, algorithm GetIDS++ is used. 

For example, to use the algorithm GetIDS++ on the dataset example.txt, you should use:

```
g++ main.cpp -o main -std=c++11 -O3
./main ./example.txt -d
```

# Usage of Approximation Algorithms

Create a folder named "output" to store the output subgraph.

```
g++ appro.cpp -o appro -std=c++11 -O3
./appro <dataset_address> <argument1> <argument2>
```

When argument1 is -s, algorithm SingleCore is used, and argument2 is blank. When argument is -m, algorithm MultiCore is used, and argument2 is epsilon. When argument1 is -a, algorithm AllCore is used, and argument2 is blank.

For example, to use the algorithm MultiCore on the dataset example.txt with an epsilon of 0.5, you should use:

```
g++ appro.cpp -o appro -std=c++11 -O3
./appro ./example.txt -m 0.5
```