# ems2
Improved algorithms for findinding edit distance based motifs

```
Usage: ./ems [OPTIONS] <input-sequence-file>
	-s <version>  Possible versions: 1, 2, 2m, 2p
	-l <l>        Length (l) of (l,d) motif
	-d <d>        Maximum edit distance (d) of (l,d) motif
	-t <int>      Number of threads
```

The option `-t` is applicable to version `2p` only. Other versions are single threaded. 

## Build

```
$ cd src
$ make
```


## Run

```
$ cd src
$ ems -s 2p -l 16 -d 3 -t 8 ../test/planted_l16_d3.txt
```

## References

- Soumitra Pal, Sanguthevar Rajasekaran, "Improved algorithms for finding edit distance based motifs", IEEE International Conference on Bioinformatics and Biomedicine (BIBM) 2015, pp. 537-542, doi:10.1109/BIBM.2015.7359740 
