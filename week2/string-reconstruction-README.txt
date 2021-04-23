just stitched three previous bits of code together

1. create debruijn path from kmers input_data.txt > a1.out
2. create eulerian path from graph a1.out > a2.out
3. stitch dna string from eulerian path a2.out


e.g.
./string-reconstruction-build-fragment-graph.py ~/Downloads/dataset_203_7.txt > a1.out
./string-reconstruction-create-eulerian-path.py a1.out > a2.out
./string-reconstruction-stitch-dna-from-eulerian-path-kmers.py a2.out > a3.out
