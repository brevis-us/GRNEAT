
# Gene Regulatory Network Evolution Through Augmenting Topologies

# Quickstart:

This code was originally made in Eclipse, and works well there with the standard "Run" mechanism.

Go to: java/fun.fun.grn.grneat.grn.grneat.fun.grn.grneat.evolver/Evolver.java and press Run

# To change experiments: 

Go to the main() function in Evolver.java

Change line:
			e.evaluator = new IntertwinedSpirals( args );
			
To the evaluator of your choice.

The main function was designed for command line usage in a .jar, but by swapping the evaluator class that you initialize you can easily change the problem within the IDE of your choice.

# Command line usage

Install:

`./install $HOME/bin`

Usage:
```
grneat-random --filePath=test2.grn --nInput=3 --nOutput=3 --nRegulatory=10 --randomSeed=17
grneat-display --filePath=test2.grn
grneat-random --filePath=test2.grn --nInput=3 --nOutput=3 --nRegulatory=10 --randomSeed=17
grneat-mutate --inputPath=test3.grn --outputPath=test3_child.grn --pAddMutation=0.3 --addMutationMaxSize=5 --delMutationMinSize=3 --pDelMutation=0.3 --pChangeMutation=0.3 --randomSeed=12245
grneat-distance --pPath=test3.grn --qPath=test3_child.grn
```

# Leiningen:

```
:repositories [["brevis" "http://dl.bintray.com/kephale/brevis"]]
```

```
:dependencies [[us.brevis/GRNEAT "0.0.1"]]
```

# Citing:

Cussat-Blanc, S., Harrington, K., and Pollack, J. (2015) Gene Regulatory Network Evolution Through Augmenting Topologies. IEEE Transactions on Evolutionary Computation 19(6), pp. 823 - 837.
http://dx.doi.org/10.1109/TEVC.2015.2396199
