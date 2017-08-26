package fun.grn.grneat.operators;

import java.util.Random;

import fun.grn.grneat.evolver.GRNGenome;


public abstract class GRNCrossoverOperator {
	public String name="SuperCrossover!";
	
	public abstract GRNGenome reproduce(GRNGenome parent1, GRNGenome parent2, Random rng);
}
