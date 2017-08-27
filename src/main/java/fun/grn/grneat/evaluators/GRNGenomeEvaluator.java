package fun.grn.grneat.evaluators;

import java.util.Vector;

import fun.grn.grneat.evolver.GRNGene;
import fun.grn.grneat.evolver.GRNGenome;
import fun.grn.grneat.grn.GRNModel;
import fun.grn.grneat.grn.GRNProtein;

public class GRNGenomeEvaluator {
	public int generation=0;// some fun.grn.grneat.evaluators are dynamic problems dependent on generation number
	public boolean nonCacheable = false;// some fun.grn.grneat.evaluators cannot reuse fitness (i.e. dynamic problems)

	public int numGRNInputs=0;
	public int numGRNOutputs=0;
	public static int numEvaluations=0;

	public String name="SuperClass";
	
	public double evaluate(GRNGenome aGenome) {
		return 0;
	}
	
	public static GRNModel buildGRNFromGenome(GRNGenome aGenome) {
		Vector<GRNProtein> prots=new Vector<GRNProtein>();
		for (GRNGene gi : aGenome.getInputGenes()) {
			prots.add(gi.getProtein());
		}
		for (GRNGene go : aGenome.getOutputGenes()) {
			prots.add(go.getProtein());
		}
		for (GRNGene gr : aGenome.getRegulatoryGenes()) {
			prots.add(gr.getProtein());
		}
		GRNModel p=new GRNModel(prots, aGenome.getBeta(), aGenome.getDelta());
		return p;
	}

}
