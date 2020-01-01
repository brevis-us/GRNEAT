package fun.grn.grneat.cli;

import fun.grn.grneat.evolver.GRNGenome;
import fun.grn.grneat.operators.GRNAddGeneMutationOperator;
import fun.grn.grneat.operators.GRNDeleteGeneMutationOperator;
import fun.grn.grneat.operators.GRNGeneMutationOperator;
import fun.grn.grneat.operators.GRNMutationOperator;
import picocli.CommandLine;

import java.util.Random;
import java.util.concurrent.Callable;

@CommandLine.Command(description = "Measure the distance between two genomes.",
         name = "grneat-distance", mixinStandardHelpOptions = true, version = "grneat 0.0.5")
public class DistanceGenomes implements Callable<Integer> {

    @CommandLine.Option(names = "--pPath", description = "The location of the P genome GRN.")
    private String pPath;

    @CommandLine.Option(names = "--qPath", description = "The location of the Q genome GRN.")
    private String qPath;

    @CommandLine.Option(names = "--compareDynamicsCoeff", description = "Boolean: use the beta/delta coefficients in distance measure.")
    private boolean compareDynamicsCoeff;

    public static void main(String... args) throws Exception {
        //args = new String[]{"--randomSeed", "17", "--filePath", "test.grn", "--nInput", "3", "--nOutput", "3", "--nRegulatory", "10"};

        int exitCode = new CommandLine(new DistanceGenomes()).execute(args);
        System.exit(exitCode);
    }


    public Integer call() throws Exception {
        GRNGenome p = GRNGenome.loadFromFile(pPath);
        GRNGenome q = GRNGenome.loadFromFile(qPath);

        double distance = p.distanceTo(q, compareDynamicsCoeff);

        System.out.println(distance);
        return 0;
    }
}
