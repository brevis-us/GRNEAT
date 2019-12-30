package fun.grn.grneat.cli;

import fun.grn.grneat.evolver.GRNGene;
import fun.grn.grneat.evolver.GRNGenome;
import fun.grn.grneat.grn.GRNProtein;
import fun.grn.grneat.operators.GRNAddGeneMutationOperator;
import fun.grn.grneat.operators.GRNDeleteGeneMutationOperator;
import fun.grn.grneat.operators.GRNGeneMutationOperator;
import fun.grn.grneat.operators.GRNMutationOperator;
import picocli.CommandLine;

import java.util.Random;
import java.util.concurrent.Callable;

@CommandLine.Command(description = "Mutate a GRN genome.",
         name = "grneat-mutate", mixinStandardHelpOptions = true, version = "grneat 0.0.5")
class MutateGenome implements Callable<Integer> {

    @CommandLine.Option(names = "--inputPath", description = "The location of source GRN.")
    private String inputPath;

    @CommandLine.Option(names = "--outputPath", description = "The location of output GRN.")
    private String outputPath;

    // Add gene operation
    @CommandLine.Option(names = "--pAddMutation", description = "Probability of adding genes.")
    private double pAddMutation;

    @CommandLine.Option(names = "--addMutationMaxSize", description = "Max size of add mutation.")
    private int addMutationMaxSize;

    // Delete gene operation
    @CommandLine.Option(names = "--pDelMutation", description = "Probability of deleting genes.")
    private double pDelMutation;

    @CommandLine.Option(names = "--delMutationMinSize", description = "Minimum size to delete.")
    private int delMutationMinSize;

    // Mutate gene operation
    @CommandLine.Option(names = "--pChangeMutation", description = "Probability of changing genes.")
    private double pChangeMutation;

    @CommandLine.Option(names = "--randomSeed", description = "Random seed to use for generating this genome.")
    private long randomSeed;

    public static void main(String... args) throws Exception {
        //args = new String[]{"--randomSeed", "17", "--filePath", "test.grn", "--nInput", "3", "--nOutput", "3", "--nRegulatory", "10"};

        int exitCode = new CommandLine(new MutateGenome()).execute(args);
        System.exit(exitCode);
    }


    public Integer call() throws Exception {
        Random rng = new Random(randomSeed);

        GRNGenome g = GRNGenome.loadFromFile(inputPath);

        GRNMutationOperator operator;
        GRNGenome child = null;
        while( child == null ) {
            double r = rng.nextDouble();
            if( r < pAddMutation ) {
                operator = new GRNAddGeneMutationOperator(addMutationMaxSize, pAddMutation);
            } else if( r < pAddMutation + pDelMutation ) {
                operator = new GRNDeleteGeneMutationOperator(Math.max(delMutationMinSize, g.getInputGenes().size()+g.getOutputGenes().size()+1), pDelMutation);
            } else if( r < pAddMutation + pDelMutation + pChangeMutation ) {
                operator = new GRNGeneMutationOperator(pChangeMutation);
            } else {
                throw new Exception("No mutation applied");
            }
            child = operator.cloneAndMutate(g, rng);
        }

        child.writeToFile(outputPath);

        System.out.println(outputPath);
        return 0;
    }
}
