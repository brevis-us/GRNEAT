package fun.grn.grneat.cli;

import fun.grn.grneat.evolver.GRNGene;
import fun.grn.grneat.evolver.GRNGenome;
import fun.grn.grneat.grn.GRNProtein;
import picocli.CommandLine;

import java.util.Random;
import java.util.concurrent.Callable;

@CommandLine.Command(description = "Generate a random GRN genome.",
         name = "grneat-random", mixinStandardHelpOptions = true, version = "grneat 0.0.5")
class RandomGenome implements Callable<Integer> {

    @CommandLine.Option(names = "--filePath", description = "The location to write the GRN.")
    private String filepath;

    @CommandLine.Option(names = "--nInput", description = "Number of GRN inputs.")
    private int nInput;

    @CommandLine.Option(names = "--nOutput", description = "Number of GRN outputs.")
    private int nOutput;

    @CommandLine.Option(names = "--nRegulatory", description = "Number of GRN regulatory proteins.")
    private int nRegulatory;

    @CommandLine.Option(names = "--randomSeed", description = "Random seed to use for generating this genome.")
    private long randomSeed;

    public static void main(String... args) throws Exception {
        //args = new String[]{"--randomSeed", "17", "--filePath", "test.grn", "--nInput", "3", "--nOutput", "3", "--nRegulatory", "10"};

        int exitCode = new CommandLine(new RandomGenome()).execute(args);
        System.exit(exitCode);
    }


    public Integer call() throws Exception {
        Random rng = new Random(randomSeed);

        GRNGenome g = new GRNGenome();
        for (int i=0; i<nInput; i++) {
            GRNGene gene=GRNGene.generateRandomGene(GRNProtein.INPUT_PROTEIN, i, rng);
//					System.err.print(gene+" ");
            g.addGene(gene);
        }
        for (int i=0; i<nOutput; i++) {
            g.addGene(GRNGene.generateRandomGene(GRNProtein.OUTPUT_PROTEIN, i, rng));
        }
        for (int i=0; i<nRegulatory; i++) {
            g.addGene(GRNGene.generateRandomRegulatoryGene(rng));
        }

        g.setBeta(rng.nextDouble()*(g.getBetaMax()-g.getBetaMin())+g.getBetaMin());
        g.setDelta(rng.nextDouble()*(g.getDeltaMax()-g.getDeltaMin())+g.getDeltaMin());

        g.writeToFile(filepath);

        System.out.println(filepath);
        return 0;
    }
}
