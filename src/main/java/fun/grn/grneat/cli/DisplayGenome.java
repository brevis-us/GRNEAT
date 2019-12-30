package fun.grn.grneat.cli;

import fun.grn.grneat.evolver.GRNGene;
import fun.grn.grneat.evolver.GRNGenome;
import fun.grn.grneat.grn.GRNProtein;
import picocli.CommandLine;

import java.util.Random;
import java.util.concurrent.Callable;

@CommandLine.Command(description = "Generate a random GRN genome.",
         name = "grneat-display", mixinStandardHelpOptions = true, version = "grneat 0.0.5")
class DisplayGenome implements Callable<Integer> {

    @CommandLine.Option(names = "--filePath", description = "The location to write the GRN.")
    private String filepath;

    public static void main(String... args) throws Exception {
        //args = new String[]{"--filePath", "test.grn"};

        int exitCode = new CommandLine(new DisplayGenome()).execute(args);
        System.exit(exitCode);
    }


    public Integer call() throws Exception {
        GRNGenome g = GRNGenome.loadFromFile(filepath);

        System.out.println("input_proteins:");
        for( GRNGene el : g.getInputGenes() ) {
            System.out.println("\t" + el.toString() );
        }

        System.out.println("output_proteins:");
        for( GRNGene el : g.getOutputGenes() ) {
            System.out.println("\t" + el.toString() );
        }

        System.out.println("regulatory_proteins:");
        for( GRNGene el : g.getRegulatoryGenes() ) {
            System.out.println("\t" + el.toString() );
        }

        System.out.println("beta: " + g.getBeta());
        System.out.println("delta: " + g.getDelta());

        return 0;
    }
}
