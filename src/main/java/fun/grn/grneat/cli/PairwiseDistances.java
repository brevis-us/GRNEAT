package fun.grn.grneat.cli;

import fun.grn.grneat.evolver.GRNGenome;
import picocli.CommandLine;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

@CommandLine.Command(description = "Measure all v all pairwise distances between a directory of GRNs.",
         name = "grneat-pairwise-distances", mixinStandardHelpOptions = true, version = "grneat 0.0.5")
class PairwiseDistances implements Callable<Integer> {

    @CommandLine.Option(names = "--directory", description = "Directory of GRNs.")
    private String directory;

    @CommandLine.Option(names = "--compareDynamicsCoeff", description = "Boolean: use the beta/delta coefficients in distance measure.")
    private boolean compareDynamicsCoeff;

    public static void main(String... args) throws Exception {
        //args = new String[]{"--randomSeed", "17", "--filePath", "test.grn", "--nInput", "3", "--nOutput", "3", "--nRegulatory", "10"};

        int exitCode = new CommandLine(new PairwiseDistances()).execute(args);
        System.exit(exitCode);
    }


    public Integer call() throws Exception {
        List<GRNGenome> genomes = new ArrayList<GRNGenome>();
        List<String> filenames = new ArrayList<String>();
        File listing = new File(directory);
        for( File f : listing.listFiles() ) {
            if( f.getName().endsWith(".grn") ) {

                GRNGenome g = GRNGenome.loadFromFile(f.getAbsolutePath());
                filenames.add( f.getName() );
                genomes.add(g);
            }
        }

        double[][] distance = new double[genomes.size()][genomes.size()];

        System.out.println(genomes.size());

        for( int k = 0; k < genomes.size(); k++ ) {
            System.out.print("\"" + filenames.get(k) + "\"");
            if( k < genomes.size()-1 ) System.out.print("\t");
        }
        System.out.println();

        for( int p = 0; p < genomes.size(); p++ ) {
            for( int q = 0; q < genomes.size(); q++ ) {
                double d = genomes.get(p).distanceTo(genomes.get(q), compareDynamicsCoeff);
                distance[p][q] = d;
                System.out.print(d);
                if( q < genomes.size()-1 )
                    System.out.print("\t");
            }
            System.out.println();
        }

        return 0;
    }
}
