package cli;

import workflows.negativeWGMs.NegativeWGM_PostprocessingGFListFDR;
import workflows.negativeWGMs.NegativeWGM_PostprocessingGFRangeFDR;

import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.PrintStream;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command(name = "negPP", description = "Post-Processing of negative WGMs", mixinStandardHelpOptions = true)
public class Neg_PostProc_cmd implements Runnable {

    @Option(names = {"-sim_full", "--neg_sim_h0full_lrt"}, description = "Directory of simulated LRTs under H0=Full", paramLabel = "SIMULATED_H0FULL", required = true)
    private String negatives_H0Full_outputDir;

    @Option(names = {"-sim_neg", "--neg_sim_h0neg_lrt"}, description = "Directory of simulated LRTs under H0=Neg", paramLabel = "SIMULATED_H0NEG", required = true)
    private String negatives_H0Rm_outputDir;

    @Option(names = {"-obs", "--neg_obs_lrt"}, description = "Directory of observed LRTs", paramLabel = "OBSERVED_LRT", required = true)
    private String origLRTDir;

    @Option(names = {"-num_w", "--num_wgms"}, description = "Number of tested WGMs", paramLabel = "NUMBER_WGM", required = true)
    private int wgdnumber;

    @Option(names = {"-fdr", "--fdr_threshold"}, description = "FDR threshold", paramLabel = "FDR_THRESHOLD", required = true)
    private double fdr;

    @Option(names = {"-gf_list", "--gene_family_rank_list"}, description = "Comma-separated list of GF ranks", paramLabel = "GF_RANK_LIST")
    private String listGFs;

    @Option(names = {"-gf_range", "--gene_family_rank_range"}, arity = "4", description = "Space-separated start and end ranks of top and bottom GFs", paramLabel = "GF_RANK_RANGE")
    private String[] rangeGFs;

    @Option(names = {"-out", "--output_file"}, description = "Output filepath", paramLabel = "OUTPUT_FILE", required = true)
    private String outputFile;


    @Override
    public void run(){

        // If both list and range are provided, print an error message and exit
        if (listGFs != null && rangeGFs != null) {
            System.err.println("Error: can't provide both a GF list and a GF range at the same time.");
            System.exit(1);
        }

        // If neither list or range are provided, print an error message and exit
        if (listGFs == null && rangeGFs == null) {
            System.err.println("Error: must provide either a GF list or a GF range.");
            System.exit(1);
        }

        // Either use LIST of GF ranks...
        if (listGFs != null) {

            try (PrintStream output = new PrintStream(new FileOutputStream(outputFile))) {
                // Redirect stdout to the output file
                System.setOut(output);
                
                // Create and execute the NegativeWGM_PostprocessingGFListFDR instance
                NegativeWGM_PostprocessingGFListFDR lrtDist = new NegativeWGM_PostprocessingGFListFDR(negatives_H0Full_outputDir, negatives_H0Rm_outputDir, origLRTDir, listGFs, wgdnumber, fdr);
                lrtDist.execute();
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                // Restore stdout
                System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
            }

        // ... or use RANGE of GF ranks
        } else if (rangeGFs != null) {

            String startTopGF = rangeGFs[0];
            String endTopGF = rangeGFs[1];
            String startBottomGF = rangeGFs[2];
            String endBottomGF = rangeGFs[3];

            try (PrintStream output = new PrintStream(new FileOutputStream(outputFile))) {
                // Redirect stdout to the output file
                System.setOut(output);
                
                // Create and execute the NegativeWGM_PostprocessingGFRangeFDR instance
                NegativeWGM_PostprocessingGFRangeFDR lrtDist = new NegativeWGM_PostprocessingGFRangeFDR(negatives_H0Full_outputDir, negatives_H0Rm_outputDir, origLRTDir, startTopGF, endTopGF, startBottomGF, endBottomGF, wgdnumber, fdr);
                lrtDist.execute();
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                // Restore stdout
                System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
            }
        }

        // Exit if neither LIST nor RANGE provided
        else {
            System.out.println("Error: Please either provide a comma-separated list of GF ranks or a range of GF ranks.");
        }
    }
}
