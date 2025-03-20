package cli;

import workflows.negativeWGMs.NegativeWGM_H0Full;

import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.PrintStream;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

@Command(name = "negSimFull", description = "Simulated LRT with H0Full on negative WGMs", mixinStandardHelpOptions = true)
public class Neg_SimH0Full_cmd implements Runnable{
    @Option(names = {"-t", "--tree"}, description = "Species tree filepath", paramLabel = "SPECIES_TREE", required = true)
    private String treeFile;

    @Option(names = {"-w", "--wgm"}, description = "WGM list filepath", paramLabel = "WGM_LIST", required = true)
    private String wgdFile;

    @Option(names = {"-gf", "--gene_family_rank"}, description = "Gene family rank (starts from 0)", paramLabel = "CURRENT_GF_RANK", required = true)
    private int gfNumber;

    @Option(names = {"-r", "--ranking"}, description = "Gene family ranking filepath", paramLabel = "LAMBDA_RANKING", required = true)
    private String combinedOutputFile;

    @Option(names = {"-obs", "--observed_lr"}, description = "Negative observed LR current gene family filepath", paramLabel = "OBSERVED_LRT", required = true)
    private String lrtOriginalFile;

    @Option(names = {"-out", "--output_file"}, description = "Output filepath", paramLabel = "OUTPUT_FILE", required = true)
    private String outputFile;


    @Override
    public void run(){

        try (PrintStream output = new PrintStream(new FileOutputStream(outputFile))) {
            // Redirect stdout to the output file
            System.setOut(output);
            
            // Create and execute the NegativeWGM_H0Full instance
            NegativeWGM_H0Full lrtDist = new NegativeWGM_H0Full(treeFile, wgdFile, gfNumber, combinedOutputFile, lrtOriginalFile);
            lrtDist.execute();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            // Restore stdout
            System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
        }
    }
}
