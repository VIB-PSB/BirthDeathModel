package cli;
import picocli.CommandLine;
import picocli.CommandLine.Command;

@Command(name = "cLRTs_for_WGMs",
        description = "Complementary LRTs for WGM inference",
        subcommands = {Pos_ObsLR_cmd.class, Pos_SimH0Full_cmd.class, Pos_SimH0Rm_cmd.class, Pos_PostProc_cmd.class,  Neg_ObsLR_cmd.class, Neg_SimH0Full_cmd.class, Neg_SimH0Neg_cmd.class, Neg_PostProc_cmd.class},
        mixinStandardHelpOptions = true,
        version = "1.0")
public class cLRTs_for_WGMs_CLI implements Runnable {

    public static void main(String[] args) {
        int exitCode = new CommandLine(new cLRTs_for_WGMs_CLI()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public void run() {
        System.out.println("Welcome to cLRTs_for_WGMs command line tool! Please run '-h' for usage.");
    }
}
