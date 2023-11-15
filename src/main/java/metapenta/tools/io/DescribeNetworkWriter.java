package metapenta.tools.io;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class DescribeNetworkWriter {

    private StringBuilder metabolitesBuilder = new StringBuilder();
    private StringBuilder reactionsBuilder = new StringBuilder();
    private StringBuilder sMatrixBuilder = new StringBuilder();
    private StringBuilder eMetabolitesFileBuilder = new StringBuilder();
    private StringBuilder cMetabolitesFileBuilder = new StringBuilder();
    private StringBuilder reversibleRxnBuilder = new StringBuilder();

    private StringBuilder irreversibleRxnBuilder = new StringBuilder();
    private StringBuilder upBoundFileBuilder = new StringBuilder();
    private StringBuilder loBoundsFileBuilder = new StringBuilder();
    private String prefix;
    private String metabolitesFileName;
    private String reactionsFile;
    private String sMatrixFile;
    private String eMetabolitesFile;
    private String cMetabolitesFile;
    private String reversibleRxnFile;
    private String irreversibleRxnFile;
    private String upBoundFile;
    private String loBoundsFile;
    private StringUtils stringUtils;

    public DescribeNetworkWriter(String prefix){
        this.prefix = prefix;
        this.metabolitesFileName = prefix + "_compounds.txt";
        this.reactionsFile = prefix + "_reactions.txt";
        this.sMatrixFile = prefix + "_s_matrix.txt";
        this.eMetabolitesFile = prefix + "_extracellularMetabolites.txt";
        this.cMetabolitesFile = prefix + "_cytosolic_Metabolites.txt";
        this.reversibleRxnFile = prefix + "_reversible_rxn.txt";
        this.upBoundFile = prefix + "_upperbounds_on_fluxes.txt";
        this.loBoundsFile = prefix + "_lowerbounds_on_fluxes.txt";
        this.irreversibleRxnFile = prefix + "_irreversible_reactions";
        this.stringUtils = new StringUtils();
    }

    public void WriteInMetabolites(String metabolite){
        stringUtils
                .SetString(metabolite)
                .addSingleQuotes()
                .addEmptySpace()
                .addBreakLine();

        metabolitesBuilder.append(stringUtils.GetString());
    }

    public void writeMetabolites(List<String> ids){
        for(String metaboliteId: ids){
         WriteInMetabolites(metaboliteId);
        }
    }

    public void writeEMetabolites(List<String> ids){
        for(String metaboliteId: ids){
            WriteEMetabolite(metaboliteId);
        }
    }

    public void writeCMetabolites(List<String> ids){
        for(String metaboliteId: ids){
            WriteCMetabolite(metaboliteId);
        }
    }

    public void writeReactions(List<String> ids){
        for(String reactionId: ids){
            WriteInReactions(reactionId);
        }
    }

    public void writeReversibleReactions(List<String> ids){
        for(String reactionId: ids){
            WriteInReversibleReactions(reactionId);
            writeBoundsForReversibleReaction(reactionId);
        }
    }


    public void writeIrreversibleReactions(List<String> ids){
        for(String reactionId: ids){
            writeIrreversibleReactions(reactionId);
            writeBoundsForReversibleReaction(reactionId);
        }
    }


    public void WriteInReactions(String reactionID){
        stringUtils
                .SetString(reactionID)
                .addSingleQuotes()
                .addBreakLine();

        reactionsBuilder.append(stringUtils.GetString());
    }

    public void WriteInSMatrix(String metaboliteID, double stoichiometry){
        stringUtils.SetString(metaboliteID)
                .addSingleQuotes()
                .addEmptySpace()
                .addDouble(stoichiometry)
                .addBreakLine();

        sMatrixBuilder.append(stringUtils.GetString());
    }

    public void WriteEMetabolite(String metaboliteID){
        stringUtils
                .SetString(metaboliteID)
                .addSingleQuotes()
                .addBreakLine();

        eMetabolitesFileBuilder.append(stringUtils.GetString());
    }

    public void WriteCMetabolite(String metaboliteID){
        stringUtils
                .SetString(metaboliteID)
                .addSingleQuotes()
                .addBreakLine();

        cMetabolitesFileBuilder.append(stringUtils.GetString());
    }

    public void WriteImUpBound(String s){
        upBoundFileBuilder.append(s);
    }

    public void WriteInReversibleReactions(String reactionID){
        stringUtils
                .SetString(reactionID)
                .addSingleQuotes()
                .addBreakLine();

        reversibleRxnBuilder.append(stringUtils.GetString());
    }

    public void writeIrreversibleReactions(String reactionID){
        stringUtils
                .SetString(reactionID)
                .addSingleQuotes()
                .addBreakLine();

        irreversibleRxnBuilder.append(stringUtils.GetString());
    }

    public void writeBoundsForReversibleReaction(String reactionName){
        stringUtils.
                SetString(reactionName).
                addSingleQuotes().
                addEmptySpace().
                addInt(-1000).
                addBreakLine();

        WriteInLowBound(stringUtils.GetString());

        stringUtils.
                SetString(reactionName);
        stringUtils.
                addSingleQuotes().
                addEmptySpace().
                addInt(1000).
                addBreakLine();

        WriteImUpBound(stringUtils.GetString());
    }

    public void writeBoundsForUnreversibleReaction(String reactionId){
        StringUtils stringUtils = new StringUtils();
        stringUtils.
                SetString(reactionId).
                addSingleQuotes().
                addEmptySpace().
                addInt(0).
                addBreakLine();
        WriteInLowBound(stringUtils.GetString());

        stringUtils.
                SetString(reactionId).
                addSingleQuotes().
                addEmptySpace().
                addInt(1000).
                addBreakLine();

        WriteImUpBound(stringUtils.GetString());
    }

    public void WriteInLowBound(String s){
        loBoundsFileBuilder.append(s);
    }

    public void Write() throws Exception{
        Files.write(Paths.get(metabolitesFileName), metabolitesBuilder.toString().getBytes());
        Files.write(Paths.get(reactionsFile), reactionsBuilder.toString().getBytes());
        Files.write(Paths.get(sMatrixFile), sMatrixBuilder.toString().getBytes());
        Files.write(Paths.get(eMetabolitesFile), eMetabolitesFileBuilder.toString().getBytes());
        Files.write(Paths.get(cMetabolitesFile), cMetabolitesFileBuilder.toString().getBytes());
        Files.write(Paths.get(reversibleRxnFile), reversibleRxnBuilder.toString().getBytes());
        Files.write(Paths.get(upBoundFile), upBoundFileBuilder.toString().getBytes());
        Files.write(Paths.get(loBoundsFile), loBoundsFileBuilder.toString().getBytes());
    }
}
