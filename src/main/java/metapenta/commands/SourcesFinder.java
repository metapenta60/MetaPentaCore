package metapenta.commands;

import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import metapenta.model.MetabolicNetwork;
import metapenta.model.MetabolicNetworkXMLLoader;
import metapenta.model.Metabolite;
/**
 * Class to test sources finder method 
 * @author Valerie Parra Cortés
 */
public class SourcesFinder {
	/**
	 * The main method of class
	 * @param args[0] the path of the XML file of metabolic network
	 * @param args[1] the path of output file
	 *  @throws Exception if exists any error of I/O
	 */
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		network.makeNet();
		List<Metabolite> sinks = network.findSources();
		try (PrintStream out = new PrintStream(args[1])) {
			out.print("{sources:[");
			for (int i = 0; i < sinks.size(); i++) {				
				String s = (i>0)?",":"";
				out.print(s+sinks.get(i).toString());
			}
			out.print("]}");			
		}		
	}
}
