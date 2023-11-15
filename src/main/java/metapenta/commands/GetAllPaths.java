package metapenta.commands;
import java.util.ArrayList;
import java.util.List;
import metapenta.model.MetabolicNetwork;
import metapenta.model.MetabolicNetworkXMLLoader;

/**
 * Class to test get all paths of metabolic networks
 * @author Valerie Parra Cortés
 */
public class GetAllPaths {
	/**
	 * The main method of class
	 * @param args[0] the path of the XML file 
	 * @param args[1] The initial metabolites separated by comma
	 * @param args[2] Target metabolite
	 * @param args[3] fileName1 Metabolic Network In CSV	
	 * @throws Exception if exists any error of I/O
	 */
	public static void main(String[] args) throws Exception {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		String[] initialMetabolites = args[1].split(",");
		List<String> im= new ArrayList<String>();
		network.makeNet();
		for (int i = 0; i < initialMetabolites.length; i++) {
			im.add(initialMetabolites[i]);			
		}
		network.getAllPaths(im, args[2], args[3]);
	}
}
