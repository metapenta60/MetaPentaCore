package metapenta.commands;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;

import metapenta.model.MetabolicNetwork;
import metapenta.model.MetabolicNetworkXMLLoader;
/**
 * Class to test connected components method
 * @author Valerie Parra Cort�s
 */
public class ConnectedComponents {

	/**
	 * Class to test the conncected components of the metabolic networkd
	 * @param args [0] the path of the XML file and args[1] out file in csv
	 * @throws Exception if exists any error of I/O
	 */
	public static void main(String[] args) throws Exception {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		network.makeNet();
		Map<String, Integer> connectedComponents = network.connectedComponents();
		StringBuilder csv = new StringBuilder();
		Set<String> reactionsSet = connectedComponents.keySet();
		for (String reaction : reactionsSet) {
			csv.append(reaction+","+connectedComponents.get(reaction)+"\n");
		}
		Files.write(Paths.get(args[1]), csv.toString().getBytes());
	}
}
