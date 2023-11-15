package metapenta.commands;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import metapenta.model.MetabolicNetwork;
import metapenta.model.MetabolicNetworkXMLLoader;
import metapenta.model.Reaction;

/**
 * Returns the common reaction of two different metabolic networks
 * @param args[0] Path of metabolic network 1
 * @param args[0] Path of metabolic network 2
 * @author Valerie Parra Cortés
 */
public class GetCommonReactions {
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network1 = loader.loadNetwork(args[0]);
		MetabolicNetwork network2 = loader.loadNetwork(args[1]);
		List<Reaction> commonReations = network1.commonReactions(network2);
		StringBuilder reactions = new StringBuilder("{\"commonReactions\":[");
		for (int i = 0; i <commonReations.size(); i++) {
			reactions.append(commonReations.get(i).toString());
			reactions.append((i==commonReations.size()-1)?"":",\n");
		}		
		
		reactions.append("]}");
		Files.write(Paths.get(args[2]), reactions.toString().getBytes());
	}

}
