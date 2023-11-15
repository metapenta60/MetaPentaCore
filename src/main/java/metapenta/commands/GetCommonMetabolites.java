package metapenta.commands;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import metapenta.model.MetabolicNetwork;
import metapenta.model.MetabolicNetworkXMLLoader;
import metapenta.model.Metabolite;

/**
 * Returns the common metabolites of two differente metabolic networks
 * @param args[0] path of XML file metbolic network 1
 * @param args[0] path of XML file metbolic network 2
 * @author vparrac
 */
public class GetCommonMetabolites {
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network1 = loader.loadNetwork(args[0]);
		MetabolicNetwork network2 = loader.loadNetwork(args[1]);
		
		List<Metabolite> commonMetabolites = network1.commonMetabolites(network2);
		for (Metabolite metabolite : commonMetabolites) {
			System.out.println(metabolite);
		}
		
		StringBuilder commonMetabolitesJSON = new StringBuilder("{\"commonReactions\":[");
		for (int i = 0; i <commonMetabolites.size(); i++) {
			commonMetabolitesJSON.append(commonMetabolites.get(i).toString());
			commonMetabolitesJSON.append((i==commonMetabolites.size()-1)?"":",\n");
		}		
		commonMetabolitesJSON.append("]}");
		Files.write(Paths.get(args[2]), commonMetabolitesJSON.toString().getBytes());
	}
}
