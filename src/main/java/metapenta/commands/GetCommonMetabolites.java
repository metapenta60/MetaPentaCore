package metapenta.commands;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import metapenta.model.MetabolicNetwork;
import metapenta.model.Metabolite;
import metapenta.tools.io.MetabolicNetworkXMLFileLoader;

public class GetCommonMetabolites {
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLFileLoader loader = new MetabolicNetworkXMLFileLoader();
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
