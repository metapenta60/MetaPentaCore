package metapenta.commands;

import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import metapenta.model.MetabolicNetwork;
import metapenta.model.Reaction;
import metapenta.tools.io.MetabolicNetworkXMLFileLoader;

/**
 * Method that gived a Enzyme ID print a JSON with the information of reactions
 * where Enzyme is required
 * File example output
 * {
  "reactions": [
    {
      "id": "R_ME2",
      "name": "Malic enzyme (NADP)",
      "reversible": false,
      "reactants": [
        {
          "metaboliteId": "M_mal__L_c",
          "metaboliteName": "L-Malate",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_nadp_c",
          "metaboliteName": "Nicotinamide adenine dinucleotide phosphate",
          "stoichiometry": "1.0"
        }
      ],
      "products": [
        {
          "metaboliteId": "M_co2_c",
          "metaboliteName": "CO2 CO2",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_nadph_c",
          "metaboliteName": "Nicotinamide adenine dinucleotide phosphate - reduced",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_pyr_c",
          "metaboliteName": "Pyruvate",
          "stoichiometry": "1.0"
        }
      ]
    }
  ]
}

 * @author vparrac
 *
 */
public class GetReactionsCatalyzedBy {
	/**
	 * The main method of class
	 * @param args[0] the path of the XML file of metabolic network
	 * @param args[1] The enzyme to find the reactions
	 * @param args[2] Name of file out
	 * @throws Exception if exists any error of I/O
	 */
	public static void main(String[] args) throws IOException  {
		MetabolicNetworkXMLFileLoader loader = new MetabolicNetworkXMLFileLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		List<Reaction> result = network.getReactionsCatalyzedBy(args[1]);
		try (PrintStream out = new PrintStream(args[2])) {
			out.print("{reactions:[");
			for (int i = 0; i < result.size(); i++) {
				if(i==result.size()-1) {
					out.print(result.get(i).toString());
				}else {
					out.print(result.get(i).toString()+",");
				}
			}
			out.print("]}");
		}
	}
}


