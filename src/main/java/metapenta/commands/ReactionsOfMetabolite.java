package metapenta.commands;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

import metapenta.model.MetabolicNetwork;
import metapenta.model.Reaction;
import metapenta.tools.io.MetabolicNetworkXMLFileLoader;

/**
 * Method that gived a metabolite ID print a JSON with the information of reactions
 * where metabolites is Substrate or Product
 * Output file example:
 * {
  "isSubstrate": [
    {
      "id": "R_ACKr",
      "name": "Acetate kinase",
      "reversible": true,
      "reactants": [
        {
          "metaboliteId": "M_ac_c",
          "metaboliteName": "Acetate",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_atp_c",
          "metaboliteName": "ATP C10H12N5O13P3",
          "stoichiometry": "1.0"
        }
      ],
      "products": [
        {
          "metaboliteId": "M_actp_c",
          "metaboliteName": "Acetyl phosphate",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_adp_c",
          "metaboliteName": "ADP C10H12N5O10P2",
          "stoichiometry": "1.0"
        }
      ]
    }
  ],
  "isProduct": [
    {
      "id": "R_ACt2r",
      "name": "Acetate reversible transport via proton symport",
      "reversible": true,
      "reactants": [
        {
          "metaboliteId": "M_ac_e",
          "metaboliteName": "Acetate",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_h_e",
          "metaboliteName": "H+",
          "stoichiometry": "1.0"
        }
      ],
      "products": [
        {
          "metaboliteId": "M_ac_c",
          "metaboliteName": "Acetate",
          "stoichiometry": "1.0"
        },
        {
          "metaboliteId": "M_h_c",
          "metaboliteName": "H+",
          "stoichiometry": "1.0"
        }
      ]
    }
  ]
}
 @author Valerie Parra Cort�s
 *
 */
public class ReactionsOfMetabolite {
	/**
	 * The main method of class
	 * @param args[0] the path of the XML file 
	 * @param args[1] A metabolite to find the reactions
	 * @param args[2] Name of file out
	 * @throws Exception if exists any error of I/O
	 */
	public static void main(String[] args) throws Exception {
		MetabolicNetworkXMLFileLoader loader = new MetabolicNetworkXMLFileLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		//network.makeGraph();
		Map<String,List<Reaction>> reactions= network.getReactionOfMetabolite(args[1]);
		try (PrintStream out = new PrintStream(args[2])) {			
			List<Reaction> reactionsS = reactions.get("Substrates");
			out.print("{");
			if(!reactionsS.isEmpty()) {
				out.print("\"isSubstrate\":[");
				for (int i = 0; i < reactionsS.size(); i++) {
					if(i==reactionsS.size()-1) {
						out.print(reactionsS.get(i)+"");
					}
					else {
						out.print(reactionsS.get(i)+", ");
					}										
				}				
				out.print("]");
			}
			List<Reaction> products = reactions.get("Products");
			if(!products.isEmpty()) {
				out.print(",");
				out.print("\"isProduct\":[");
				for (int i = 0; i < products.size(); i++) {
					if(i==products.size()-1) {
						out.print(products.get(i)+"");
					}
					else {
						out.print(products.get(i)+", ");
					}										
				}				
				out.print("]");
			}						
			out.print("}");		
			
		}	
	}
}
