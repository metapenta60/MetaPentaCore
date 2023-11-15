package metapenta.commands;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import metapenta.model.MetabolicNetwork;
import metapenta.model.MetabolicNetworkXMLLoader;
import metapenta.model.Metabolite;
import metapenta.model.Reaction;
import metapenta.model.ReactionComponent;


public class GapFindGams {

	/**
	 * Class to test the generation of gap find GAMS script 
	 * @param args[0] xml file of metabolic network
	 * @throws Exception if exists any error of I/O
	 */
	
	
	public static void main(String[] args) throws IOException {
		MetabolicNetworkXMLLoader loader = new MetabolicNetworkXMLLoader();
		MetabolicNetwork network = loader.loadNetwork(args[0]);
		List<Metabolite> metabolites = network.getMetabolitesAsList();
		PrintStream compounds = new PrintStream("compounds.txt") ;
		PrintStream extracellular_compounds = new PrintStream("extracellular_compounds.txt") ;
		PrintStream cytosolic_compounds = new PrintStream("cytosolic_compounds.txt") ;
		compounds.println("/");
		cytosolic_compounds.println("/");
		extracellular_compounds.print("/\n");
		for (int i = 0; i < metabolites.size(); i++) {
			Metabolite currentMetabolite = metabolites.get(i);
			compounds.println("'"+currentMetabolite.getId()+"["+currentMetabolite.getCompartment()+"]'");
			if(currentMetabolite.getCompartment().equals("e")) {
				extracellular_compounds.print("'"+currentMetabolite.getId()+"[e]'\n");
			}
			if(currentMetabolite.getCompartment().equals("c")) {
				cytosolic_compounds.print("'"+currentMetabolite.getId()+"[c]'\n");
			}
		}
		
		List<Reaction> reactions_list = network.getReactionsAsList();
		PrintStream lowerbound_on_fluxes = new PrintStream("lowerbound_on_fluxes.txt") ;
		PrintStream upperbound_on_fluxes = new PrintStream("upperbound_on_fluxes.txt") ;
		PrintStream reactions = new PrintStream("reactions.txt") ;	
		PrintStream reversible_reactions = new PrintStream("reversible_reactions.txt") ;
		PrintStream s_matrix = new PrintStream("S_matrix.txt") ;		
		lowerbound_on_fluxes.print("/\n");
		upperbound_on_fluxes.print("/\n");
		reactions.print("/\n");
		s_matrix.print("/\n");
		reversible_reactions.print("/\n");
		for (int i = 0; i < reactions_list.size(); i++) {
			Reaction currentReaction = reactions_list.get(i);
			int lo = (int) currentReaction.getLowerBoundFlux();
			int up = (int) currentReaction.getUpperBoundFlux();
			lowerbound_on_fluxes.print("'"+currentReaction.getId()+"' "+lo+"\n");
			upperbound_on_fluxes.print("'"+currentReaction.getId()+"' "+up+"\n");
			reactions.print("'"+currentReaction.getId()+"'\n");
			if(currentReaction.isReversible()) {
				reversible_reactions.print("'"+currentReaction.getId()+"'\n");
			}
			
			List<ReactionComponent> reactants =currentReaction.getReactants();
			for (int j = 0; j < reactants.size(); j++) {
				ReactionComponent rc = reactants.get(j);		
				Integer s = (int) rc.getStoichiometry();
				s_matrix.print("'"+rc.getMetabolite().getId()+"["+rc.getMetabolite().getCompartment()+"]'.'"+currentReaction.getId()+"' -"+s+"\n");
			}
			List<ReactionComponent> products =currentReaction.getProducts();
			for (int j = 0; j < products.size(); j++) {
				ReactionComponent rc = products.get(j);
				Integer s = (int) rc.getStoichiometry();
				s_matrix.print("'"+rc.getMetabolite().getId()+"["+rc.getMetabolite().getCompartment()+"]'.'"+currentReaction.getId()+"' "+s+"\n");
			}
			
		}
				
		lowerbound_on_fluxes.print("/");
		lowerbound_on_fluxes.close();
		upperbound_on_fluxes.print("/");
		upperbound_on_fluxes.close();
		reactions.print("/");
		reactions.close();
		reversible_reactions.print("/");
		reversible_reactions.close();
		s_matrix.print("/");
		s_matrix.close();
		compounds.print("/");
		compounds.close();
		extracellular_compounds.print("/");
		extracellular_compounds.close();
		cytosolic_compounds.print("/");
		cytosolic_compounds.close();
	}
}
