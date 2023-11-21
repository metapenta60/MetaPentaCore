package metapenta.model;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import metapenta.petrinet.Edge;
import metapenta.petrinet.Place;
import metapenta.petrinet.Transition;
import metapenta.tools.io.DescribeNetworkWriter;

/**
 * Represents a metabolic network of reactions on metabolites
 * @author Jorge Duitama
 */
public class MetabolicNetwork {
	private Map<String,GeneProduct> geneProducts = new TreeMap<String,GeneProduct>();
	private String name;

	/**
	 * Metabolites of  the metabolic Network
	 */
	private Map<String,Metabolite> metabolites = new TreeMap<String,Metabolite>();
	private Set<String> compartments = new TreeSet<String>();
	private Map<String,Reaction> reactions = new TreeMap<String,Reaction>();

	public void setName(String name) {
		this.name = name;
	}
	
	
	/**
	 * Represents the infinite distance between 2 metabolites
	 */
	public static final int INFINITE=100000;
	/**
	 * String with the value of comma (,) to generate the csv
	 */
	public static final String COMMA=",";
	/**F
	 * Adds a new gene product that can catalyze reactions
	 * @param product New gene product
	 */
	public void addGeneProduct(GeneProduct product) {
		geneProducts.put(product.getId(), product);
	}
	/**
	 * Adds a new metabolite. If a metabolite with the given name is already added, it 
	 * @param metabolite New metabolite
	 */
	public void addMetabolite(Metabolite metabolite) {
		metabolites.put(metabolite.getId(), metabolite);
		compartments.add(metabolite.getCompartment());
	}
	/**
	 * Adds a new reaction
	 * @param r New reaction between metabolites
	 */
	public void addReaction(Reaction r) {
		reactions.put(r.getId(),r);
	}
	/**
	 * Returns the gene product with the given id
	 * @param id of the product to search
	 * @return GeneProduct with the given id
	 */
	public GeneProduct getGeneProduct (String id) {
		return geneProducts.get(id);
	}
	/**
	 * Returns the metabolite with the given id
	 * @param id of the metabolite to search
	 * @return Metabolite with the given id
	 */
	public Metabolite getMetabolite (String id) {
		return metabolites.get(id);
	}


	public Reaction getReaction(String id) {
		return reactions.get(id);
	}



	/**
	 * @return List of metabolites in the network
	 */
	public List<Metabolite> getMetabolitesAsList() {
		return new ArrayList<Metabolite>(metabolites.values());
	}


	/**
	 * @return List of metabolites in the network
	 */
	public List<String> getMetabolitesAsListString() {
		List<String> metabolitesString = new ArrayList<String>();
		for (String string : metabolites.keySet()) {
			metabolitesString.add(string);
		}

		return metabolitesString;
	}

	/**
	 * @return List of reactions in the network
	 */
	public List<Reaction> getReactionsAsList () {
		return new ArrayList<Reaction>(reactions.values());
	}
	
	public MetabolicNetwork() {
		
	}
	
	private int[][] graph;


	//-------------------------------------------------------------------
	//-----------------------Metabolic Network --------------------------
	//---------------------------As graph--------------------------------

	/**
	 * Method that find the reactions where a metabolite  
	 * @param metaboliteKeyName the id of the metabolite
	 * @return A map of the where reaction where metabolite is a substrate and where the metabolite is a product 
	 */
	public Map<String,List<Reaction>> getReactionOfMetabolite(String metaboliteKeyName) {
		Map<String,List<Reaction>> reaction= new TreeMap<>();	
		List<Reaction> rsubstrates= new ArrayList<>();
		List<Reaction> rproducts=new ArrayList<>();
		Set<String> keys=reactions.keySet();		
		for (String key : keys) {
			Reaction rea= reactions.get(key);
			List<ReactionComponent> substrates= rea.getReactants();
			List<ReactionComponent> products= rea.getProducts();			
			for (int i = 0; i < substrates.size(); i++) {
				if(substrates.get(i).getMetabolite().getId().equals(metaboliteKeyName)) {
					rsubstrates.add(rea);					
				}
			}
			for (int i = 0; i < products.size(); i++) {
				if(products.get(i).getMetabolite().getId().equals(metaboliteKeyName)) {
					rproducts.add(rea);					
				}
			}
		}		
		reaction.put("Substrates", rsubstrates);
		reaction.put("Products", rproducts);
		return reaction;
	}
	/**
	 * Find the reactions where the enzyme is the catalyst
	 * @param enzymeName the id of Enzume
	 * @return a mapa whith the reactions
	 */

	public List<Reaction> getReactionsCatalyzedBy(String enzymeName){
		List<Reaction> cata= new ArrayList<Reaction>();
		Set<String> keys=reactions.keySet();		
		for (String key : keys) {
			Reaction rea= reactions.get(key);
			List<GeneProduct> enzymes=rea.getEnzymes();
			for (GeneProduct enzyme : enzymes) {				
				if(enzyme.getName().equals(enzymeName)) {
					cata.add(rea);
					break;
				}
			}			
		}
		return cata;
	}

	public void writeGraph(String nameOfFile,int[][] graph) throws Exception {
		try (PrintStream out = new PrintStream(nameOfFile)) {
			for (int i = 0; i < graph.length; i++) {
				for (int j = 0; j < graph.length; j++) {
					out.print(graph[i][j]);
				}
				out.print("\n");
			}
		}
	}



	public List<String> getReactionIds(){
		List<String> reactionIds = new ArrayList();
		Set<String> keys = reactions.keySet();

		for(String key: keys) {
			reactionIds.add(key);
		}

		return reactionIds;
	}

	public List<Metabolite> commonMetabolites(MetabolicNetwork mn2){
		List<Metabolite> commonMetabolites = new ArrayList<Metabolite>();
		List<Metabolite> metabolitesNetwork2 = mn2.getMetabolitesAsList();		
		for (Metabolite metabolite : metabolitesNetwork2) {	
			Metabolite place = this.metabolites.get(metabolite.getId());
			if(place!=null) {
				commonMetabolites.add(metabolite);
			}
		}		
		return commonMetabolites;
	}

	public List<Reaction> commonReactions(MetabolicNetwork mn2){
		List<Reaction> commonReactions = new ArrayList<Reaction>();
		List<Reaction> metabolitesNetwork2 = mn2.getReactionsAsList();		
		for (Reaction reaction : metabolitesNetwork2) {
			Reaction transicion = reactions.get(reaction.getId());
			if(transicion!=null) {
				commonReactions.add(reaction);
			}
		}		
		return commonReactions;
	}

	
	public Map<String,String> reactiosnMetboliteString(String metabolite) {
		Map<String,List<Reaction>> reactions= getReactionOfMetabolite(metabolite);
		Map<String, String> mapa = new TreeMap<String, String>();
		List<Reaction> reactionsS = reactions.get("Substrates");
		String isSubstrate ="";
		for (int i = 0; i < reactionsS.size(); i++) {
			isSubstrate+=reactionsS.get(i).getName()+"\n";
		}
		List<Reaction> productsS = reactions.get("Products");
		String isProduct ="";
		for (int i = 0; i < productsS.size(); i++) {
			isProduct+=productsS.get(i).getName()+"\n";
		}
		
		mapa.put("Substrates", isSubstrate);
		mapa.put("Products", isProduct);
		
		return mapa;
		
	}
	public void printInAFileReactionsOfMetabolite(String metabolite, String fileName) throws FileNotFoundException {
		Map<String,List<Reaction>> reactions= getReactionOfMetabolite(metabolite);
		try (PrintStream out = new PrintStream(fileName)) {			
			List<Reaction> reactionsS = reactions.get("Substrates");
			out.print("{");
			if(!reactionsS.isEmpty()) {
				out.print("\"isSubstrate\":[");
				for (int i = 0; i < reactionsS.size(); i++) {
					if(i==reactionsS.size()-1) {
						out.print(reactionsS.get(i)+"");
					}
					else {
						out.print(reactionsS.get(i)+", \n");
					}										
				}				
				out.print("]\n");
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
						out.print(products.get(i)+", \n");
					}										
				}				
				out.print("] \n");
			}						
			out.print("}");		
			
		}	
	}
	 

	
	public String getEnzymesAsString(String reaction) {
		List<GeneProduct> enzymes = reactions.get(reaction).getEnzymes();
		String enzymesString ="";
		for (int i = 0; i < enzymes.size(); i++) {
			enzymesString+=enzymes.get(0).getName()+"\n";
		}
		return enzymesString;
	}

}