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

	//------------------------------------------------
	//.----------------Petri Net----------------------
	//------------------------------------------------

	/**
	 * A map of Integer to transition. The map represents the Transitions of Petri net 
	 */
	private  Map<Integer, Transition<Metabolite, Reaction>> transitions = new TreeMap<Integer, Transition<Metabolite,Reaction>>();

	/**
	 * A map of a String (that represents the id of the metabolite) to a Integer.
	 * This map represents the places of Petri net
	 */
	private  Map<String,Place<Metabolite, Reaction>> places = new TreeMap<String, Place<Metabolite,Reaction>>();


	private  Map<Integer,Place<Metabolite, Reaction>> placesbyNumber = new TreeMap<Integer, Place<Metabolite,Reaction>>();

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
	
	public MetabolicNetwork(Map<String,Reaction> reactions) {
		int numberMetabolites=1,numberTransition=1;		
		placesbyNumber = new TreeMap<Integer, Place<Metabolite,Reaction>>();
		transitions= new TreeMap<Integer, Transition<Metabolite, Reaction>>();		
		places = new TreeMap<String,Place<Metabolite, Reaction>>();
		Set<String> keysReaction = reactions.keySet();					
		for (String key : keysReaction) {			
			Reaction rea = reactions.get(key);
			Transition< Metabolite, Reaction> transition = new Transition< Metabolite, Reaction>(numberTransition, rea);
			numberTransition++;			
			List<ReactionComponent> reactantsC=rea.getReactants();
			List<ReactionComponent> productsC=rea.getProducts();			
			for (ReactionComponent rc : reactantsC) {			
				Metabolite meta = rc.getMetabolite();							
				Place< Metabolite, Reaction> currentPlace = places.get(meta.getId());
				
				if(currentPlace==null) {					
					Place< Metabolite, Reaction> nm = new Place< Metabolite, Reaction>(meta, numberMetabolites);
					places.put(meta.getId(),nm);
					placesbyNumber.put(numberMetabolites, nm);
					numberMetabolites++;
					currentPlace = nm;
			}				
				
				
				Edge<Place< Metabolite, Reaction>> placeCurrentMetabolite= new Edge<Place< Metabolite, Reaction>>(rc.getStoichiometry(),currentPlace);
				Edge<Transition< Metabolite, Reaction>> edgeTransition= new Edge<Transition< Metabolite, Reaction>>(rc.getStoichiometry(), transition);								
				transition.addPlaceIn(placeCurrentMetabolite);
				currentPlace.addOutTransition(edgeTransition);
				
			}			
			for (ReactionComponent rc : productsC) {
				Metabolite meta = rc.getMetabolite();				
				Place< Metabolite, Reaction> currentPlace = places.get(meta.getId());
				if(currentPlace==null) {					
					Place< Metabolite, Reaction> nm = new Place< Metabolite, Reaction>(meta, numberMetabolites);
					places.put(meta.getId(),nm);
					placesbyNumber.put(numberMetabolites, nm);
					numberMetabolites++;	
					currentPlace = nm;
				}
				
				Edge<Place< Metabolite, Reaction>> placeCurrentMetabolite= new Edge<Place< Metabolite, Reaction>>(rc.getStoichiometry(),currentPlace);
				Edge<Transition< Metabolite, Reaction>> edgeTransition= new Edge<Transition< Metabolite, Reaction>>(rc.getStoichiometry(), transition);
				transition.addPlaceOut(placeCurrentMetabolite);
				currentPlace.addInTransition(edgeTransition);
			}			
			transitions.put(transition.getNumber(), transition);						
		}	
	}

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

	//-------------------------------------------------------------------
	//-----------------------Metabolic Network --------------------------
	//--------------------------As Petri Net-----------------------------


	/**
	 * Create the Petri Net that represent the metabolic network
	 */
	public void  makeNet() {	
		int numberMetabolites=1,numberTransition=1;
		placesbyNumber = new TreeMap<Integer, Place<Metabolite,Reaction>>();
		transitions= new TreeMap<Integer, Transition<Metabolite, Reaction>>();		
		places = new TreeMap<String,Place<Metabolite, Reaction>>();


		Set<String> keysReaction = reactions.keySet();					
		for (String key : keysReaction) {			
			Reaction rea = reactions.get(key);
			Transition< Metabolite, Reaction> transition = new Transition< Metabolite, Reaction>(numberTransition, rea);
			numberTransition++;			
			List<ReactionComponent> reactantsC=rea.getReactants();
			List<ReactionComponent> productsC=rea.getProducts();			
			for (ReactionComponent rc : reactantsC) {			
				Metabolite meta = rc.getMetabolite();							
				Place< Metabolite, Reaction> currentPlace = places.get(meta.getId());

				if(currentPlace==null) {					
					Place< Metabolite, Reaction> nm = new Place< Metabolite, Reaction>(meta, numberMetabolites);
					places.put(meta.getId(),nm);
					placesbyNumber.put(numberMetabolites, nm);
					numberMetabolites++;
					currentPlace = nm;
				}				


				Edge<Place< Metabolite, Reaction>> placeCurrentMetabolite= new Edge<Place< Metabolite, Reaction>>(rc.getStoichiometry(),currentPlace);
				Edge<Transition< Metabolite, Reaction>> edgeTransition= new Edge<Transition< Metabolite, Reaction>>(rc.getStoichiometry(), transition);								
				transition.addPlaceIn(placeCurrentMetabolite);
				currentPlace.addOutTransition(edgeTransition);

			}			
			for (ReactionComponent rc : productsC) {
				Metabolite meta = rc.getMetabolite();				
				Place< Metabolite, Reaction> currentPlace = places.get(meta.getId());
				if(currentPlace==null) {					
					Place< Metabolite, Reaction> nm = new Place< Metabolite, Reaction>(meta, numberMetabolites);
					places.put(meta.getId(),nm);
					placesbyNumber.put(numberMetabolites, nm);
					numberMetabolites++;	
					currentPlace = nm;
				}

				Edge<Place< Metabolite, Reaction>> placeCurrentMetabolite= new Edge<Place< Metabolite, Reaction>>(rc.getStoichiometry(),currentPlace);
				Edge<Transition< Metabolite, Reaction>> edgeTransition= new Edge<Transition< Metabolite, Reaction>>(rc.getStoichiometry(), transition);
				transition.addPlaceOut(placeCurrentMetabolite);
				currentPlace.addInTransition(edgeTransition);
			}			
			transitions.put(transition.getNumber(), transition);						
		}				
	}

	/**
	 * Find the transition that can be triggered 
	 * @param transitions the transitions to evaluate
	 * @param metabolitesVisited a matrix where the rows are the metabolites, the first column indicate if exists the metabolite 
	 * the second column is the reaction that allows to obtain the metabolite and the third column the distance of a initial reaction 
	 * @return
	 */
	public  ArrayList<Transition< Metabolite, Reaction>> transitionsThaCanBeTriggered(List<Edge<Transition<Metabolite, Reaction>>> transitions, int[][] metabolitesVisited){
		ArrayList<Transition< Metabolite, Reaction>> transitionsThaCanBeTriggered= new ArrayList<>();
		for (int j = 0; j < transitions.size(); j++) {
			List<Edge<Place< Metabolite, Reaction>>> inEdges = transitions.get(j).getObject().getInPlaces();
			//The greatest distance that is not infinite.
			boolean canBeTriggered=true;
			for (int k = 0; k <inEdges.size(); k++) {
				int number=places.get(inEdges.get(k).getObject().getObject().getId()).getMetaboliteNumber();				
				if(metabolitesVisited[number][2]==INFINITE) {
					canBeTriggered=false;
				}
			}
			if(canBeTriggered) {
				transitionsThaCanBeTriggered.add(transitions.get(j).getObject());
			}
			canBeTriggered=false;
		}
		return transitionsThaCanBeTriggered;
	}
	/**
	 * Find the max distance of a metabolite in a specific transition
	 * @param transition the transition to be evaluated
	 * @param metabolitesVisited a matrix where the rows are the metabolites, the first column indicate if exists the metabolite 
	 * the second column is the reaction that allows to obtain the metabolite and the third column the distance of a initial reaction
	 * @return the max distance of a metabolite in the transition
	 */
	public int findMaximunDistanceOfTheInletMetabolitesOfATransition(Transition< Metabolite, Reaction> transition,int[][] metabolitesVisited) {
		List<Edge<Place< Metabolite, Reaction>>> edgesIn=transition.getInPlaces();
		int max=-1;
		for (int i = 0; i < edgesIn.size(); i++) {			
			if(metabolitesVisited[places.get(edgesIn.get(i).getObject().getObject().getId()).getMetaboliteNumber()][2]>max) {
				max=metabolitesVisited[places.get(edgesIn.get(i).getObject().getObject().getId()).getMetaboliteNumber()][2];
			}
		}		
		return max;		
	}



	/**
	 * Find all paths between a List of initial metabolties to a target metabolite
	 * @param first List o initial metabolites
	 * @param last target metabolite
	 * @param file file Where the paths will be printed
	 * @throws FileNotFoundException 
	 */
	public void getAllPaths(List<String> first, String last, String file) throws FileNotFoundException {
		int[] metabolitesVisited=new int[places.size()+1];
		for (int j = 0; j <first.size(); j++) {					
			metabolitesVisited[places.get(first.get(j)).getMetaboliteNumber()]=1;														
		}
		ArrayList<Metabolite> first_metabolites= new ArrayList<Metabolite>();
		String pathsString = "";		
		for (String m : first) {
			ArrayList<ArrayList<Integer>> paths = new ArrayList<ArrayList<Integer>>();
			ArrayList<Integer> path = new ArrayList<Integer>();
			Metabolite meta = places.get(m).getObject();	
			first_metabolites.add(meta);
			ArrayList<Metabolite> ancestors = new ArrayList<Metabolite>();
			searchPaths(meta, metabolitesVisited, path, paths, last,ancestors);
			ArrayList<ArrayList<Integer>> p_paths = new ArrayList<ArrayList<Integer>>();
			for (ArrayList<Integer> p : paths) {
				List<Edge<Place< Metabolite, Reaction>>> i_reaction= transitions.get(p.get(0)).getInPlaces();
				boolean allAtBegining=true;
				for (Edge<Place< Metabolite, Reaction>> i : i_reaction) {
					Metabolite mi = i.getObject().getObject();
					if(!first_metabolites.contains(mi)) {
						allAtBegining=false;
						break;
					}
				}
				if(allAtBegining) {
					if(!p_paths.contains(p)) {
						p_paths.add(p);
						String pv="Path: ";				
						for (int i = 0; i < p.size(); i++) {
							String r_s= transitions.get(p.get(0)).getObject().getName();
							pv+=(i==p.size()-1)?r_s+"\n":r_s+",";					
						}						
						pathsString+=pv;
					}
				}
			}
		}
		try (PrintStream out = new PrintStream(file)) {
			out.println(pathsString);
		}
	}

	@SuppressWarnings("unchecked")
	private void searchPaths(Metabolite m,int[] metabolitesVisited, ArrayList<Integer> path, ArrayList<ArrayList<Integer>> paths, String target, ArrayList<Metabolite> ancestors) {	
		if(ancestors.contains(m)) {
			return;
		}
		ancestors.add(m);
		int[] previousState = Arrays.copyOf(metabolitesVisited, metabolitesVisited.length);
		List<Transition<Metabolite, Reaction>> transitions= transitionsThaCanBeTriggeredAllPaths(places.get(m.getId()).getInTransitions(),metabolitesVisited);
		if(transitions.size()==0||m.getId().equals(target)) {			
			ArrayList<Integer> path_f = (ArrayList<Integer>) path.clone();
			paths.add(path_f);
			return;
		}
		for (Transition<Metabolite, Reaction> transition : transitions) {			
			List<Edge<Place<Metabolite, Reaction>>> out=transition.getOutPlaces();
			for (Edge<Place<Metabolite, Reaction>> edge : out) {
				Metabolite mp = edge.getObject().getObject();
				int number = places.get(mp.getId()).getMetaboliteNumber();					
				metabolitesVisited[number]=1;
			}			
		}
		for (Transition<Metabolite, Reaction> transition : transitions) {	
			path.add(transition.getNumber());
			List<Edge<Place<Metabolite, Reaction>>> out=transition.getOutPlaces();
			for (Edge<Place<Metabolite, Reaction>> edge : out) {
				Metabolite mp = edge.getObject().getObject();
				searchPaths(mp, metabolitesVisited, path, paths, target, ancestors);
			}
			path.remove(path.size()-1);
		}		

		ancestors.remove(m);
		metabolitesVisited= previousState;		
	}


	public List<Transition<Metabolite, Reaction>> transitionsThaCanBeTriggeredAllPaths(List<Edge<Transition< Metabolite, Reaction>>> list, int[] metabolitesVisited){
		ArrayList<Transition<Metabolite, Reaction>> transitionsThaCanBeTriggered= new ArrayList<>();
		for (int j = 0; j < list.size(); j++) {
			List<Edge<Place<Metabolite, Reaction>>> inEdges = list.get(j).getObject().getInPlaces();		
			boolean canBeTriggered=true;
			for (int k = 0; k <inEdges.size(); k++) {
				int number=places.get(inEdges.get(k).getObject().getObject().getId()).getMetaboliteNumber();				
				if(metabolitesVisited[number]!=1) {
					canBeTriggered=false;
					break;
				}
			}
			if(canBeTriggered) {
				transitionsThaCanBeTriggered.add(list.get(j).getObject());
			}
			canBeTriggered=false;
		}
		return transitionsThaCanBeTriggered;		
	}


	/**
	 * This method initializes the array metabolitesVisited the initial distance infinite 
	 * @param metabolitesVisited the array to initialize
	 * @param first a list with initial metabolites
	 */

	public void initialValuesOfShortestPath(int[][] metabolitesVisited,List<String> first) {
		for (int j= 0; j < metabolitesVisited.length; j++) {
			metabolitesVisited[j][2]=INFINITE;
		}
		for (int j = 0; j <first.size(); j++) {		
			System.out.println(places.get(first.get(j)).getMetaboliteNumber());
			metabolitesVisited[places.get(first.get(j)).getMetaboliteNumber()][2]=0; //Distance
			metabolitesVisited[places.get(first.get(j)).getMetaboliteNumber()][0]=1; //is there that metabolite?
			metabolitesVisited[places.get(first.get(j)).getMetaboliteNumber()][1]=-1; //The last transition, -1 if it is not assigned								
		}
	}


	public int[][] shortestPathByCatalystsNumber(List<String> first, String last, String fileName1,String fileName2, String fileName3) throws Exception{
		int[][] graph = null;		
		int[] reactionsVisited; 
		//Row 1. Moles Row 2. The reaction 3. The <<Distance>>	
		int[][] metabolitesVisited=new int[places.size()+1][3]; //Places.size()+1 for enumeration
		int[] pqea;		
		int distanceBestPath=INFINITE+1;		
		// Find for each metabolite the transitions that can be triggered				
		reactionsVisited= new int[transitions.size()+1]; 
		metabolitesVisited=new int[places.size()+1][3];			
		pqea= new int[places.size()+1];			
		//Initial values		
		initialValuesOfShortestPath(metabolitesVisited, first);
		PriorityQueue<MetabolitesP> pq = new PriorityQueue<>();	
		for (String i : first) {
			pq.add(new MetabolitesP(metabolites.get(i),0));
		}						
		while (!pq.isEmpty()) { 
			Metabolite mp = pq.poll().getMetabolite();			
			int n = places.get(mp.getId()).getMetaboliteNumber();
			if(mp.getId().equals(last)) {
				break;
			}			
			if(pqea[n]==0) {
				pqea[n]++;	

				List<Transition< Metabolite, Reaction>> transitions= transitionsThaCanBeTriggered(places.get(mp.getId()).getInTransitions(),metabolitesVisited);
				for (int j = 0; j < transitions.size(); j++) {
					reactionsVisited[transitions.get(j).getNumber()]=1;
					List<Edge<Place< Metabolite, Reaction>>> inEdgesOft = transitions.get(j).getOutPlaces();						
					//The max previous distance
					int maxValue= findMaximunDistanceOfTheInletMetabolitesOfATransition(transitions.get(j),metabolitesVisited)+transitions.get(j).getObject().getEnzymes().size();										
					//Update the distances
					for (int k = 0; k < inEdgesOft.size(); k++) {
						Metabolite meta = inEdgesOft.get(k).getObject().getObject();
						int number = places.get(meta.getId()).getMetaboliteNumber();
						if( (maxValue) < metabolitesVisited[number][2]) {
							metabolitesVisited[number][2] = maxValue;
							metabolitesVisited[number][1] = transitions.get(j).getNumber();
						}							
						MetabolitesP mp1 = new MetabolitesP(meta,(int) metabolitesVisited[number][2]);
						pq.add(mp1);							
						metabolitesVisited[number][0]=1;								
					}
				}	
			}
		}
		//Makes the metabolic Pathway
		int numberLast =places.get(last).getMetaboliteNumber();
		if(metabolitesVisited[numberLast][2]<distanceBestPath) {
			graph=metabolicPathway(metabolitesVisited, last);	
			distanceBestPath=metabolitesVisited[numberLast][2];
		}			
		printMetabolicNetworkInCSV(graph,fileName1,first,last);
		printReactionsGraphInCSV(graph,fileName2);
		printCatalystOfMethabolicPath(graph, fileName3);
		return graph;
	}	

	/**
	 * 
	 * @param first Lista metabolitos iniciales
	 * @param last Met�bolito a producir
	 * @param fileName1 Metabolic Network In CSV
	 * @param fileName2 Reactions Graph In CSV
	 * @param fileName3 Methabolic Path
	 * @throws FileNotFoundException 
	 */
	public int[][] shortestPathByMetabolitesNumber(List<String> first, String last, String fileName1,String fileName2, String fileName3) throws FileNotFoundException {
		int[][] graph = null;		
		int[] reactionsVisited; 
		//Row 1. Moles Row 2. The reaction 3. The <<Distance>>	
		int[][] metabolitesVisited=new int[places.size()+1][3]; //Places.size()+1 for enumeration
		int[] pqea;		
		int distanceBestPath=INFINITE+1;		
		// Find for each metabolite the transitions that can be triggered				
		reactionsVisited= new int[transitions.size()+1]; 
		metabolitesVisited=new int[places.size()+1][3];			
		pqea= new int[places.size()+1];			
		//Initial values		
		initialValuesOfShortestPath(metabolitesVisited, first);
		PriorityQueue<MetabolitesP> pq = new PriorityQueue<>();	
		for (String i : first) {
			pq.add(new MetabolitesP(metabolites.get(i),0));
		}						
		while (!pq.isEmpty()) { 
			Metabolite mp = pq.poll().getMetabolite();			
			int n = places.get(mp.getId()).getMetaboliteNumber();
			if(mp.getId().equals(last)) {
				break;
			}			
			if(pqea[n]==0) {
				pqea[n]++;	
				List<Transition< Metabolite, Reaction>> transitions= transitionsThaCanBeTriggered(places.get(mp.getId()).getOutTransitions(),metabolitesVisited);
				for (int j = 0; j < transitions.size(); j++) {
					reactionsVisited[transitions.get(j).getNumber()]=1;
					List<Edge<Place< Metabolite, Reaction>>> inEdgesOft = transitions.get(j).getOutPlaces();						
					//The max distance
					int maxValue= findMaximunDistanceOfTheInletMetabolitesOfATransition(transitions.get(j),metabolitesVisited)+1;
					//Update the distances
					for (int k = 0; k < inEdgesOft.size(); k++) {
						Metabolite meta = inEdgesOft.get(k).getObject().getObject();
						int number = places.get(meta.getId()).getMetaboliteNumber();
						if((maxValue)<metabolitesVisited[number][2]) {
							metabolitesVisited[number][2]=maxValue;
							metabolitesVisited[number][1]=transitions.get(j).getNumber();
						}							
						MetabolitesP mp1 = new MetabolitesP(meta,(int) metabolitesVisited[number][2]);
						pq.add(mp1);							
						metabolitesVisited[number][0]=1;								
					}
				}	
			}
		}
		//Makes the metabolic Pathway
		int numberLast =places.get(last).getMetaboliteNumber();
		if(metabolitesVisited[numberLast][2]<distanceBestPath) {
			graph=metabolicPathway(metabolitesVisited, last);	
			distanceBestPath=metabolitesVisited[numberLast][2];
		}			
		printMetabolicNetworkInCSV(graph, fileName1, first,last);
		printReactionsGraphInCSV(graph, fileName2);
		printCatalystOfMethabolicPath(graph, fileName3);
		return graph;
	}	

	
	public void printShortestPath(String fileName1, String fileName2, String fileName3, List<String> first, String last) throws FileNotFoundException {
		printMetabolicNetworkInCSV(graph, fileName1, first,last);
		printReactionsGraphInCSV(graph, fileName2);
		printCatalystOfMethabolicPath(graph, fileName3);
	}
	

	/**
	 * 
	 * @param first Lista metabolitos iniciales
	 * @param last Met�bolito a producir
	 * @throws Exception
	 */
	public Set<String> shortestPathByMetabolitesNumber(List<String> first, String last){
		Set<String> reactionss = new TreeSet<String>();
		int[] reactionsVisited; 
		//Row 1. Moles Row 2. The reaction 3. The <<Distance>>	
		int[][] metabolitesVisited=new int[places.size()+1][3]; //Places.size()+1 for enumeration
		int[] pqea;		
		int distanceBestPath=INFINITE+1;		
		// Find for each metabolite the transitions that can be triggered				
		reactionsVisited= new int[transitions.size()+1]; 
		metabolitesVisited=new int[places.size()+1][3];			
		pqea= new int[places.size()+1];			
		//Initial values		
		initialValuesOfShortestPath(metabolitesVisited, first);
		PriorityQueue<MetabolitesP> pq = new PriorityQueue<>();	
		for (String i : first) {
			pq.add(new MetabolitesP(metabolites.get(i),0));
		}						
		while (!pq.isEmpty()) { 
			Metabolite mp = pq.poll().getMetabolite();			
			int n = places.get(mp.getId()).getMetaboliteNumber();
			if(mp.getId().equals(last)) {
				break;
			}			
			if(pqea[n]==0) {
				pqea[n]++;	
				List<Transition< Metabolite, Reaction>> transitions= transitionsThaCanBeTriggered(places.get(mp.getId()).getOutTransitions(),metabolitesVisited);
				for (int j = 0; j < transitions.size(); j++) {
					reactionsVisited[transitions.get(j).getNumber()]=1;
					List<Edge<Place< Metabolite, Reaction>>> inEdgesOft = transitions.get(j).getOutPlaces();						
					//The max distance
					int maxValue= findMaximunDistanceOfTheInletMetabolitesOfATransition(transitions.get(j),metabolitesVisited)+1;
					//Update the distances
					for (int k = 0; k < inEdgesOft.size(); k++) {
						Metabolite meta = inEdgesOft.get(k).getObject().getObject();
						int number = places.get(meta.getId()).getMetaboliteNumber();
						if((maxValue)<metabolitesVisited[number][2]) {
							metabolitesVisited[number][2]=maxValue;
							metabolitesVisited[number][1]=transitions.get(j).getNumber();
						}							
						MetabolitesP mp1 = new MetabolitesP(meta,(int) metabolitesVisited[number][2]);
						pq.add(mp1);							
						metabolitesVisited[number][0]=1;								
					}
				}	
			}
		}
		//Makes the metabolic Pathway
		int numberLast =places.get(last).getMetaboliteNumber();
		if(metabolitesVisited[numberLast][2]<distanceBestPath) {
			 reactionss = reactionsMetabolicPathway(metabolitesVisited, last);	
			distanceBestPath=metabolitesVisited[numberLast][2];
			graph=metabolicPathway(metabolitesVisited, last);	
		}				
		return reactionss;
	}	

	/**
	 * Print the adjacency matrix of a graph 
	 * @param nameOfFile the path of the file
	 * @param graph the graph 
	 * @throws Exception in errors of I/O
	 */

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

	/**
	 * Wtrite
	 * @param metabolitesVisited
	 * @param last
	 * @return
	 */
	public int[][] metabolicPathway(int[][] metabolitesVisited, String last) {
		int[][] graph = new int [transitions.size()+1][transitions.size()+1];
		int[] transitionsVisited=  new int[this.transitions.size()+1];
		boolean isTheBegining=true;

		if(metabolitesVisited[places.get(last).getMetaboliteNumber()][1]==-1) {
			System.err.println("En los metabolitos de entrada ya est� el metabolito final");
		}

		else if(metabolitesVisited[places.get(last).getMetaboliteNumber()][2]==INFINITE) {
			System.err.println("No es posible llegar al metabolito final");
		}

		else {
			ArrayDeque<Integer> transitionss = new ArrayDeque<>();
			Transition< Metabolite, Reaction> currentTransition=transitions.get(metabolitesVisited[places.get(last).getMetaboliteNumber()][1]);
			do {				
				if(transitionsVisited[currentTransition.getNumber()]!=0) {	
					if(!transitionss.isEmpty()) {
						currentTransition=transitions.get(transitionss.poll());}
					else {
						currentTransition=null;
					}
					continue;
				}

				transitionsVisited[currentTransition.getNumber()]++;
				List<Edge<Place< Metabolite, Reaction>>> edgesIn= currentTransition.getInPlaces();		
				int nInitialMetabolites=0;
				for (int i = 0; i < edgesIn.size(); i++) {
					int meta = places.get(edgesIn.get(i).getObject().getObject().getId()).getMetaboliteNumber();
					int transition =metabolitesVisited[meta][1];
					if(transition!=-1) {
						graph[transition][currentTransition.getNumber()]=1;
						transitionss.add(transition);

					}
					else {
						nInitialMetabolites++;
					}
				}
				if(nInitialMetabolites==edgesIn.size()&&isTheBegining) {
					System.err.println("S�lo es necesaria una transici�n para llegar al metabolito \n"+currentTransition.getObject().getName());
					break;
				}				
				isTheBegining=false;
				if(!transitionss.isEmpty()) {
					currentTransition=transitions.get(transitionss.poll());
				}				
			}  while(currentTransition!=null);
		}
		return graph;
	}



	public Set<String> reactionsMetabolicPathway(int[][] metabolitesVisited, String last) {
		Set<String> reactions = new TreeSet<String>();
		int[][] graph = new int [transitions.size()+1][transitions.size()+1];
		int[] transitionsVisited=  new int[this.transitions.size()+1];
		boolean isTheBegining=true;

		if(metabolitesVisited[places.get(last).getMetaboliteNumber()][1]==-1) {
			System.err.println("En los metabolitos de entrada ya est� el metabolito final");
		}

		else if(metabolitesVisited[places.get(last).getMetaboliteNumber()][2]==INFINITE) {
			System.err.println("No es posible llegar al metabolito final");
		}

		else {
			ArrayDeque<Integer> transitionss = new ArrayDeque<>();
			Transition< Metabolite, Reaction> currentTransition=transitions.get(metabolitesVisited[places.get(last).getMetaboliteNumber()][1]);
			do {				
				if(transitionsVisited[currentTransition.getNumber()]!=0) {	
					if(!transitionss.isEmpty()) {
						currentTransition=transitions.get(transitionss.poll());						
					}
					else {
						currentTransition=null;
					}
					continue;
				}
				transitionsVisited[currentTransition.getNumber()]++;
				reactions.add(currentTransition.getObject().getId());
				List<Edge<Place< Metabolite, Reaction>>> edgesIn= currentTransition.getInPlaces();		
				int nInitialMetabolites=0;
				for (int i = 0; i < edgesIn.size(); i++) {
					int meta = places.get(edgesIn.get(i).getObject().getObject().getId()).getMetaboliteNumber();
					int transition =metabolitesVisited[meta][1];
					if(transition!=-1) {
						graph[transition][currentTransition.getNumber()]=1;
						transitionss.add(transition);
						reactions.add(transitions.get(transition).getObject().getId());	
						
					}
					else {
						nInitialMetabolites++;
					}
				}
				if(nInitialMetabolites==edgesIn.size()&&isTheBegining) {
					System.err.println("S�lo es necesaria una transici�n para llegar al metabolito \n"+currentTransition.getObject().getName());
					break;
				}				
				isTheBegining=false;
				if(!transitionss.isEmpty()) {
					currentTransition=transitions.get(transitionss.poll());
				}				
			}  while(currentTransition!=null);
		}
		return reactions;
	}

	public void printReactionsGraphInCSV(int[][] graph, String fileName) throws FileNotFoundException{		
		try (PrintStream out = new PrintStream(fileName)) {
			out.println("source,target,interaction,directed,symbol,value");
			for (int i = 0; i < graph.length; i++) {
				for (int j = 0; j < graph.length; j++) {
					if(graph[i][j]!=0) {
						String nameTransition1= transitions.get(i).getObject().getId();
						String nameTransition2= transitions.get(j).getObject().getId();						
						out.println(nameTransition1+COMMA+nameTransition2+",PP,TRUE,abc"+i+j+",1.234");						
					}						
				}

			}
		}		
	}




	/**
	 * Buscar varias rutas
	 * jduitama
	 * @param graph
	 * @param fileName
	 * @throws Exception
	 */
	private void printCatalystOfMethabolicPath(int[][] graph, String fileName) throws FileNotFoundException{	
		Set<GeneProduct> allCatalyst = new TreeSet<GeneProduct>();
		try (PrintStream out = new PrintStream(fileName)) {
			out.println("Catalys of metabolic pathway");
			for (int i = 0; i < graph.length; i++) {
				for (int j = 0; j < graph.length; j++) {
					if(graph[i][j]!=0) {
						List<GeneProduct> catalysReaction1= transitions.get(i).getObject().getEnzymes();
						List<GeneProduct> catalysReaction2= transitions.get(j).getObject().getEnzymes();
						allCatalyst.addAll(catalysReaction1);
						allCatalyst.addAll(catalysReaction2);
					}						
				}
			}		
		}		
	}


	/**
	 * Print the graph entered by parameter in a CSV
	 * @param graph The graph to be printed, the graph represents the adyacent matrix of reactions
	 * @param fileName the path of CSV
	 * @param first The ids of initial metabolites
	 * @param last the final metabolite
	 * @throws FileNotFoundException 
	 * @throws Exception in errors of I/O
	 */
	private void printMetabolicNetworkInCSV(int[][] graph, String fileName, List<String> first, String last) throws FileNotFoundException{
		int[] metabolitesVisited= new int[places.size()];
		int[] transitionsVisited= new int[transitions.size()];
		try (PrintStream out = new PrintStream(fileName)) {
			out.println("source,target,interaction,directed,symbol,value");
			for (int i = 0; i < graph.length; i++) {
				for (int j = 0; j < graph.length; j++) {
					if(graph[i][j]!=0) {
						printATransitionInCSV(transitions.get(i),out,i,j,metabolitesVisited,transitionsVisited,first,last);
						printATransitionInCSV(transitions.get(j),out,i,j,metabolitesVisited,transitionsVisited,first,last);			
					}			
				}
			}
		}
	}

	/**
	 * Print the graph entered by parameter in a CSV
	 * @throws Exception in errors of I/O
	 */
	public void printAllMetabolicNetworkInCSV(String fileName) throws IOException{
		try (PrintStream out = new PrintStream(fileName)) {
			out.println("source,target,interaction,directed,symbol,value");
			Set<Integer> keys= transitions.keySet();
			for (Integer integer : keys) {				
				printATransitionInCSVToAllNetwork(transitions.get(integer),out);
			}	
		}
	}


	/**
	 * Method that print the transition in the CSV entered by parameter 
	 * @param t the transiton to be printed
	 * @param out the Stream of the CSV
	 */
	private void printATransitionInCSVToAllNetwork(Transition< Metabolite, Reaction> t, PrintStream out) {
		List<Edge<Place< Metabolite, Reaction>>> metaInTransition = t.getInPlaces();		
		for (Edge<Place< Metabolite, Reaction>> edge : metaInTransition) {
			String nameOfMeta = edge.getObject().getObject().getId();		
			out.println(nameOfMeta+COMMA+t.getObject().getId()+",PP,TRUE,abc"+",1.234");
		}
		metaInTransition = t.getOutPlaces();		
		for (Edge<Place< Metabolite, Reaction>> edge : metaInTransition) {
			String nameOfMeta = edge.getObject().getObject().getId();		
			out.println(t.getObject().getId()+COMMA+nameOfMeta+",PP,TRUE,abc"+",1.234");
		}
	}	

	/**
	 * Method that print the transition in the CSV entered by parameter 
	 * @param t the transiton to be printed
	 * @param out the Stream of the CSV
	 * @param i A identifier
	 * @param j A identifier
	 * @param metabolitesVisited a matrix where the rows are the metabolites, the first column indicate if exists the metabolite 
	 * the second column is the reaction that allows to obtain the metabolite and the third column the distance of a initial reaction
	 * @param transitionsVisited the transition that were printed 
	 * @param first The initial ids of metabolites
	 * @param last The final metabolite
	 */
	private void printATransitionInCSV(Transition< Metabolite, Reaction> t, PrintStream out,int i, int j,int[] metabolitesVisited,int[] transitionsVisited, List<String> first, String last) {
		List<Edge<Place< Metabolite, Reaction>>> metaInTransition = t.getInPlaces();
		if(transitionsVisited[t.getNumber()]==0) {
			transitionsVisited[t.getNumber()]++;
			for (Edge<Place< Metabolite, Reaction>> edge : metaInTransition) {
				String nameOfMeta = edge.getObject().getObject().getId();
				if(first.contains(nameOfMeta)) {
					nameOfMeta="TT"+nameOfMeta;
				}
				else if(last.equalsIgnoreCase(nameOfMeta)) {
					nameOfMeta="OO"+nameOfMeta;
				}
				out.println(nameOfMeta+COMMA+t.getObject().getId()+",PP,TRUE,abc"+i+j+",1.234");
			}
			metaInTransition = t.getOutPlaces();		
			for (Edge<Place< Metabolite, Reaction>> edge : metaInTransition) {
				String nameOfMeta = edge.getObject().getObject().getId();
				if(first.contains(nameOfMeta)) {
					nameOfMeta="TT"+nameOfMeta;
				}
				else if(last.equalsIgnoreCase(nameOfMeta)) {
					nameOfMeta="OO"+nameOfMeta;
				}
				out.println(t.getObject().getId()+COMMA+nameOfMeta+",PP,TRUE,abc"+i+j+",1.234");
			}	
		}
	}	

	/**
	 * Find the sinks of the net
	 * @return List with sinks metabolites
	 */
	public List<Metabolite> findSinks() {	
		List<Metabolite> sinks= new ArrayList<Metabolite>();
		for (String key : metabolites.keySet()) {
			Metabolite m = metabolites.get(key);
			if(places.get(m.getId()).getOutTransitions().size()==0) {
				sinks.add(m);
			}
		}
		return sinks;
	}

	/**
	 * Find the sinks of the net
	 * @return List with sinks metabolites
	 */
	public List<Metabolite> findSources() {	
		List<Metabolite> sinks= new ArrayList<Metabolite>();
		for (String key : metabolites.keySet()) {
			Metabolite m = metabolites.get(key);
			if(places.get(m.getId()).getInTransitions().size()==0) {
				sinks.add(m);
			}
		}
		return sinks;
	}


	//-------------------------------------------------
	//---------------Inner auxiliar classes------------
	//-------------------------------------------------


	/**
	 * Construct the correspondent sub-net of the comparment 
	 * @param prefixOut the output prefix
	 */


	public void describeNet(String prefixOut) throws Exception{
		DescribeNetworkWriter gfw = new DescribeNetworkWriter(prefixOut);

		List<String> metaboliteIds = getMetaboliteIds();
		gfw.writeMetabolites(metaboliteIds);

		List<String> eMetaboliteIds = getEMetaboliteIds();
		gfw.writeEMetabolites(eMetaboliteIds);

		List<String> cMetaboliteIds = getCMetaboliteIds();
		gfw.writeCMetabolites(cMetaboliteIds);

		List<String> reactionIds = getReactionIds();
		gfw.writeReactions(reactionIds);

		List<String> reversibleReactionsIds = getReversibleReactionsIds();
		gfw.writeReversibleReactions(reversibleReactionsIds);

		List<String> irreversibleReactionsIds = getIrreversibleReactionsIds();
		gfw.writeIrreversibleReactions(irreversibleReactionsIds);

		Set<Integer> keys = transitions.keySet();
		for (Integer key: keys){
			Transition<Metabolite,Reaction> transition = transitions.get(key);
			Reaction reaction = transition.getObject();
			List<ReactionComponent> reactants = reaction.getReactants();
			for (ReactionComponent reactant: reactants){
				Metabolite reactantMetabolite = reactant.getMetabolite();
;				double reactantStoichiometry = - reactant.getStoichiometry();
				gfw.WriteInSMatrix(reactantMetabolite.getId(), reactantStoichiometry);
			}

			List<ReactionComponent> products = reaction.getProducts();
			for (ReactionComponent product: products){
				Metabolite productMetabolite = product.getMetabolite();
				gfw.WriteInSMatrix(productMetabolite.getId(), product.getStoichiometry());
			}
		}

		gfw.Write();
	}

	public List<String> getMetaboliteIds(){
		List<String> metabolitesIds = new ArrayList();
		Set<String> metabolitesKeys = places.keySet();

		for(String key: metabolitesKeys) {
			Metabolite metabolite = places.get(key).getObject();
			metabolitesIds.add(metabolite.getId());
		}

		return metabolitesIds;
	}

	public List<String> getEMetaboliteIds(){
		List<String> metabolitesIds = new ArrayList();
		Set<String> metabolitesKeys = places.keySet();

		for(String key: metabolitesKeys) {
			Metabolite metabolite = places.get(key).getObject();
			if (metabolite.getCompartment().equals("e")){
				metabolitesIds.add(metabolite.getId());
			}
		}

		return metabolitesIds;
	}

	public List<String> getCMetaboliteIds(){
		List<String> metabolitesIds = new ArrayList();
		Set<String> metabolitesKeys = places.keySet();

		for(String key: metabolitesKeys) {
			Metabolite metabolite = places.get(key).getObject();
			if (!metabolite.getCompartment().equals("e")){
				metabolitesIds.add(metabolite.getId());
			}
		}

		return metabolitesIds;
	}

	public List<String> getReactionIds(){
		List<String> reactionIds = new ArrayList();
		Set<String> keys = reactions.keySet();

		for(String key: keys) {
			reactionIds.add(key);
		}

		return reactionIds;
	}


	public List<String> getReversibleReactionsIds(){
		List<String> reactionIds = new ArrayList();
		Set<Integer> keys = transitions.keySet();

		for(Integer key: keys) {
			Reaction reaction = transitions.get(key).getObject();
			if(reaction.isReversible()){
				reactionIds.add(reaction.getId());
			}
		}

		return reactionIds;
	}

	public List<String> getIrreversibleReactionsIds(){
		List<String> reactionIds = new ArrayList();
		Set<Integer> keys = transitions.keySet();

		for(Integer key: keys) {
			Reaction reaction = transitions.get(key).getObject();
			if(!reaction.isReversible()){
				reactionIds.add(reaction.getId());
			}
		}

		return reactionIds;
	}

	public Map<String, Integer> connectedComponents() {				
		int[] metabolitesVisited = new int[places.size()+1];
		int[] transitionsVisited = new int[transitions.size()+1];
		int conectedComponentId =1;
		int tiempo=0;

		Set<String> keys = places.keySet();

		for(String node: keys) {
			if(metabolitesVisited[places.get(node).getMetaboliteNumber()]==0) {
				visitNode(metabolitesVisited, transitionsVisited, node, tiempo, conectedComponentId);
				conectedComponentId++;
			}			
		}
		Map<String, Integer> connectedComponents = new TreeMap<String, Integer>();
		Set<Integer> transitionsKeys = transitions.keySet();

		for (Integer integer : transitionsKeys) {
			connectedComponents.put(transitions.get(integer).getObject().getName(), transitionsVisited[integer]);
		}				
		return connectedComponents;
	}

	private void visitNode(int[] metabolitesVisited, int[] transitionsVisited, String node,  int tiempo, int conectedComponentId) {
		metabolitesVisited[places.get(node).getMetaboliteNumber()] = 1;		
		tiempo = tiempo+1;					
		List<Edge<Transition< Metabolite, Reaction>>> neighbours = places.get(node).getInTransitions();
		for (int i = 0; i < neighbours.size(); i++) {
			Edge<Transition< Metabolite, Reaction>> reaction = neighbours.get(i);
			transitionsVisited[reaction.getObject().getNumber()]=conectedComponentId;
			List<Edge<Place< Metabolite, Reaction>>> in= reaction.getObject().getInPlaces();
			List<Edge<Place< Metabolite, Reaction>>> out=reaction.getObject().getOutPlaces();
			for (int j = 0; j < in.size(); j++) {
				String node_1=in.get(j).getObject().getObject().getId();
				if(metabolitesVisited[places.get(node_1).getMetaboliteNumber()]==0) {
					visitNode(metabolitesVisited, transitionsVisited, node_1, tiempo, conectedComponentId);
				}				
			}
			for (int j = 0; j < out.size(); j++) {
				String node_2=out.get(j).getObject().getObject().getId();
				if(metabolitesVisited[places.get(node_2).getMetaboliteNumber()]==0) {
					visitNode(metabolitesVisited, transitionsVisited, node_2, tiempo, conectedComponentId);
				}
			}

		}		
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

	public Map<Integer, Transition<Metabolite, Reaction>> getTransitions() {
		return transitions;
	}

	public Map<String, Place<Metabolite, Reaction>> getPlaces() {
		return places;
	}

	public Map<Integer,Place<Metabolite, Reaction>> getPlacesbyNumber() {
		return placesbyNumber;
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
	 
	public void printsSinksInAFile(String filename) throws FileNotFoundException {		
		List<Metabolite> sinks = findSinks();
		try (PrintStream out = new PrintStream(filename)) {
			out.print("{sinks:[");
			for (int i = 0; i < sinks.size(); i++) {
				if(i==sinks.size()-1) {
					out.print(sinks.get(i).toString());
				}
				else {
					out.print(sinks.get(i).toString()+",");
				}				
			}
			out.print("]}");			
		}		
	}
	
	public void printsSourcesInAFile(String filename) throws FileNotFoundException {		
		List<Metabolite> sinks = findSources();
		try (PrintStream out = new PrintStream(filename)) {
			out.print("{sinks:[");
			for (int i = 0; i < sinks.size(); i++) {
				if(i==sinks.size()-1) {
					out.print(sinks.get(i).toString());
				}
				else {
					out.print(sinks.get(i).toString()+",");
				}				
			}
			out.print("]}");			
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
	
	public void printConnectedComponents(String filename) throws IOException {		
		Map<String, Integer> connectedComponents = connectedComponents();
		StringBuilder csv = new StringBuilder();
		Set<String> reactionsSet = connectedComponents.keySet();
		for (String reaction : reactionsSet) {
			csv.append(reaction+","+connectedComponents.get(reaction)+"\n");
		}
		Files.write(Paths.get(filename), csv.toString().getBytes());
	}
	
}