package metapenta.petrinet;
import java.util.ArrayList;
import java.util.List;

import metapenta.model.GeneProduct;
/**
 * Class that represent a transition
 * @author Valerie Parra
 */
public class Transition <PlaceObjectClass, TransitionsObjectClass>{
	/**
	 * The id of the transition
	 */
	private int transitionNumber;
	/**
	 * The list of edges in
	 */
	private List<Edge<Place<PlaceObjectClass, TransitionsObjectClass>>> inPlaces;
	/**
	 * The list of edges out
	 */
	private List<Edge<Place<PlaceObjectClass, TransitionsObjectClass>>> outPlaces;

	
	private TransitionsObjectClass object;
	
	/**
	 * The constructor of class
	 * @param gp the list of the enzymes of the reaction
	 * @param number the id (number)
	 * @param name of the reaction
	 * @param id (the string)
	 */
	public Transition(int number, TransitionsObjectClass object) {
		this.inPlaces = new ArrayList<>();
		this.outPlaces = new ArrayList<>();		
		this.object = object;
		this.transitionNumber=number;

	}	
	/**
	 * Add a edges in the list of edges in
	 * @param edge to add
	 */
	public void addPlaceIn(Edge<Place<PlaceObjectClass, TransitionsObjectClass>> edge) {
		inPlaces.add(edge);
	}	
	/**
	 * Add a edges in the list of edges out
	 * @param edge to add
	 */
	public void addPlaceOut(Edge<Place<PlaceObjectClass, TransitionsObjectClass>> edge) {
		outPlaces.add(edge);
	}
	
	/**	
	 * @return The in edges of the transition
	 */
	public List<Edge<Place<PlaceObjectClass, TransitionsObjectClass>>> getInPlaces() {
		return inPlaces;
	}

	/**
	 * @return The in edges of the transition
	 */
	public List<Edge<Place<PlaceObjectClass, TransitionsObjectClass>>> getAllPLaces() {
		ArrayList allPlaces = new ArrayList<>();
		for (Edge<Place<PlaceObjectClass, TransitionsObjectClass>> place: inPlaces){
			allPlaces.add(place);
		}
		for (Edge<Place<PlaceObjectClass, TransitionsObjectClass>> place: outPlaces){
			allPlaces.add(place);
		}

		return allPlaces;
	}

	/**	
	 * @return The out edges of the transition
	 */
	
	public List<Edge<Place<PlaceObjectClass, TransitionsObjectClass>>> getOutPlaces() {
		return outPlaces;
	}	

	/**
	* Set the number (id) of the transition
	* @param number the new number (id) of transition
	**/
	public void setNumber(int number) {
		this.transitionNumber = number;
	}
	
	/**
	* @return the number (id) of this transition
	**/	
	public int getNumber() {
		return transitionNumber;
	}	
	
	public TransitionsObjectClass getObject() {
		return object;
	}
}