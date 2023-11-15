package metapenta.petrinet;
/**
 * Represents a edge of the petri net
 * @author Valerie Parra Cortés
 *
 */
public class Edge <O> {
	/**
	 * The stoichiometry of the metabolite
	 */
	private double tokens;

	/**
	 * The object of the edge
	 */
	private O object;

	
	public Edge(double tokens, O object) {		
		this.tokens = tokens;			
		this.object=object;
		
	}
	/**
	 * @return the stoichiometry of the edge
	 */
	public double getStoichiometry() {
		return tokens;
	}	
	/**
	 * This method changes the stoichiometry of the Edge  
	 * @param number
	 */
	public void setStoichiometry(int number) {
		this.tokens = number;
	}			
	/**
	* @return the metabolite of edge
	 */
	public O getObject() {
		return object;
	}
}