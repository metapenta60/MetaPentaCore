package metapenta.model;

/**
 * Class to create the priority queue. Represents a metabolite with a priority that is the distance to
 * any of initial reactions 
 * @author Valerie Parra
 */
class MetabolitesP implements Comparable<MetabolitesP>{
	/**
	 * The metabolite 
	 */		
	private Metabolite metabolite;
	/**
	 * The priority
	 */
	private int priority;	
	/**
	 * Constructor of the class, initializes the atributes in the value of the parameters
	 * @param metabolite The metabolite
	 * @param priority the priority
	 */
	public MetabolitesP(Metabolite metabolite, int priority) {
		super();
		this.metabolite = metabolite;
		this.priority = priority;
	}
	/**
	 * Method to get the metabolite
	 * @return themetabolite
	 */
	public Metabolite getMetabolite() {
		return metabolite;
	}

	@Override
	public int compareTo(MetabolitesP o) {		
		return o.priority-this.priority;
	}		
}