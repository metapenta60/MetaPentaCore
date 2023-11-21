package metapenta.model;

/**
 * Represents a metabolite that participates in chemical reactions
 * @author Jorge Duitama
 */
public class Metabolite {
	private String id;

	private String name;

	private String compartment;

	private String chemicalFormula;

	public Metabolite(String id, String name, String compartment) {
		super();
		this.id = id;
		this.name = name;
		this.compartment = compartment;
	}

	public String getChemicalFormula() {
		return chemicalFormula;
	}

	public void setChemicalFormula(String chemicalFormula) {
		this.chemicalFormula = chemicalFormula;
	}

	public String getId() {
		return id;
	}

	public String getName() {
		return name;
	}

	public String getCompartment() {
		return compartment;
	}

}
