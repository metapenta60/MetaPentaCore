package metapenta.tools.io;

import metapenta.model.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class MetabolicNetworkXMLFileLoader {
    protected static final String ELEMENT_MODEL = "model";
    protected static final String ELEMENT_LISTPRODUCTS = "fbc:listOfGeneProducts";
    protected static final String ELEMENT_LISTMETABOLITES = "listOfSpecies";
    protected static final String ELEMENT_LISTREACTIONS = "listOfReactions";
    protected static final String ELEMENT_GENEPRODUCT = "fbc:geneProduct";
    protected static final String ATTRIBUTE_FBCID = "fbc:id";
    protected static final String ATTRIBUTE_FBCNAME = "fbc:name";
    protected static final String ATTRIBUTE_FBCLABEL = "fbc:label";
    protected static final String ELEMENT_METABOLITE = "species";
    protected static final String ATTRIBUTE_ID = "id";
    protected static final String ATTRIBUTE_NAME = "name";
    protected static final String ATTRIBUTE_COMPARTMENT = "compartment";
    protected static final String ATTRIBUTE_FBCFORMULA = "fbc:chemicalFormula";
    protected static final String ELEMENT_REACTION = "reaction";
    protected static final String ATTRIBUTE_REVERSIBLE = "reversible";
    protected static final String ELEMENT_GENEASSOC = "fbc:geneProductAssociation";
    protected static final String ELEMENT_LISTREACTANTS = "listOfReactants";
    protected static final String ELEMENT_LISTMETABPRODUCTS = "listOfProducts";
    protected static final String ELEMENT_METABREF = "speciesReference";
    protected static final String ATTRIBUTE_STOICHIOMETRY = "stoichiometry";
    protected static final String ELEMENT_GENEPRODUCTREF = "fbc:geneProductRef";
    protected static final String ATTRIBUTE_FBC_LOWERBOUND = "fbc:lowerFluxBound";
    protected static final String ATTRIBUTE_FBC_UPPERBOUND = "fbc:upperFluxBound";

    protected void loadGeneProducts(Element listElem, MetabolicNetwork network) throws IOException {
        NodeList offspring = listElem.getChildNodes();
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_GENEPRODUCT.equals(elem.getNodeName())) {
                    String id = elem.getAttribute(ATTRIBUTE_FBCID);
                    if(id==null || id.length()==0) throw new IOException("Every gene product should have an id");
                    String name = elem.getAttribute(ATTRIBUTE_FBCNAME);
                    if(name==null || name.length()==0) name = id;
                    String label = elem.getAttribute(ATTRIBUTE_FBCLABEL);


                    GeneProduct product = new GeneProduct(id, name);
                    if(label!=null) product.setLabel(label);
                    network.addGeneProduct (product);
                }
            }
        }
    }

    protected void loadMetabolites(Element listElem, MetabolicNetwork network) throws IOException {
        NodeList offspring = listElem.getChildNodes();
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_METABOLITE.equals(elem.getNodeName())) {
                    String id = elem.getAttribute(ATTRIBUTE_ID);
                    if(id==null || id.length()==0) throw new IOException("Every metabolite should have an id");
                    String name = elem.getAttribute(ATTRIBUTE_NAME);
                    if(name==null || name.length()==0) throw new IOException("Invalid name for metabolite with id "+id);
                    String compartment = elem.getAttribute(ATTRIBUTE_COMPARTMENT);
                    if(compartment==null || compartment.length()==0) throw new IOException("Invalid compartment for metabolite with id "+id);
                    String formula = elem.getAttribute(ATTRIBUTE_FBCFORMULA);
                    Metabolite metabolite = new Metabolite(id, name, compartment);
                    if(formula!=null) metabolite.setChemicalFormula(formula);
                    network.addMetabolite(metabolite);
                }
            }
        }

    }

    protected void loadReactions(Element listElem, MetabolicNetwork network) throws IOException {
        NodeList offspring = listElem.getChildNodes();
        int k=0;
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_REACTION.equals(elem.getNodeName())) {
                    String id = elem.getAttribute(ATTRIBUTE_ID);
                    if(id==null || id.length()==0) throw new IOException("Every reaction should have an id");
                    String name = elem.getAttribute(ATTRIBUTE_NAME);

                    if(name==null || name.length()==0) throw new IOException("Invalid name for reaction with id "+id);
                    String reversibleStr = elem.getAttribute(ATTRIBUTE_REVERSIBLE);
                    List<GeneProduct> enzymes = new ArrayList<GeneProduct>();
                    List<ReactionComponent> reactants = new ArrayList<ReactionComponent>();
                    List<ReactionComponent> products = new ArrayList<ReactionComponent>();
                    String lbCode = elem.getAttribute(ATTRIBUTE_FBC_LOWERBOUND);
                    String ubCode = elem.getAttribute(ATTRIBUTE_FBC_UPPERBOUND);

                    NodeList offspring2 = elem.getChildNodes();
                    for(int j=0;j<offspring2.getLength();j++) {
                        Node node2 = offspring2.item(j);
                        if (node2 instanceof Element){
                            Element elem2 = (Element) node2;

                            if(ELEMENT_GENEASSOC.equals(elem2.getNodeName())) {
                                enzymes = loadEnzymes(id, elem2,network);
                            }
                            if(ELEMENT_LISTREACTANTS.equals(elem2.getNodeName())) {
                                reactants = loadReactionComponents(id, elem2,network);
                            }
                            if(ELEMENT_LISTMETABPRODUCTS.equals(elem2.getNodeName())) {
                                products = loadReactionComponents(id, elem2,network);
                            }
                        }
                    }
                    if(reactants.size()==0) throw new IOException("No reactants found for reaction "+id);
                    if(products.size()==0) {
                        //System.err.println("No products found for reaction "+id);
                        continue;
                    }

                    Reaction r = new Reaction(id, name, reactants, products);
                    if("true".equals(reversibleStr)) r.setReversible(true);
                    r.setEnzymes(enzymes);
                    //TODO: Load bounds
                    if(lbCode!=null && lbCode.contains("cobra_0")) r.setLowerBoundFlux(0);
                    if(ubCode!=null && ubCode.contains("cobra_0")) r.setUpperBoundFlux(0);
                    network.addReaction(r);
                }
            }
        }

    }

    private List<GeneProduct> loadEnzymes(String reactionId, Element listElem, MetabolicNetwork network) throws IOException {
        List<GeneProduct> answer = new ArrayList<GeneProduct>();
        NodeList offspring = listElem.getChildNodes();
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_GENEPRODUCTREF.equals(elem.getNodeName())) {
                    String enzymeId = elem.getAttribute(ELEMENT_GENEPRODUCT);
                    if(enzymeId==null || enzymeId.length()==0) throw new IOException("Invalid enzyme for reaction "+reactionId);
                    GeneProduct enzyme = network.getGeneProduct(enzymeId);
                    if(enzyme==null) throw new IOException("Enzyme "+enzymeId+" not found for reaction "+reactionId);
                    answer.add(enzyme);
                } else {
                    answer.addAll(loadEnzymes(reactionId, elem, network));
                }
            }
        }


        return answer;
    }

    private List<ReactionComponent> loadReactionComponents(String reactionId, Element listElem, MetabolicNetwork network) throws IOException {
        List<ReactionComponent> answer = new ArrayList<ReactionComponent>();
        NodeList offspring = listElem.getChildNodes();
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_METABREF.equals(elem.getNodeName())) {
                    String metabId = elem.getAttribute(ELEMENT_METABOLITE);
                    if(metabId==null || metabId.length()==0) throw new IOException("Invalid metabolite association for reaction "+reactionId);
                    Metabolite m = network.getMetabolite(metabId);
                    if(m==null) throw new IOException("Metabolite "+metabId+" not found for reaction "+reactionId);

                    String stchmStr = elem.getAttribute(ATTRIBUTE_STOICHIOMETRY);
                    if(stchmStr==null || stchmStr.length()==0) throw new IOException("Absent stoichiometry for metabolite "+metabId+" in reaction "+reactionId);
                    double stoichiometry;
                    try {
                        stoichiometry = Double.parseDouble(stchmStr);
                    } catch (NumberFormatException e) {
                        throw new IOException("Invalid stoichiometry "+stchmStr+" for metabolite "+metabId+" in reaction "+reactionId,e);
                    }
                    ReactionComponent component = new ReactionComponent(m, stoichiometry);
                    answer.add(component);
                }
            }
        }
        return answer;


    }


    public MetabolicNetwork loadNetwork (String filename) throws IOException {
        InputStream is = null;
        MetabolicNetwork mn = null;
        Document doc;
        try {
            is = new FileInputStream(filename);
            DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            doc = documentBuilder.parse(new InputSource(is));
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            if(is!=null) {
                try {
                    is.close();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        Element rootElement = doc.getDocumentElement();
        NodeList offspring = rootElement.getChildNodes();
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_MODEL.equals(elem.getNodeName())) {
                    mn = loadModel(elem);
                }
            }
        }
        if(mn != null) {
            return mn;
        }

        throw new IOException("Malformed XML file. The element "+ELEMENT_MODEL+" could not be found");
    }

    private MetabolicNetwork loadModel(Element modelElem) throws IOException {
        MetabolicNetwork answer = new MetabolicNetwork();
        NodeList offspring = modelElem.getChildNodes();
        for(int i=0;i<offspring.getLength();i++){
            Node node = offspring.item(i);
            if (node instanceof Element){
                Element elem = (Element)node;
                if(ELEMENT_LISTPRODUCTS.equals(elem.getNodeName())) {
                    loadGeneProducts (elem, answer);
                } else if(ELEMENT_LISTMETABOLITES.equals(elem.getNodeName())) {
                    loadMetabolites (elem, answer);
                } else if(ELEMENT_LISTREACTIONS.equals(elem.getNodeName())) {
                    loadReactions (elem, answer);
                }
            }
        }
        return answer;
    }

}