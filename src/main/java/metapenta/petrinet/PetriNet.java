package metapenta.petrinet;

import metapenta.model.Metabolite;
import metapenta.model.Reaction;
import metapenta.model.ReactionComponent;
import metapenta.tools.io.DescribeNetworkWriter;

import java.util.*;

public class PetriNet implements IPetriNet {
    private Map<String, Place<Metabolite>> places;
    private Map<String, Transition<Reaction>> transitions;

    public PetriNet(){
        this.places = new TreeMap<>();
        this.transitions = new TreeMap<>();
    }
    public Map<String, Transition<Reaction>> getTransitions() {
        return transitions;
    }

    public Transition getTransition(String key){
        return transitions.get(key);
    }

    public Map<String, Place<Metabolite>> getPlaces() {
        return places;
    }

    public void AddPlace(String id, Place<Metabolite> place){
        this.places.put(id, place);
    }

    public void AddTransition(String id, Transition<Reaction> transition){
        this.transitions.put(id, transition);
    }

    public void setTransitions(Map<String, Transition<Reaction>> transitions) {
        this.transitions = transitions;
    }

    @Override
    public List<Place<?>> getSources() {
        return getPlacesByStatus(Place.SOURCE);
    }

    @Override
    public List<Place<?>> getSinks() {
        return getPlacesByStatus(Place.SINK);
    }


    public void addPlace(String key, Place place){
        places.put(key, place);
    }

    public Place getPlace(String key){
        return places.get(key);
    }


    private List<Place<?>> getPlacesByStatus(String status){
        List<Place<?>> sources = new ArrayList<>();
        for (String key : places.keySet()) {
            if (places.get(key).isStatus(status)){
                sources.add(places.get(key));
            }
        }
        return sources;
    }


    public void describeMetabolicNetwork(String prefixOut) throws Exception {
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

        Set<String> keys = transitions.keySet();
        for (String key: keys){
            Transition<Reaction> transition = transitions.get(key);
            Reaction reaction = transition.getObject();
            List<ReactionComponent> reactants = reaction.getReactants();
            for (ReactionComponent reactant: reactants){
                Metabolite reactantMetabolite = reactant.getMetabolite();
                double reactantStoichiometry = - reactant.getStoichiometry();
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
        Set<String> keys = transitions.keySet();

        for(String key: keys) {
            Reaction reaction = transitions.get(key).getObject();
            reactionIds.add(reaction.getId());
        }

        return reactionIds;
    }


    public List<String> getReversibleReactionsIds(){
        List<String> reactionIds = new ArrayList();
        Set<String> keys = transitions.keySet();

        for(String key: keys) {
            Reaction reaction = transitions.get(key).getObject();
            if(reaction.isReversible()){
                reactionIds.add(reaction.getId());
            }
        }

        return reactionIds;
    }

    public List<String> getIrreversibleReactionsIds(){
        List<String> reactionIds = new ArrayList();
        Set<String> keys = transitions.keySet();

        for(String key: keys) {
            Reaction reaction = transitions.get(key).getObject();
            if(!reaction.isReversible()){
                reactionIds.add(reaction.getId());
            }
        }

        return reactionIds;
    }



}
