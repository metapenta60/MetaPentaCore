package metapenta.petrinet2;

import metapenta.model.Metabolite;

import java.util.ArrayList;
import java.util.List;

public class Place<O> {
    public static final String SINK = "SINK";
    public static final String SOURCE = "SOURCE";
    private String Id;
    private String label;
    private O object;
    List<Edge<?>> edgesIn;
    List<Edge<?>> edgesOut;

    public Place(String id, String label, O object) {
        Id = id;
        this.label = label;
        this.object = object;
        this.edgesIn = new ArrayList<>();
        this.edgesOut = new ArrayList<>();
    }
    public void setLabel(String label) {
        this.label = label;
    }
    public void setId(String ID) {
        this.Id = ID;
    }
    public void setObject(O object){
        this.object = object;
    }

    public O getObject() {
       return object;
    }
    public String getLabel() {
        return label;
    }

    public List<Edge<?>> getEdgesIn() {
        return edgesIn;
    }

    public List<Edge<?>> getEdgesOut() {
        return edgesOut;
    }

    public String getID() {
        return Id;
    }
    public O getAttributes() {
        return object;
    }
    public void AddEdgeIn(Edge transition){
        this.edgesIn.add(transition);
    }
    public void AddEdgeOut(Edge transition){
        this.edgesOut.add(transition);
    }

    public boolean isStatus(String status){
        return switch (status) {
            case SINK -> isSource();
            case SOURCE -> isSink();
            default -> false;
        };
    }

    public boolean isSource(){
        return edgesIn.isEmpty();
    }

    public boolean isSink(){
        return edgesOut.isEmpty();
    }

}
