package metapenta.petrinet2;

import java.util.ArrayList;
import java.util.List;

public class Transition<A> {
    private String ID;
    private String label;
    private A object;
    List<Edge> edgesIn;
    List<Edge> edgesOut;

    public Transition(String id, String label, A object) {
        ID = id;
        this.label = label;
        this.object = object;
        this.edgesIn = new ArrayList<>();
        this.edgesOut = new ArrayList<>();
    }

    public String getID() {
        return ID;
    }

    public String getLabel() {
        return label;
    }

    public A getObject() {
        return object;
    }

    public void setObject(A attributes){
        this.object = attributes;
    }

    public void setID(String ID){
        this.ID = ID;
    }

    public void setLabel(String label){
        this.label = label;
    }

    public void AddEdgeIn(Edge<?> edge){
        this.edgesIn.add(edge);
    }

    public void AddEdgesIn(List<Edge<?>> edges){
        for(Edge edge: edges){
            this.edgesIn.add(edge);
        }
    }
    public void AddEdgeOut(Edge<?> transition){
        this.edgesOut.add(transition);
    }

    public void AddEdgesOut(List<Edge> edgesOut){
        for(Edge edge: edgesOut){
            this.edgesOut.add(edge);
        }
    }

}
