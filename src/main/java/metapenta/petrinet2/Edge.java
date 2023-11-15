package metapenta.petrinet2;

public class Edge<Target>{
    private Target target;
    private double weight;
    public Edge(Target target, double weight) {
        this.target = target;
        this.weight = weight;
    }
    public Target getTarget() {
        return target;
    }
    public double getWeight() {
        return weight;
    }
    public void setTarget(Target target){
         this.target = target;
    }
    public void setWeight(int weight){
        this.weight = weight;
    }

}
