package metapenta.model;

import metapenta.petrinet.Edge;
import metapenta.petrinet.Transition;
import metapenta.petrinet.PetriNet;
import metapenta.petrinet.Place;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.linear.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;

import java.util.*;

public class FluxBalanceAnalysis {
    private RealMatrix stequiometryMatrix;
    private Map<String, Integer> rowsMetabolites = new HashMap<>();
    private Map<String, Integer> columnReactions = new HashMap<>();
    private ArrayList<double[]> reactionsBounds = new ArrayList();
    private PetriNet petriNet;

    private int growthReactionIndex;
    private static final double LOWER_LIMIT_NO_REVERSIBLE_RXN = 0;
    private static final double UPPER_LIMIT_RXN = 1000;
    private static final double LOWER_LIMIT_REVERSIBLE_RXN = -1000;


    // Mejor que reciba una metabolic network.
    FluxBalanceAnalysis(PetriNet petriNet, String growthReactionID) {
        int actualColum = 0;
        Set<String> reactionsKeySet = petriNet.getTransitions().keySet();
        for (String reaction : reactionsKeySet) {
            Reaction r = (Reaction) petriNet.getTransition(reaction).getObject();
            columnReactions.put(r.getId(), actualColum);

            if (r.getId().equals(growthReactionID)) {
                growthReactionIndex = actualColum;
            }

            actualColum++;
        }

        int actualRow = 0;
        Set<String> metabolitesKeySet = petriNet.getPlaces().keySet();
        for (String metabolite : metabolitesKeySet) {
            rowsMetabolites.put(petriNet.getPlace(metabolite).getID(), actualRow);
            actualRow++;
        }

        this.petriNet = petriNet;
        loadStequiometryMatrix();
    }

    private void loadStequiometryMatrix() {
        double[][] data = new double[petriNet.getPlaces().size()][petriNet.getTransitions().size()];

        Set<String> reactionsKeys = columnReactions.keySet();

        for (String reaction : reactionsKeys) {
            Integer columnReaction = columnReactions.get(reaction);
            Reaction r = (Reaction) petriNet.getTransition(reaction).getObject();

            double[] bounds = getReactionsBounds(r);
            reactionsBounds.add(bounds);

            List<ReactionComponent> reactants = r.getReactants();

            for (ReactionComponent reactant : reactants) {
                String metaboliteId = reactant.getMetabolite().getId();
                Integer rowMetabolite = rowsMetabolites.get(metaboliteId);
                data[rowMetabolite][columnReaction] = -reactant.getStoichiometry();
            }

            List<ReactionComponent> products = r.getProducts();

            for (ReactionComponent product : products) {
                String metaboliteId = product.getMetabolite().getId();
                int rowMetabolite = rowsMetabolites.get(metaboliteId);
                data[rowMetabolite][columnReaction] = product.getStoichiometry();
            }
        }

        stequiometryMatrix = MatrixUtils.createRealMatrix(data);
    }

    private LinearObjectiveFunction getObjectiveFunction() {
        double[] reactionVector = new double[columnReactions.size()];
        reactionVector[growthReactionIndex] = 1;

        return new LinearObjectiveFunction(reactionVector, 0);
    }

    private LinearConstraintSet getConstraints() {
        Collection constraints = new ArrayList();
        for (int i = 0; i < rowsMetabolites.size(); i++) {
            double[] dmdt = stequiometryMatrix.getRow(i);

            LinearConstraint steadyStateConstraint = new LinearConstraint(dmdt, Relationship.EQ, 0);
            Collection boundsConstraints = getLinearConstraintForBounds(reactionsBounds.get(i), i);

            constraints.add(steadyStateConstraint);
            constraints.addAll(boundsConstraints);

        }

        return new LinearConstraintSet(constraints);
    }

    private Collection getLinearConstraintForBounds(double[] bounds, int reactionNumber){
        double[] reactionVector = new double[columnReactions.size()];
        reactionVector[reactionNumber] = 1;

        LinearConstraint lowerBoundConstraint = new LinearConstraint(reactionVector, Relationship.GEQ, bounds[0]);
        LinearConstraint upperBoundConstraint = new LinearConstraint(reactionVector, Relationship.LEQ, bounds[1]);

        Collection boundsConstraints = new ArrayList();
        boundsConstraints.add(lowerBoundConstraint);
        boundsConstraints.add(upperBoundConstraint);


        // Rea x Metabolito , [flujos], optimizar flujo [i - Indice donde esta la funcion de crecimiento] y las condiciones S*v = 0 <--- SS

        return boundsConstraints;
    }

    private double[] getReactionsBounds(Reaction r) {
        if (r.isReversible()) {
            return new double[]{LOWER_LIMIT_REVERSIBLE_RXN, UPPER_LIMIT_RXN};
        }

        return new double[]{LOWER_LIMIT_NO_REVERSIBLE_RXN, UPPER_LIMIT_RXN};
    }

    public void Optimize() {

        LinearObjectiveFunction lof = getObjectiveFunction();
        LinearConstraintSet linearConstraints = getConstraints();


        double[] initialB = {1, 1, 1, 1, 1, 1};

        MultivariateOptimizer optimizer = new SimplexSolver();

        PointValuePair solution = optimizer.optimize(
                new MaxIter(2000),
                lof,
                linearConstraints,
                GoalType.MAXIMIZE);


        double[] optimalB = solution.getPoint();
        System.out.println(Arrays.toString(optimalB));
    }


    public static void main(String[] args) {
        Metabolite a = new Metabolite("a", "a", "e");
        Metabolite b = new Metabolite("b", "b", "e");
        Metabolite c = new Metabolite("c", "c", "e");
        Metabolite d = new Metabolite("d", "d", "e");
        Metabolite e = new Metabolite("e", "e", "e");

        List<Reaction> reactions = new ArrayList<>();


        //R1 Reaction
        ReactionComponent r1a = new ReactionComponent(a,1);
        ReactionComponent r1d = new ReactionComponent(d,1);
        List<ReactionComponent> reactantsR1 = new ArrayList();
        reactantsR1.add(r1a);
        reactantsR1.add(r1d);

        ReactionComponent r1b = new ReactionComponent(b, 1);
        List<ReactionComponent> productsR1 = new ArrayList();
        productsR1.add(r1b);

        Reaction r1 = new Reaction("r1", "r1", reactantsR1, productsR1);

        //R2 Reaction
        ReactionComponent r2b = new ReactionComponent(b,2);
        List<ReactionComponent> reactantsR2 = new ArrayList();
        reactantsR2.add(r2b);

        ReactionComponent r2e = new ReactionComponent(e,2);
        List<ReactionComponent> productsR2 = new ArrayList();
        productsR2.add(r2e);

        Reaction r2 = new Reaction("r2", "r2", reactantsR2, productsR2);

        //R3 Reaction
        ReactionComponent r3e = new ReactionComponent(e,1);
        List<ReactionComponent> reactantsR3 = new ArrayList();
        reactantsR3.add(r3e);

        ReactionComponent r3a = new ReactionComponent(a,1);
        List<ReactionComponent> productsR3 = new ArrayList();
        productsR3.add(r3a);

        Reaction r3 = new Reaction("r3", "r3", reactantsR3, productsR3);

        //R4 Reaction
        ReactionComponent r4b = new ReactionComponent(b,1);
        ReactionComponent r4e = new ReactionComponent(e,1);
        List<ReactionComponent> reactantsR4 = new ArrayList();
        reactantsR4.add(r4b);
        reactantsR4.add(r4e);

        ReactionComponent r4c = new ReactionComponent(c,2);
        List<ReactionComponent> productsR4 = new ArrayList();
        productsR4.add(r4c);

        Reaction r4 = new Reaction("r4", "r4", reactantsR4, productsR4);


        //R5 Reaction
        ReactionComponent r5a = new ReactionComponent(a,1);
        List<ReactionComponent> productsR5 = new ArrayList();
        productsR5.add(r5a);

        Reaction r5 = new Reaction("r5", "r5", new ArrayList<>(), productsR5);
        r5.setReversible(true);


        //R6 Reaction
        ReactionComponent r6d = new ReactionComponent(d,1);
        List<ReactionComponent> productsR6 = new ArrayList();
        productsR6.add(r6d);

        Reaction r6 = new Reaction("r6", "r6", new ArrayList<>(), productsR6);
        r6.setReversible(true);

        //R7 Reaction
        ReactionComponent r7e = new ReactionComponent(e,1);
        List<ReactionComponent> productsR7 = new ArrayList();
        productsR6.add(r7e);

        Reaction r7 = new Reaction("r7", "r7", new ArrayList<>(), productsR7);
        r7.setReversible(true);

        // Create and load PetriNet
        PetriNet pn = new PetriNet();

        pn.addPlace("a", new Place("a", "a", a));
        pn.addPlace("b", new Place("b", "b", b));
        pn.addPlace("c", new Place("c", "a", c));
        pn.addPlace("d", new Place("d", "d", d));
        pn.addPlace("e", new Place("e", "e", e));

        reactions.add(r1);
        reactions.add(r2);
        reactions.add(r3);
        reactions.add(r4);
        reactions.add(r5);
        reactions.add(r6);
        reactions.add(r7);

        for (Reaction r: reactions) {
            Transition t = new Transition(r.getId(), r.getName(), r);

            for (ReactionComponent reactant: r.getReactants()) {
                Metabolite metabolite = reactant.getMetabolite();

                Place<Metabolite> place = pn.getPlace(metabolite.getId());
                Edge<Place> edge = new Edge<>(place, reactant.getStoichiometry());

                t.AddEdgeIn(edge);
            }

            for (ReactionComponent reactant: r.getProducts()) {
                Metabolite metabolite = reactant.getMetabolite();

                Place<Metabolite> place = pn.getPlace(metabolite.getId());
                Edge<Place> edge = new Edge<>(place, reactant.getStoichiometry());

                t.AddEdgeOut(edge);
            }

            pn.AddTransition(t.getID(), t);
        }

        FluxBalanceAnalysis fba = new FluxBalanceAnalysis(pn, "r4");
        fba.Optimize();

    }

}