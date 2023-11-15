package metapenta.commands;

import metapenta.model.*;

public class DescribeMetabolicNetwork {
    public static void main(String[] args) throws Exception {
        IMetaPenta network = new MetaPenta(args[0]);
        network.describeMetabolicNetwork(args[1]);
    }
}
