package metapenta.model;

import java.io.IOException;
import java.util.ArrayList;

public interface IMetaPenta {
     void describeMetabolicNetwork(String outFilePrefix) throws Exception;
     void getSinks(String outFilePrefix) throws Exception;
     ArrayList<Double> fluxVector(String outFilePrefix) throws Exception;
}
