package metapenta.petrinet;

import java.util.List;

public interface IPetriNet {
   public List<Place<?>> getSources();
   public List<Place<?>> getSinks();
}
