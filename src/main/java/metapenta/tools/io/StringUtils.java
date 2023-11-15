package metapenta.tools.io;

public class StringUtils {

    private String string;

    public StringUtils SetString(String string){
        this.string = string;
        return this;
    }


    public StringUtils addSingleQuotes(){
        this.string = "'" + string + "'";
        return this;
    }

    public StringUtils addBreakLine(){
        this.string = string + "\n";
        return this;
    }

    public StringUtils addString(String string){
        this.string = this.string + string;
        return this;
    }

    public StringUtils addDouble(double d){
        this.string = this.string + d;
        return this;
    }

    public StringUtils addInt(int i){
        this.string = this.string + i;
        return this;
    }


    public StringUtils addEmptySpace(){
        this.string = string + " ";
        return this;
    }

    public String GetString(){
        return this.string;
    }

}
