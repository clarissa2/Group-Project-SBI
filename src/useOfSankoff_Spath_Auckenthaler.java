import java.io.IOException;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Enumeration;
/**
 * Project implementation for different uses of sankoff
 * Sequence Bioinformatics Group Project, WS 22/23
 * useOfSankoff_Spath_Auckenthaler
 * Authors: Vincent Spath, Clarissa Auckenthaler
 */
public class useOfSankoff_Spath_Auckenthaler {

    /**
     * run and use of sankoffs small parsimony algorithm
     * @param args commandline arguments
     * @throws IOException if arguments are incorrect or sequences not found
     */
    public static void main(String[] args) throws IOException {
        System.out.println("useOfSankoff_Spath_Auckenthaler");
        //initialize
        String[] characterStates= new String[] {"low", "high", "unknown"};
        String newick = "((Jaguarundi,Puma), Cheetah), Pallas))";
        newick= newick.replace(" ","");
        Double inf = Double.POSITIVE_INFINITY;

        double [][]weightMatrix = new double[][]{{0, 1, 1}, {1, 0, 1}, {inf, inf, 0}};

        Dictionary<String, String> character_dic = new Hashtable<String, String>();
        character_dic.put("Puma", "high");
        character_dic.put("Jaguarundi", "low");
        character_dic.put("Cheetah", "high");
        character_dic.put("Pallas", "low");

        String[] DentalStates= new String[] {"28", "30", "unknown"};
        double [][]dentalMatrix = new double[][]{{0, 1, 1}, {1, 0, 1}, {inf, inf, 0}};
        Dictionary<String, String> characterDental_dic = new Hashtable<String, String>();
        characterDental_dic.put("Puma", "30");
        characterDental_dic.put("Jaguarundi", "30");
        characterDental_dic.put("Cheetah", "30");
        characterDental_dic.put("Pallas", "28");


        //read in from argument
        if (args.length == 3) {
            newick = String.valueOf(args[0]);
            //weightMatrix=;
            //characterStates= ;

        }
        //else {
          //  throw new IOException("Usage: useOfSankoff_Spath_Auckenthaler [newick, matrix]");
        //}

        System.out.println("NewickTree:"+newick);
        //ArrayList<Node> nodes_newick = new ArrayList<Node>();

        Node X = new Node(0);
        Node Pallas = new Node(1);
        Node Y = new Node(2);
        Node Cheetah = new Node (3);
        Node Z= new Node(4);
        Node Jaguarundi= new Node(5);
        Node Puma= new Node(6);

        X.setName("X");
        X.addChild(Pallas);
        X.addChild(Y);

        Pallas.setName("Pallas");
        Pallas.setParent(X);
        Pallas.setCharacterStateDetal(characterDental_dic.get(Pallas.getName()));
        Pallas.setCharacterStateWeight(character_dic.get(Pallas.getName()));
        Pallas.initScores(characterStates,character_dic);

        Y.setName("Y");
        Y.addChild(Cheetah);
        Y.addChild(Z);
        Y.setParent(X);

        Cheetah.setName("Cheetah");
        Cheetah.setParent(Y);
        Cheetah.setCharacterStateWeight(character_dic.get(Cheetah.getName()));
        Cheetah.initScores(characterStates,character_dic);


        Z.setName("Z");
        Z.setParent(Y);
        Z.addChild(Jaguarundi);
        Z.addChild(Puma);

        Puma.setName("Puma");
        Puma.setParent(Z);
        Puma.setCharacterStateWeight(character_dic.get(Puma.getName()));
        Puma.initScores(characterStates,character_dic);

        Jaguarundi.setName("Jaguarundi");
        Jaguarundi.setParent(Z);
        Jaguarundi.setCharacterStateWeight(character_dic.get(Jaguarundi.getName()));
        Jaguarundi.initScores(characterStates,character_dic);


        Z.updateScores(Z.getLeftChild(),Z.getRightChild(),weightMatrix);
        System.out.println("New Z-Scores:"+Arrays.toString(Z.getScores()));
        System.out.println("New Puma-Scores:"+Arrays.toString(Puma.getScores()));
        Y.updateScores(Y.getLeftChild(),Y.getRightChild(),weightMatrix);
        System.out.println("New Y-Scores:"+Arrays.toString(Y.getScores()));

        X.updateScores(X.getLeftChild(),X.getRightChild(),weightMatrix);
        System.out.println("New X-Scores:"+Arrays.toString(X.getScores()));

        String[] p1_low_high= new String[]{};
        System.out.println(X.getRightChild().getCharacterStateWeight());
        /**
         *  int nodeID= 0;
         *  boolean first = true;
         *  Node currentNode = null;
         * for (int i=0; i<newick.length() ;i++){
         String c= Character.toString(newick.charAt(i));
         if (first== true){
         root= new Node (0);
         currentNode = root;
         first= false;
         }
         else if(c.equals("(")){
         nodeID++;
         node = new Node(nodeID);
         node.setParent(currentNode);
         currentNode.addChild(node);
         currentNode= node;
         }
         else if (c.equals(",")){
         if ()
         }

         }
         */


    }


    public static double computeJaccardIndex( SortedSet<Integer> A, SortedSet<Integer> B) {
        // from assignment07 Sequence Bioinfromatics
        Set<Integer> union = new HashSet<>(A);
        Set<Integer> intersect = new HashSet<>(A);

        union.addAll(B);
        intersect.retainAll(B);

        intersect.retainAll(union);
        double unions = union.size();
        double intersection = intersect.size();
        double jc = intersection/unions;
        return jc;
    }

}
