import java.io.IOException;
import java.util.*;

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
        Double inf = Double.POSITIVE_INFINITY;


        //set values
        String newick = "((Jaguarundi,Puma), Cheetah), Pallas))";
        newick= newick.replace(" ","");
        //sort of dollo-cost Matrix for 3 states
        double [][]weightMatrix = new double[][]{{0, 1, 1}, {1, 0, 1}, {inf, inf, 0}};

        //create dictionary of animals and their states
        Dictionary<String, String[]> characters_dic = new Hashtable<String, String[]>();
        String[] puma= new String[] {"high","30"};
        String[] jaguarundi= new String[] {"low","30"};
        String[] cheetah= new String[] {"high","30"};
        String[] pallas= new String[] {"low","28"};

        characters_dic.put("Puma", puma);
        characters_dic.put("Jaguarundi", jaguarundi);
        characters_dic.put("Cheetah", cheetah);
        characters_dic.put("Pallas", pallas);

        String[][] states= new String[][]{{"low","high","unknown"},{"28","30","unknown"}};

        //read in from arguments
        if (args.length == 3) {
            newick = String.valueOf(args[0]);


        }
        else {
          //throw new IOException("Usage: useOfSankoff_Spath_Auckenthaler newick, matrix");
        }

        System.out.println("NewickTree:"+newick);

        //1. Step make tree
        //NEED TO BE CREATED BY NEWICK
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
        Pallas.setCharacterStates(characters_dic);
        Y.setName("Y");
        Y.addChild(Cheetah);
        Y.addChild(Z);
        Y.setParent(X);
        Cheetah.setName("Cheetah");
        Cheetah.setParent(Y);
        Cheetah.setCharacterStates(characters_dic);
        Z.setName("Z");
        Z.setParent(Y);
        Z.addChild(Jaguarundi);
        Z.addChild(Puma);
        Puma.setName("Puma");
        Puma.setParent(Z);
        Puma.setCharacterStates(characters_dic);
        Jaguarundi.setName("Jaguarundi");
        Jaguarundi.setParent(Z);
        Jaguarundi.setCharacterStates(characters_dic);

        //2. Step create top down and bottom up orderd Lists of Nodes
        List<Node> list =new ArrayList<>();
        list= X.traversePreOrder(X,list);
        List<Node> newlist =new ArrayList<>();
        newlist= reverseOrder(list);


        //new List of parent Nodes
        List<Node> parents = new ArrayList<>();

        //3. Step initalise states
        parents = initScores(newlist,states,characters_dic);
        //Update states of Parent states
        updateScores(parents,weightMatrix);

        //4. Step implement top-down phase and expand it to get set of state-chages.
        //combinations of each possible choice of state for each internal node as elements of a Cartesian product
        System.out.println("Changes for i-->j states eg. low--> high");
        double[] p =X.getParsimonyScoresAtIndex(0);

        double[] pScores;
        double checkScore;


        List<Integer> g = new ArrayList<Integer>();
        Node root= list.get(0);
        int st = 0;

        for(int i=0; i<root.getParsimonyScores().size(); i++){
            System.out.println("State: "+st);
            boolean first = true;
            pScores = root.getParsimonyScoresAtIndex(i);
            g = indexOfMin(pScores);
            for (Node node:list){
                if(first&& node.getLeftChild() != null && node.getRightChild() != null){
                    topDown(g,node,pScores,st);
                    first= false;
                }
                else if (node.getLeftChild() != null && node.getRightChild() != null){
                    pScores= node.getParsimonyScoresAtIndex(st);
                    g = indexOfMin(pScores);
                    topDown(g,node,pScores,st);
                }
            }
            //Get Sets
            for (Node node:list){
                if(node.getLowchange().isEmpty()==false){
                    System.out.println("Nodechange low--> high "+(node.getLowchange()));
                    printChanges(node.getLowchange(),node);

                }
                if(node.getHighchange().isEmpty()==false){
                    System.out.println("Nodechange high--> low "+node.getHighchange());
                    printChanges(node.getHighchange(),node);

                }
            }

            //clean up for new state values
            for (Node node:list){
                List<Node>x= new ArrayList<Node>();
                node.setHighchange(x);
                node.setLowchange(new ArrayList<Node>());
            }

            st+=1;
        }


    }


    public static void printChanges(List<Node> nodes,Node parent){
        for(Node node:nodes){
            System.out.print(parent.getName()+" --> ");
            System.out.print(node.getName());
            System.out.println();
        }
    }
    public static void topDown( List<Integer> g, Node node, double[] pScores, int st){
        Double inf = Double.POSITIVE_INFINITY;
        double checkScore;
        Node childLeft= node.getLeftChild();
        System.out.println("Left child: " + childLeft.getName());
        Node childRight= node.getRightChild();
        System.out.println("Right child: " + childRight.getName());
        double[] S_left= childLeft.getParsimonyScoresAtIndex(st);
        System.out.println("S_left:"+Arrays.toString(S_left));
        double[] S_right= childRight.getParsimonyScoresAtIndex(st);
        System.out.println("S_right:"+Arrays.toString(S_right));

        List<Integer> possibles_left = new ArrayList<Integer>();
        List<Integer> possibles_right = new ArrayList<Integer>();
        for (int i=0; i<S_left.length; i++) {
            if(S_left[i]!= inf){
                possibles_left.add(i);
            }
            if(S_right[i]!= inf){
                possibles_right.add(i);
            }
        }

        for(int j = 0; j<g.size(); j++){
            int index= g.get(j);
            checkScore=pScores[index];
            System.out.println("Checkscore: "+checkScore+" at index "+index);
            for(int l :possibles_left){
                for (int k:possibles_right){
                    if(index!=l || index!=k){
                        if(checkScore== S_left[l]+S_right[k]+1){
                            System.out.println("CHANGE IN STATES");
                            System.out.println("Node: "+ node.getName()+" can be derived by: "+childLeft.getName()+ " Score left= "+S_left[l]+" index left: "+l+" "+childRight.getName()+" Score right= "+S_right[k]+" index right: "+k);
                            if (index==0){
                                if(l==1){
                                    System.out.println("Low to high "+index+" "+node.getName()+l+" "+childLeft.getName());
                                    node.addLowchange(childLeft);
                                }
                                else if(k==1){
                                    System.out.println("Low to high "+index+" "+node.getName()+" "+k+" "+childRight.getName());
                                    node.addLowchange(childRight);
                                }
                            }
                            if (index==1){
                                if(l==0){
                                    System.out.println("high to low "+index+" "+node.getName()+l+" "+childLeft.getName());
                                    node.addHighchange(childLeft);
                                }
                                else if(k==0){
                                    System.out.println("high to low "+index+" "+node.getName()+" "+k+" "+childRight.getName());
                                    node.addHighchange(childRight);
                                }
                            }
                        }
                    }
                    else if(l==k){

                        if(checkScore== S_left[l]+S_right[k]) {
                            System.out.println("NO CHANGE IN STATES");
                            System.out.println("Node: " + node.getName() + " can be derived by: " + childLeft.getName() + " Score left= " + S_left[l] + " index left: " + l + " " + childRight.getName() + " Score right= " + S_right[k] + " index right: " + k);
                        }
                        else if(checkScore== S_left[l]+S_right[k]+2){
                            System.out.println("NO CHANGE IN STATES");
                            System.out.println("Node: " + node.getName() + " can be derived by: " + childLeft.getName() + " Score left= " + S_left[l] + " index left: " + l + " " + childRight.getName() + " Score right= " + S_right[k] + " index right: " + k);
                        }
                    }
                    else{
                        if(checkScore== S_left[l]+S_right[k]+1){
                            System.out.println("CHANGE IN STATES");
                            System.out.println("Node: "+ node.getName()+" can be derived by: "+childLeft.getName()+ " Score left= "+S_left[l]+" index left: "+l+" "+childRight.getName()+" Score right= "+S_right[k]+" index right: "+k);
                        }
                    }
                }
            }
            System.out.println();

        }
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

    public static List<Node> reverseOrder(List<Node>list){
        List<Node> newList = new ArrayList<>();
        for (int i = list.size()-1; i >= 0; i--) {
            newList.add(list.get(i));
        }
        return newList;
    }

    public static List<Node> initScores(List<Node> newlist, String[][] states, Dictionary dic) {
        Double inf = Double.POSITIVE_INFINITY;
        List<Node> parentsList = new ArrayList<>();
        for(Node node:newlist){
            String[] animalStates = node.getCharacterStates();
            ArrayList<double[]> sc = new ArrayList<double[]>();
            //init states
            if(node.getCharacterStates()!= null){
                for (int i=0; i<animalStates.length; i++){
                    String [] cStates= states[i];
                    double [] scores = new double[cStates.length];
                    for(int j=0; j<cStates.length; j++){
                        if(animalStates[i].equals(cStates[j])){
                            scores[j] = 0;
                        }
                        else{
                            scores[j]= inf;
                        }
                    }
                    sc.add(scores);

                    System.out.println("Initialisation-Scores " + node.getName()+" Characterstate "+(i+1)+": "+ Arrays.toString(scores));
                }
                node.setParsimonyScores(sc);
            }
            else {
                parentsList.add(node);
            }
        }
        return parentsList;
    }

    public static void updateScores(List<Node>newlist, double [][]weightMatrix){
        /**
         * Update the scores, using Sankoffs small parsimoies algorithm
         */
        for(Node node:newlist){
            //if (node.getParent()!= null) {
                //Node parent = node.getParent();
                Node childLeft = node.getLeftChild();
                Node childRight = node.getRightChild();
                ArrayList<double[]> PS = new ArrayList<double[]>();

                for (int k = 0; k < childLeft.getCharacterStates().length; k++) {
                    double[] sL = childLeft.getParsimonyScoresAtIndex(k);
                    double[] sR = childRight.getParsimonyScoresAtIndex(k);
                    ArrayList<Double> low_vectorLeft = new ArrayList<Double>();
                    ArrayList<Double> low_vectorRight = new ArrayList<Double>();
                    ArrayList<Double> high_vectorLeft = new ArrayList<Double>();
                    ArrayList<Double> high_vectorRight = new ArrayList<Double>();
                    ArrayList<Double> unknown_vectorLeft = new ArrayList<Double>();
                    ArrayList<Double> unknown_vectorRight = new ArrayList<Double>();
                    low_vectorLeft.add(sL[0]+weightMatrix[0][0]);
                    low_vectorLeft.add(sL[1]+weightMatrix[0][1]);
                    low_vectorLeft.add(sL[2]+weightMatrix[0][2]);
                    low_vectorRight.add(sR[0]+weightMatrix[0][0]);
                    low_vectorRight.add(sR[1]+weightMatrix[0][1]);
                    low_vectorRight.add(sR[2]+weightMatrix[0][2]);
                    high_vectorLeft.add(sL[0]+weightMatrix[1][0]);
                    high_vectorLeft.add(sL[1]+weightMatrix[1][1]);
                    high_vectorLeft.add(sL[2]+weightMatrix[1][2]);
                    high_vectorRight.add(sR[0]+weightMatrix[1][0]);
                    high_vectorRight.add(sR[1]+weightMatrix[1][1]);
                    high_vectorRight.add(sR[2]+weightMatrix[1][2]);
                    unknown_vectorLeft.add(sL[0]+weightMatrix[2][0]);
                    unknown_vectorLeft.add(sL[1]+weightMatrix[2][1]);
                    unknown_vectorLeft.add(sL[2]+weightMatrix[2][2]);
                    unknown_vectorRight.add(sR[0]+weightMatrix[2][0]);
                    unknown_vectorRight.add(sR[1]+weightMatrix[2][1]);
                    unknown_vectorRight.add(sR[2]+weightMatrix[2][2]);
                    double[] S = new double[3];
                    S[0] = ((Collections.min(low_vectorLeft) + Collections.min(low_vectorRight)));
                    S[1] = ((Collections.min(high_vectorLeft) + Collections.min(high_vectorRight)));
                    S[2] = ((Collections.min(unknown_vectorLeft) + Collections.min(unknown_vectorRight)));
                    PS.add(S);
                    System.out.println("Updated Scores " + node.getName() + " "+"Characterstate "+(k+1) + ": " + Arrays.toString(S));

                }
                node.setParsimonyScores(PS);
           // }
        }
    }
    //helper function to get minimum score in array
    public static double minimum(double [] array) {
        double min = array[0];

        for (int i = 0; i < array.length; i++) {
            if (array[i] < min) {
                min = array[i];
            }
        }
        return min;
    }

    public static List<Integer> indexOfMin(double[] array){
        List<Integer> listOfMin = new ArrayList<Integer>();
        int index =0;
        double min = array[index];
        boolean first = true;

        for (int i = 1; i < array.length; i++){
            if (first= true){
                if (array[i] <= min){
                    listOfMin.add(index,index);
                    first= false;
                }
            }
            if (array[i] <= min){
                min = array[i];
                index = i;
                listOfMin.add(index,index);
            }
        }
        return listOfMin;
    }

    //helper function to get absolute value in array
    public static double[] absolut(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = Math.abs(array[i]);
        }
        return array;
    }

}
