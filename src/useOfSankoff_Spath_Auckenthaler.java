
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
        String newick = "(((Jaguarundi,Puma),Cheetah),Pallas)";
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

        // this start is the root node
        Node start = parseNewick(newick, characters_dic);

        //1. Step make tree
        //NEED TO BE CREATED BY NEWICK
       /* Node X = new Node("X");
        Node Pallas = new Node("Pallas");
        Node Y = new Node("Y");
        Node Cheetah = new Node ("Cheetah");
        Node Z= new Node("Z");
        Node Jaguarundi= new Node("Jaguarundi");
        Node Puma= new Node("Puma");*/



        /*X.addChild(Pallas);
        X.addChild(Y);
        Pallas.setParent(X);
        Pallas.setCharacterStates(characters_dic);
        Y.addChild(Cheetah);
        Y.addChild(Z);
        Y.setParent(X);
        Cheetah.setParent(Y);
        Cheetah.setCharacterStates(characters_dic);
        Z.setParent(Y);
        Z.addChild(Jaguarundi);
        Z.addChild(Puma);
        Puma.setParent(Z);
        Puma.setCharacterStates(characters_dic);
        Jaguarundi.setParent(Z);
        Jaguarundi.setCharacterStates(characters_dic);*/

        //2. Step create top down and bottom up ordered Lists of Nodes
        List<Node> list =new ArrayList<>();
        list= start.traversePreOrder(start,list);
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
        System.out.println("Changes for i-->j & i' -->j'");
        double[] pScores;

        List<Integer> g = new ArrayList<Integer>();
        Node root= list.get(0);
        //character help - iterator, to solve first all vor one character state and then for the other one..
        int st = 0;
        HashSet<HashSet<List<String>>> uniqueValuesC0= new HashSet<>();
        HashSet<HashSet<List<String>>> uniqueValuesC1= new HashSet<>();

        HashSet<List<String>> c0= new HashSet<>();
        HashSet<List<String>> c1= new HashSet<List<String>>();

        HashSet<List<String>>finalSetsC1= new HashSet<>();
        HashSet<List<String>>finalSetsC2= new HashSet<>();


        //change values of the loop, if the dictionary has more than two elements as values in string array
        //Only compare two characteristics with each other NOT MORE!
        //If more create array of index combinations!!!!!!!!!
        for(int i=0; i<root.getParsimonyScores().size(); i++){
            //System.out.println("Character (0== weight; 1 ==dental): "+st);
            boolean first = true;
            pScores = root.getParsimonyScoresAtIndex(i);
            g = indexOfMin(pScores);
            List<String>changePos0 = new ArrayList<>();
            List<String>changePos1 = new ArrayList<>();
            for (Node node:list){
                if(first&& node.getLeftChild() != null && node.getRightChild() != null){
                    if(i==0){c0=topDown(g,node,pScores,st,changePos0,changePos1);}
                    else if (i==1){c1=topDown(g,node,pScores,st,changePos0,changePos1);}
                    first= false;}
                else if (node.getLeftChild() != null && node.getRightChild() != null){
                    pScores= node.getParsimonyScoresAtIndex(st);
                    g = indexOfMin(pScores);
                    if(i==0){c0= topDown(g,node,pScores,st,changePos0,changePos1);}
                    else if (i==1){c1= topDown(g,node,pScores,st,changePos0,changePos1);}}
            }
            //set P() of the compared characteristics
            uniqueValuesC0= new HashSet<>();
            uniqueValuesC1= new HashSet<>();
            uniqueValuesC0.add(c0);
            uniqueValuesC1.add(c1);

            //update Character
            st+=1;
        }
        // get final sets for further processing
        //eg. weight characteristic
        uniqueValuesC0.forEach(unique -> {unique.forEach(uniqueVal ->{finalSetsC1.add(uniqueVal);});});
        //eg Dental characteristic
        uniqueValuesC1.forEach(unique -> {unique.forEach(uniqueVal ->{finalSetsC2.add(uniqueVal);});});
        System.out.println("Characteristic alpha state changes(i-->j): ");
        System.out.println(finalSetsC1);
        System.out.println("Characteristic beta state changes(i'-->j'): ");
        System.out.println(finalSetsC2);

        finalSetsC1.forEach(e->{

        });
    }

    public static String getChange(Node node,boolean low){
        String changes;
        StringBuilder sb = new StringBuilder(node.getName());
        sb.append(" --> ");
        Node x;
        if(low){
            x=node.getlChange();
            sb.append(x.getName());
            changes= (sb.toString());
        }
        else{
           x=node.getrChange();
           sb.append(x.getName());
           changes= (sb.toString());
        }
        return changes;
    }


    public static HashSet<List<String>> topDown(List<Integer> g, Node node, double[] pScores, int st,List<String>changePos0,List<String>changePos1){
        Double inf = Double.POSITIVE_INFINITY;
        double checkScore;
        Node childLeft= node.getLeftChild();
        Node childRight= node.getRightChild();
        double[] S_left= childLeft.getParsimonyScoresAtIndex(st);
        double[] S_right= childRight.getParsimonyScoresAtIndex(st);
        List<Integer> possibles_left = new ArrayList<Integer>();
        List<Integer> possibles_right = new ArrayList<Integer>();
        for (int i=0; i<S_left.length; i++) {
            if(S_left[i]!= inf){
                possibles_left.add(i);}
            if(S_right[i]!= inf){
                possibles_right.add(i);}
        }
        String x="";
        HashSet<List<String>> c= new HashSet<>();
        for(int j = 0; j<g.size(); j++){
            int index= g.get(j);
            checkScore=pScores[index];
            //System.out.println("Checkscore: "+checkScore+" at index "+index);
            for(int l :possibles_left){
                for (int k:possibles_right){
                    if(index!=l || index!=k){
                        if(checkScore== S_left[l]+S_right[k]+1){
                            //System.out.println("CHANGE IN STATES");
                            //System.out.println("Node: "+ node.getName()+" can be derived by: "+childLeft.getName()+ " Score left= "+S_left[l]+" index left: "+l+" "+childRight.getName()+" Score right= "+S_right[k]+" index right: "+k);
                            if (index==0){
                                if(l==1){
                                    //System.out.println("Low to high "+index+" "+node.getName()+" "+l+" "+childLeft.getName());
                                    node.setlChange(childLeft);
                                    x= getChange(node,true);
                                    changePos0.add(x);}
                                else if(k==1){
                                    //System.out.println("Low to high "+index+" "+node.getName()+" "+k+" "+childRight.getName());
                                    node.setrChange(childRight);
                                    x= getChange(node,false);
                                    changePos0.add(x);}}
                            if (index==1){
                                if(l==0){
                                    //System.out.println("high to low "+index+" "+node.getName()+" "+l+" "+childLeft.getName());
                                    node.setlChange(childLeft);
                                    x= getChange(node,true);
                                    changePos1.add(x);}
                                else if(k==0){
                                    //System.out.println("high to low "+index+" "+node.getName()+" "+k+" "+childRight.getName());
                                    node.setrChange(childRight);
                                    x= getChange(node,false);
                                    changePos1.add(x);}
                            }
                        }
                    }
                    else if(l==k){
                        if(checkScore== S_left[l]+S_right[k]) {//System.out.println("NO CHANGE IN STATES");
                            }
                        else if(checkScore== S_left[l]+S_right[k]+2){//System.out.println("NO CHANGE IN STATES");
                            }
                    }
                }
            }
            c.add(changePos0);
            c.add(changePos1);
        }
        return c;
    }


    public static double computeJaccardIndex( SortedSet<String> A, SortedSet<String> B) {
        // from assignment07 Sequence Bioinfromatics
        Set<String> union = new HashSet<>(A);
        Set<String> intersect = new HashSet<>(A);
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
                        if(animalStates[i].equals(cStates[j])){scores[j] = 0;}
                        else{scores[j]= inf;}
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

    //get parsimonial scores for each internal node
    public static void updateScores(List<Node>newlist, double [][]weightMatrix){
        /**
         * Update the scores in bottom up phase of Sankoffs small parsimoies algorithm
         */
        for(Node node:newlist){
                Node childLeft = node.getLeftChild();
                Node childRight = node.getRightChild();
                ArrayList<double[]> PS = new ArrayList<double[]>();
                for (int k = 0; k < childRight.getCharacterStates().length; k++) {
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
        }
    }

    //helper function to get minimum score in array
    public static double minimum(double [] array) {
        double min = array[0];
        for (int i = 0; i < array.length; i++) {
            if (array[i] < min) {
                min = array[i];}
        }
        return min;
    }
    //helper function to get indices of minimum score in array
    public static List<Integer> indexOfMin(double[] array){
        List<Integer> listOfMin = new ArrayList<Integer>();
        int index =0;
        double min = array[index];
        boolean first = true;
        for (int i = 1; i < array.length; i++){
            if (first= true){
                if (array[i] <= min){
                    listOfMin.add(index,index);
                    first= false;}}
            if (array[i] <= min){
                min = array[i];
                index = i;
                listOfMin.add(index,index);}}
        return listOfMin;
    }

    //helper function to get absolute value in array
    public static double[] absolut(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = Math.abs(array[i]);}
        return array;
    }


    public static Node parseNewick(String newick, Dictionary<String, String[]> characters_dic) {
        // Use a stack to keep track of the nodes as we build the tree
        Stack<Node> stack = new Stack<>();

        // Initialize some variables to keep track of the current node and its properties
        Node current = null;
        Node root = null;
        int ascii = 90;

        // Loop through the characters in the Newick string
        for (int i = 0; i < newick.length(); i++) {
            char c = newick.charAt(i);

            // If the character is an opening parenthesis, create a new node and push it onto the stack
            if (c == '(') {
                current = new Node(Character.toString((char) ascii));
                ascii--;
                if (!stack.isEmpty()) {
                    stack.peek().addChild(current);
                }
                stack.push(current);

            }
            // If the character is a comma continue with the loop
            else if (c == ',') {
                continue;
            }
            // If the character is a closing parenthesis pop node from stack
            else if (c == ')') {
                i++;
                stack.pop();
            }
            else if (Character.isLetterOrDigit(c)) { // if any other character is encountered create new node (leaf)
                String name = "";
                while (Character.isLetterOrDigit(newick.charAt(i))) {
                    name += newick.charAt(i);
                    i++;
                }
                i--;
                current = new Node(name);
                current.setCharacterStates(characters_dic);
                stack.peek().addChild(current);
            }
            if (!stack.isEmpty()) {
                root = stack.peek();
            }

        }

        // The last node on the stack should be the root of the tree
        return root;
    }

}
