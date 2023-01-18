
import java.io.IOException;
import java.util.*;
import java.io.*;

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


        // default path to newick string
        String newick = readNewick("pantherinae.newick");
        //sort of dollo-cost Matrix for 3 states
        double [][]weightMatrix = new double[][]{{0, 1, 1}, {1, 0, 1}, {inf, inf, 0}};

        // default data path
        String characters_dic_path = "data_pantherinae.csv";

        // default state path
        String states_path = "states-pantherinae.csv";


        // read arguments if exactly 3 are given, if none is given use default data
        try {
            if(args.length != 3 && args.length != 0) {
                //throw new IllegalArgumentException("Usage: useOfSankoff_Spath_Auckenthaler newick_path, data_path, states_path");
                throw new IOException("Usage: useOfSankoff_Spath_Auckenthaler newick_path data_path states_path");
            }
            if(args.length == 3) {
                newick = readNewick(args[0]);
                characters_dic_path= args[1];
                states_path = args[2];
            }
        } catch (IllegalArgumentException e) {
            System.out.println(e.getMessage());
        }
        Dictionary<String, String[]> characters_dic = read_data(characters_dic_path);
        String[][] states = read_states(states_path);
        System.out.println(states.length);


        System.out.println("NewickTree:"+newick);

        // this start is the root node
        Node start = parseNewick(newick, characters_dic);


        //2. Step create top down and bottom up ordered Lists of Nodes
        List<Node> list =new ArrayList<>();
        list= start.traversePreOrder(start,list);
        List<Node> newlist =new ArrayList<>();
        newlist= reverseOrder(list);


        //new List of parent Nodes
        List<Node> parents = new ArrayList<>();

        //3. Step initialise states
        parents = initScores(newlist,states,characters_dic);
        //Update states of Parent states
        System.out.println();
        updateScores(parents,states,weightMatrix);

        //4. Step implement top-down phase and expand it to get set of state-chages.
       List<Node> internalnodes= new ArrayList<Node>();
        for(Node node:list){
            if(node.getLeftChild()!=null && node.getRightChild()!=null){
                internalnodes.add(node);
            }
        }

        List<List<List<String>>> pSets = new ArrayList<>();
        //generate ChangeSets
        System.out.println();
        System.out.println("Detect change and type of change for each possible way and for each characteristics.");
        pSets = sankoffTopDown((internalnodes));
        System.out.println("Resulting sets from Sankoffs top down phase: ");
        System.out.println(pSets);

        //5: Step generate Lists of state changes for all characteristics and compare each of them
        List<List<String>> SetLowHigh = new ArrayList<>();
        List<List<String>> SetHighLow = new ArrayList<>();
        //filter sets for changeTypes
        pSets.forEach(p-> {
                //System.out.println(p);
            List<String>sLH = new ArrayList<>();
            List<String>sHL = new ArrayList<>();
            p.forEach(q->{q.forEach(e->{
                        if(Character.isLowerCase(e.charAt(0))){sLH.add(e);}
                        else{sHL.add(e);}});});
            SetLowHigh.add(sLH);
            SetHighLow.add(sHL);
        });

        System.out.println();
        System.out.println("Sets for each character changes between states (0 --> 1) : ");
        System.out.println(SetLowHigh);
        System.out.println();
        System.out.println("Sets for each character changes between states (1 --> 0) : ");
        System.out.println(SetHighLow);

        //6. Step: Calculate jaccard index between eacht combination of states
        //Change 0 --> 1:
        System.out.println(); System.out.println();
        System.out.println("Weighted JaccardIndex for change 0 --> 1 in both characteristics alpha & beta");
        for(int i=0; i<SetLowHigh.size()-1; i++){
            for (int j= 1; j<SetLowHigh.size(); j++){
                if(i!=j) {
                    double jc = computeJaccardIndex(SetLowHigh.get(i), SetLowHigh.get(j),i,j);
                   // System.out.println("Character A= " + i + ", Character B= " + j);
                    //System.out.println(SetLowHigh.get(i) + " " + SetLowHigh.get(j) + " " + jc);
                }
            }
        }
        System.out.println(); System.out.println();

        System.out.println("Weighted JaccardIndex for change 1 --> 0 in both characteristics alpha & beta");
        //Change 1 -->0:
        for(int i=0; i<SetHighLow.size()-1; i++){
            for (int j= 1; j<SetHighLow.size(); j++){
                if(i!=j){
                    double jc=  computeJaccardIndex(SetHighLow.get(i),SetHighLow.get(j),i,j);
                    //System.out.println("Character A= "+i +", Character B= "+j);
                    //System.out.println(SetHighLow.get(i)+" & "+SetHighLow.get(j)+" = "+jc);
                }

            }
        }

        System.out.println(); System.out.println();
        System.out.println("Weighted JaccardIndex for change 1 --> 0 for characteristic alpha and 0-->1 characteristic beta");
        //Change 1 -->0:
        for(int i=0; i<SetHighLow.size()-1; i++){
            for (int j= 1; j<SetLowHigh.size(); j++){
                if(i!=j) {
                    double jc = computeJaccardIndex(SetHighLow.get(i), SetLowHigh.get(j),i,j);
                    //System.out.println("Character A= " + i + ", Character B= " + j);
                    //System.out.println(SetHighLow.get(i) + " & " + SetLowHigh.get(j) + " = " + jc);
                }
            }
        }

        System.out.println(); System.out.println();
        System.out.println("Weighted JaccardIndex for change 0--> 1 for characteristic alpha and 1-->0 characteristic beta");
        //Change 1 -->0:
        for(int i=0; i< SetLowHigh.size()-1; i++){
            for (int j= 1; j<SetHighLow.size(); j++){
                if(i!=j) {
                    double jc = computeJaccardIndex(SetLowHigh.get(i), SetHighLow.get(j),i,j);
                    //System.out.println("Character A= " + i + ", Character B= " + j);
                    //System.out.println(SetLowHigh.get(i) + " & " + SetHighLow.get(j) + " = " + jc);
                }
            }
        }

    }


    public static List<List<List<String>>> sankoffTopDown (List<Node>internalnodes){
        List<List<List<String>>> pSets= new ArrayList<>();
        //loop over different characters to get change List for each character
        for(int i=0; i<internalnodes.get(0).getParsimonyScores().size(); i++){
            System.out.println("Characteristic:" + i);
            List<List<String>> pSetsPerC = new ArrayList<>();
            List<String> pChanges = new ArrayList<>();
            //parsimony scores for each internal node for one character in a list
            List<double[]>parsimonyScores = new ArrayList<double[]>();
            int[] lengths = new int[internalnodes.size()];
            int n= 0;
            for(Node node:internalnodes){
                parsimonyScores.add(node.getParsimonyScoresAtIndex(i));
                lengths[n]= node.getParsimonyScoresAtIndex(i).length;
                n++;
            }
            List<Double> pars= new ArrayList<Double>();
            ArrayList<ArrayList<Integer>> paths = init_td(internalnodes, i);
            for(int j=0; j<paths.size(); j++){
                int[] ind=  paths.get(j).stream().mapToInt(q->q).toArray();
                //System.out.println(Arrays.toString(ind));
                pars= helperScoresCombination(ind,parsimonyScores);
                pChanges= detectChanges(ind,pars, internalnodes,i);
                System.out.println(pChanges);
                System.out.println("");
                pSetsPerC.add(pChanges);
            }
            pSets.add(pSetsPerC);
        }

        return pSets;
    }


    public static List<String> detectChanges(int[] indices, List<Double> pars, List<Node> internalnodes,int i){
        int n= 0;
        boolean lowHigh = false;
        int changeInStates= 0;
        List<String> pChanges = new ArrayList<>();
        Double inf = Double.POSITIVE_INFINITY;
        System.out.println(Arrays.toString(indices));
        for(Node node:internalnodes){
            String x="";
            //System.out.println("N: "+n);
            Node left= node.getLeftChild();
            Node right= node.getRightChild();
            double[] sl; double[] sr; double l; double r;
            double checkscore= pars.get(n);
            int index= indices[n];

                if(left.getLeftChild()== null && right.getRightChild()==null) {
                //last internal node, left node and right node are leafe nodes
                    sl= left.getParsimonyScoresAtIndex(i);
                    sr= right.getParsimonyScoresAtIndex(i);
                    l= sl[index];
                    r= sr[index];
                    //System.out.println("Node:"+ node.getName()+" "+ "left: "+left.getName() +" "+l +"; "+right.getName()+" "+r);
                    if(l==inf&&r==inf)break;
                    if(l == inf){
                        //System.out.println("Transition Change");
                        changeInStates=1;
                        for (int p=0; p<sl.length; p++) {
                            if(sl[p]!= inf){
                                l= sl[p];
                                lowHigh= proofTransitionType(index,p);
                            }
                        }
                        x= getChangeStatements(lowHigh,r,l,changeInStates,checkscore,node,left);
                        if(x!=""){
                            pChanges.add(x);
                        }
                    }

                    if(r == inf){
                        //System.out.println("Transition Change");
                        changeInStates=1;
                        for (int p=0; p<sr.length; p++) {
                            if(sr[p]!= inf){
                                r= sr[p];
                                lowHigh= proofTransitionType(index,p);
                            }
                        }
                        x= getChangeStatements(lowHigh,r,l,changeInStates,checkscore,node,right);
                        if(x!=""){
                            pChanges.add(x);
                        }
                    }

            }
            else if(right.getRightChild()==null&& internalnodes.contains(left)){
                //right leaf is a leaf node
                l= pars.get(n+1);
                sr= right.getParsimonyScoresAtIndex(i);
                r= sr[index];
                //System.out.println("Node:"+ node.getName()+" "+ "left: "+left.getName() +" "+l +"; "+right.getName()+" "+r);
                if(l==checkscore&&r !=0)break;
                if(r == inf){
                    //System.out.println("Transition Change");
                    changeInStates=1;
                    for (int p=0; p<sr.length; p++) {
                        changeInStates=1;
                        if(sr[p]!= inf){
                            r= sr[p];
                            lowHigh= proofTransitionType(index,p);
                        }
                    }
                    x= getChangeStatements(lowHigh,r,l,changeInStates,checkscore,node,right);
                    if(x!=""){
                        pChanges.add(x);
                    }
                }

                if(index != indices[n+1]&& r==0){
                    changeInStates=1;
                    if(l+r+changeInStates==checkscore){
                        if(index==0 && indices[n+1]==1){
                            lowHigh= true;
                        }
                        if(index==1&& indices[n+1]==0){
                            lowHigh= false;
                        }
                    }
                    x= getChangeStatements(lowHigh,r,l,changeInStates,checkscore,node,left);
                    if(x!=""){
                        pChanges.add(x);
                    }
                }

            }

            else if(left.getLeftChild()==null&& internalnodes.contains(right)){
                // left leaf is a leaf node
                r= pars.get(n+1);
                sl= left.getParsimonyScoresAtIndex(i);
                l= sl[index];
                //System.out.println("Node:"+ node.getName()+" "+ "left: "+left.getName() +" "+l +"; "+right.getName()+" "+r);
                if(r==checkscore&& l !=0)break;
                if(l == inf ) {
                    //System.out.println("Transition Change");
                    changeInStates = 1;
                    for (int p = 0; p < sl.length; p++) {
                        if (sl[p] != inf) {
                            l = sl[p];
                            lowHigh= proofTransitionType(index,p);
                        }
                    }
                    x= getChangeStatements(lowHigh,r,l,changeInStates,checkscore,node,left);
                    if(x!=""){
                        pChanges.add(x);
                    }
                }

                if(index != indices[n+1]){
                    changeInStates=1;
                    if(l+r+changeInStates==checkscore){
                        if(r==0&& index==1){
                            lowHigh= false;
                        }
                        if(r==1&& index==0){
                            lowHigh= true;
                        }
                    }
                    x= getChangeStatements(lowHigh,r,l,changeInStates,checkscore,node,right);
                    if(x!=""){
                        pChanges.add(x);
                    }
                }
            }
            else if(internalnodes.contains(left)&&internalnodes.contains(right)){
                sl= left.getParsimonyScoresAtIndex(i);
                sr= right.getParsimonyScoresAtIndex(i);
                l= sl[index];
                r= sr[index];
                //System.out.println("Node:"+ node.getName()+" "+ "left: "+left.getName() +" "+l +"; "+right.getName()+" "+r);
                if(l+r!=checkscore){
                    changeInStates=1;

                    if(l+changeInStates==checkscore){
                        if(index==1){
                            lowHigh= false;
                        }
                        if(index==0){
                            lowHigh= true;
                        }
                        if(lowHigh){ x= node.getName().toLowerCase()+" ---> "+right.getName().toUpperCase();}
                        else{ x= node.getName().toUpperCase()+" ---> "+right.getName().toLowerCase();}
                        if(x!=""){
                            pChanges.add(x);
                        }
                    }
                    if(r+changeInStates==checkscore){
                        if(index==1){
                            lowHigh= false;
                        }
                        if(index==0){
                            lowHigh= true;
                        }
                        if(lowHigh){ x= node.getName().toLowerCase()+" ---> "+left.getName().toUpperCase();}
                        else{ x= node.getName().toUpperCase()+" ---> "+left.getName().toLowerCase();}
                        if(x!=""){
                            pChanges.add(x);
                        }
                    }
                }
                }
            n++;
        }
        return pChanges;

    }

    public static boolean proofTransitionType(int index, int p){
        /**
        * returns type of transition
        * int index= index of current Possible leafnode index
        * int p = changed value
        * */
        boolean lowHigh= false;
            if(index==0&&p==1){
                lowHigh= true;}
            else if(index==1&&p==0){
                lowHigh= false;
            }
        return lowHigh;
    }

    public static String getChangeStatements(boolean lowHigh, double r, double l, int changeInStates, double checkscore, Node node, Node change){
        /**
         * returns change statement for a case
         * r= value of right child
         * l= value of left child
         * changeInStates  = 1 if a transitionstate is detected
         * checkscore = score of the node at certain position
         * node= parent node
         * change = Node where change is detected
         */
        String x= "";
        //No Change statement is not required only for testing purposes
        //if(checkscore==(l+r)){}
        if(checkscore==(l+r)+changeInStates){
            if(lowHigh){ x= node.getName().toLowerCase()+" ---> "+change.getName().toUpperCase();}
            else{ x= node.getName().toUpperCase()+" ---> "+change.getName().toLowerCase();}
        }
        return x;
    }

    public static List<Double> helperScoresCombination(int[] indices,  List<double[]>parsimonyScores){
        List<Double> scores= new ArrayList<Double>();
        for(int n=0; n<indices.length; n++){
            int i = indices[n];
            double[]score = parsimonyScores.get(n);
            double s = score[i];
            scores.add(s);
        }
        return scores;
    }


    public static double computeJaccardIndex( List<String> A, List<String> B, int i , int j) {
        // from assignment07 Sequence Bioinfromatics
        System.out.println("Character A= " + i + ", Character B= " + j);
        List<String> union= new ArrayList<String>(A);
        Set<String> intersect = new HashSet<>(A);
        for(String s :B){if(!union.contains(s)){union.add(s);}}
        intersect.retainAll(B);
        intersect.retainAll(union);
        System.out.println("Union: "+ union);
        System.out.println("Intersect: "+ intersect);

        double unionSize = union.size();
        double intersectionSize = intersect.size();
        double jc = intersectionSize/unionSize;
        if(jc!= 0.0){
            System.out.println("Weighted Jaccard Index: "+ jc);
        }
        else{
            System.out.println("No Correlation");
        }
        System.out.println();
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
    public static void updateScores(List<Node>newlist,String[][] states, double [][]weightMatrix){
        /**
         * Update the scores in bottom up phase of Sankoffs small parsimoies algorithm
         */
        for(Node node:newlist){
                Node childLeft = node.getLeftChild();
                Node childRight = node.getRightChild();
                ArrayList<double[]> PS = new ArrayList<double[]>();
                //newlist.get(0).getCharacterStates().length
                for (int k = 0; k < states.length; k++) {
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

    /**
     * Function to parse a newick String and assigning the characters to the leafes
     * @param newick the newick String
     * @param characters_dic the characters associated to each leaf
     * @return the root node of the tree
     */
    public static Node parseNewick(String newick, Dictionary<String, String[]> characters_dic) {
        // Use a stack to keep track of the nodes as we build the tree
        Stack<Node> stack = new Stack<>();

        // Initialize some variables to keep track of the current node and its properties
        Node current = null;
        Node root = null;
        int ascii = 65;

        // Loop through the characters in the Newick string
        for (int i = 0; i < newick.length(); i++) {
            char c = newick.charAt(i);

            // If the character is an opening parenthesis, create a new node and push it onto the stack
            if (c == '(') {
                current = new Node(Character.toString((char) ascii));
                ascii++;
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

    /**
     * reads the file with the newick String in it
     * @param csvFile The path to the file
     * @return newick String
     */
    public static String readNewick(String csvFile) {
        String line = "";
        try {
            File file = new File(csvFile);
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            // string is in the first line
            line = br.readLine();
            line.replaceAll("\\s+","");
            System.out.println(line);
            br.close();
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }

        return line;
    }

    /**
     * reads in the character states associated with each taxon
     * format: Taxa,character1,character2
     * @param csvFile The path to the file
     * @return a dictionary of the characters, the key being the name of the taxon
     */
    public static Dictionary<String, String[]> read_data(String csvFile) {
        Dictionary<String, String[]> characters_dic = new Hashtable<String, String[]>();
        try {
            File file = new File(csvFile);
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = "";
            String[] tempArr;
            while((line = br.readLine()) != null) {
                // split the lines by comma
                tempArr = line.split(",");
                int end = tempArr.length;
                // the first column is the taxa, the rest are the characters
                characters_dic.put(tempArr[0], Arrays.copyOfRange(tempArr, 1, end));
            }
            br.close();
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }
        return characters_dic;
    }

    /**
     * reads the possible states from a file
     * format: each line: low_state,high_state,unknown
     * same order as in character state data
     * @param csvFile path to the file
     * @return matrix of existing states
     */
    public static String[][] read_states(String csvFile) {
        ArrayList<String[]> states = new ArrayList<>();
        try {
            File file = new File(csvFile);
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = "";
            String[] tempArr;
            while((line = br.readLine()) != null) {
                // split the lines by comma
                tempArr = line.split(",");
                // add to matrix
                states.add(tempArr);
            }
            br.close();
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }
        String[][] r_states = new String[states.size()][];
        for (int i = 0; i < states.size(); i++) {
            r_states[i] = states.get(i);
        }
        return r_states;
    }

    /**
     * initializes the matrix for the top down face of the sankoff algorithm
     * @param internal_nodes list of the internal nodes, starting from the root with left children first
     * @param character the index of the character in the characters_dic
     * @return a matrix of all possible paths in order of the internal nodes
     */
    public static ArrayList<ArrayList<Integer>> init_td (List<Node>internal_nodes, int character){
        ArrayList<ArrayList<Integer>> paths = new ArrayList<>();
        // get the smallest parsimony scores ot the root node
        double smallest = getMin(internal_nodes.get(0).getParsimonyScoresAtIndex(character));
        // new lists with smallest scores in the root node (indices are states)
        for (int i = 0; i < internal_nodes.get(0).getParsimonyScoresAtIndex(character).length; i++) {
            if (internal_nodes.get(0).getParsimonyScoresAtIndex(character)[i] == smallest) {
                ArrayList<Integer> new_list= new ArrayList<>();
                new_list.add(i);
                paths.add(new_list);
            }
        }
        return top_down(internal_nodes,character, paths);
    }

    /**
     * computes all paths of the top down Sankoff for the internal nodes
     * @param internal_nodes list of the internal nodes, starting from the root with left children first
     * @param character the index of the character in the characters_dic
     * @param paths the initialized matrix of all possible paths in order of the internal nodes
     * @return a matrix of all possible paths in order of the internal nodes
     */
    public static ArrayList<ArrayList<Integer>> top_down (List<Node>internal_nodes, int character, ArrayList<ArrayList<Integer>> paths) {
        // iterate over all internal nodes starting with the second one (first used for initialization
        for (int i = 1; i < internal_nodes.size(); i++) {
            Node w_node = internal_nodes.get(i);
            double[] scores_new1;
            double[] scores_new2;
            // assign the scores of the current node to scores_new1 and its neighbor to the other one
            if (w_node == w_node.getParent().getLeftChild()) {
                scores_new1 = w_node.getParent().getLeftChild().getParsimonyScoresAtIndex(character);
                scores_new2 = w_node.getParent().getRightChild().getParsimonyScoresAtIndex(character);
            }
            else {
                scores_new2 = w_node.getParent().getLeftChild().getParsimonyScoresAtIndex(character);
                scores_new1 = w_node.getParent().getRightChild().getParsimonyScoresAtIndex(character);
            }
            // assign the scores of the parent
            double[] scores_old = w_node.getParent().getParsimonyScoresAtIndex(character);
            // iterate all states to get all possible combinations
            // iterate over all possible states (always 3) of the parent
            for (Integer k = 0; k < 2; k++) {
                // iterate over all possible states of the child
                for (int state_new1 = 0; state_new1 <= 2; state_new1++) {
                    // iterate over all possible states of the other child
                    for (int state_new2 = 0; state_new2 <= 2; state_new2++) {
                        // check if combination is possible
                        if (isPath(k,state_new1,state_new2,scores_old[k],scores_new1[state_new1],scores_new2[state_new2])) {
                            ArrayList<Integer> duplicate = new ArrayList<>();
                            boolean cccr = false;
                            Integer state_save = 0;
                            // iterate existing paths
                            for (ArrayList<Integer> path : paths) {
                                // if the parent state exists in the current path
                                if (k.equals(path.get(i - 1))) {
                                    // add new state to path if nothing was assigned yet
                                    if (path.size() == i) {
                                        path.add(state_new1);
                                    }
                                    // save the path to duplicate it later and add the other possible state
                                    else {
                                        duplicate = path;
                                        state_save = state_new1;
                                        cccr = true;
                                    }
                                }
                            }
                            // duplicate the saved path and replace corresponding state
                            if (cccr) {
                                ArrayList<Integer> duplicater = new ArrayList<>(duplicate);
                                duplicater.set(i,state_save);
                                paths.add(duplicater);
                            }
                        }
                    }
                }
            }
        }
        return paths;
    }

    /**
     * computes if one step in the path (top down) is possible
     * @param old_state state of the parent
     * @param new_state1 state of child1
     * @param new_state2 state of child2
     * @param old_score score of the parent
     * @param new_score1 score of child1
     * @param new_score2 score of child2
     * @return boolean that is true if step is possible
     */
    public static boolean isPath (int old_state, int new_state1, int new_state2, Double old_score, Double new_score1, Double new_score2) {
        if (old_state == new_state1 && old_state == new_state2 && new_score1 + new_score2 == old_score) return true;
        if ((old_state == new_state1) && (old_score-1 == new_score2 + new_score1)) return true;
        if ((old_state == new_state2) && (old_score-1 == new_score1 + new_score2)) return true;
        else  return false;
    }

    /**
     * computes the smallest value in an array
     * @param array primitive double array
     * @return smallest value in array
     */
    public static double getMin (double[] array) {
        double min = Double.MAX_VALUE;
        int minIndex = -1;

        for(int i = 0; i < array.length; i++) {
            if (array[i] < min) {
                min = array[i];
                minIndex = i;
            }
        }
        return min;
    }

}
