
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


        //set values
        String newick = "(((Jaguarundi,Puma),Cheetah),Pallas)";
        newick= newick.replace(" "," ");
        //sort of dollo-cost Matrix for 3 states
        double [][]weightMatrix = new double[][]{{0, 1, 1}, {1, 0, 1}, {inf, inf, 0}};

        //create dictionary of animals and their states
        Dictionary<String, String[]> characters_dic = read_data("test_data.csv");
        /*Dictionary<String, String[]> characters_dic = new Hashtable<String, String[]>();
        String[] puma= new String[] {"high","30"};
        String[] jaguarundi= new String[] {"low","30"};
        String[] cheetah= new String[] {"high","30"};
        String[] pallas= new String[] {"low","28"};

        characters_dic.put("Puma", puma);
        characters_dic.put("Jaguarundi", jaguarundi);
        characters_dic.put("Cheetah", cheetah);
        characters_dic.put("Pallas", pallas);*/

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

        // calculate the correlation between two different characteristics
        //eg. weight characteristic
        uniqueValuesC0.forEach(unique -> {unique.forEach(uniqueVal ->{finalSetsC1.add(uniqueVal);});});
        //eg. Dental characteristic
        uniqueValuesC1.forEach(unique -> {unique.forEach(uniqueVal ->{finalSetsC2.add(uniqueVal);});});



        List<Double> jcInd = new ArrayList();
        finalSetsC1.forEach(e->{
            int i = 0;
            List<String>A= e;
            int finalI = i;
            finalSetsC2.forEach(x->{
                int j= 0;
                List<String>B=x;
                {
                    double jc= computeJaccardIndex(A,B);
                    jcInd.add(jc);
                    //System.out.println(e+" "+x+" Correlation: "+jc);
                }
            });
        });
       //System.out.println("Jaccard indexs for all kind of changes: "+ jcInd);


       List<Node> internalnodes= new ArrayList<Node>();
        for(Node node:list){
            if(node.getLeftChild()!=null && node.getRightChild()!=null){
                internalnodes.add(node);
            }
        }

        List<List<List<String>>> pSets = new ArrayList<>();
        //generate ChangeSets
        pSets = sankoffTopDown((internalnodes));
        System.out.println(pSets);

        //filter sets for changes i-->j and i'---> j'

    }


    public static List<List<List<String>>> sankoffTopDown (List<Node>internalnodes){
        List<List<List<String>>> pSets= new ArrayList<>();
        //loop over different characters to get change List for each character
        for(int i=0; i<internalnodes.get(0).getParsimonyScores().size(); i++){
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
                if(index != indices[n+1]){
                    changeInStates=1;
                    if(l+r+changeInStates==checkscore){
                        if(l==0&& index==1){
                            lowHigh= false;
                        }
                        if(l==1&& index==0){
                            lowHigh= true;
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
                    //System.out.println("Proof this case");
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
            n++;
        }
        return pChanges;

    }

    public static boolean proofTransitionType(int index, int p){
        /*
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
        /*
        returns change statement for a case
        r= value of right child
        l= value of left child
        changeInStates  = 1 if a transitionstate is detected
        checkscore = score of the node at certain position
        node= parent node
        change = Node where change is detected
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


    public static double computeJaccardIndex( List<String> A, List<String> B) {
        // from assignment07 Sequence Bioinfromatics
        Set<String> union = new HashSet<>(A);
        Set<String> intersect = new HashSet<>(A);
        union.addAll(B);
        intersect.retainAll(B);
        intersect.retainAll(union);
        union.removeAll(B);
        double unionSize = union.size();
        double intersectionSize = intersect.size();
        double jc = intersectionSize/unionSize;
        return jc;
    }

    // OBSOLETE
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

    //OBSOLETE
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
                   // System.out.println("g: "+g.get(j)+" k: "+k+" l: "+l);
                    if(index!=l || index!=k){
                        if(checkScore== S_left[l]+S_right[k]+1){
                           // System.out.println("CHANGE IN STATES");
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
    public static Dictionary<String, String[]> read_data(String csvFile) {
        Dictionary<String, String[]> characters_dic = new Hashtable<String, String[]>();
        try {
            File file = new File(csvFile);
            FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);
            String line = "";
            String[] tempArr;
            while((line = br.readLine()) != null) {
                tempArr = line.split(",");
                int end = tempArr.length;
                characters_dic.put(tempArr[0], Arrays.copyOfRange(tempArr, 1, end));
            }
            br.close();
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }
        return characters_dic;
    }
    public static ArrayList<ArrayList<Integer>> init_td (List<Node>internal_nodes, int character){
        ArrayList<ArrayList<Integer>> paths = new ArrayList<>();
        ArrayList<Integer> states = new ArrayList<>();
        double smallest = getMin(internal_nodes.get(0).getParsimonyScoresAtIndex(character));
        // new lists with smallest scores
        for (int i = 0; i < internal_nodes.get(0).getParsimonyScoresAtIndex(character).length; i++) {
            if (internal_nodes.get(0).getParsimonyScoresAtIndex(character)[i] == smallest) {
                ArrayList<Integer> new_list= new ArrayList<>();
                new_list.add(i);
                paths.add(new_list);
                states.add(i);

            }
        }

        return top_down(internal_nodes,character, paths);

    }

    public static ArrayList<ArrayList<Integer>> top_down (List<Node>internal_nodes, int character, ArrayList<ArrayList<Integer>> paths) {

        for (int i = 1; i < internal_nodes.size(); i++) {
            Node w_node = internal_nodes.get(i);
            double[] scores_new1;
            double[] scores_new2;
            if (w_node == w_node.getParent().getLeftChild()) {
                scores_new1 = w_node.getParent().getLeftChild().getParsimonyScoresAtIndex(character);
                scores_new2 = w_node.getParent().getRightChild().getParsimonyScoresAtIndex(character);
            }
            else {
                scores_new2 = w_node.getParent().getLeftChild().getParsimonyScoresAtIndex(character);
                scores_new1 = w_node.getParent().getRightChild().getParsimonyScoresAtIndex(character);
            }
            double[] scores_old = w_node.getParent().getParsimonyScoresAtIndex(character);

            for (Integer k = 0; k < 2; k++) {
                for (int state_new1 = 0; state_new1 <= 2; state_new1++) {
                    for (int state_new2 = 0; state_new2 <= 2; state_new2++) {

                        if (isPath(k,state_new1,state_new2,scores_old[k],scores_new1[state_new1],scores_new2[state_new2])) {


                            ArrayList<Integer> duplicate = new ArrayList<>();
                            boolean cccr = false;
                            Integer state_save = 0;
                            // iterate paths
                            for (ArrayList<Integer> path : paths) {
                                if (k.equals(path.get(i - 1))) {
                                    if (path.size() == i) {
                                        path.add(state_new1);
                                    }
                                    else {
                                        duplicate = path;
                                        state_save = state_new1;
                                        cccr = true;
                                    }
                                }
                            }
                            if (cccr) {
                                ArrayList<Integer> duplicater = new ArrayList<>(duplicate);
                                duplicater.set(i,state_save);
                                paths.add(duplicater);
                            }

                        }
                    }}

            }
        }
        return paths;
    }

    public static boolean isPath (int old_state, int new_state1, int new_state2, Double old_score, Double new_score1, Double new_score2) {
        if (old_state == new_state1 && old_state == new_state2 && new_score1 + new_score2 == old_score) return true;
        if ((old_state == new_state1) && (old_score-1 == new_score2 + new_score1)) return true;
        if ((old_state == new_state2) && (old_score-1 == new_score1 + new_score2)) return true;
        else  return false;
    }

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
