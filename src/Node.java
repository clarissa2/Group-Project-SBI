import java.util.*;

public class Node {
    private Node left;
    private Node right;
    private Node parent;
    private final int nodeID;
    private String name;
    private String characterStateWeight;
    private double[] scores;
    private double S_low;
    private double S_high;
    private double S_unknown;
    private String[] characterStates;
    private ArrayList<double[]> parsimonyScores;

    private double [] S_state1;//eg. [S_low,S_28,S_...]
    private double [] S_state2;//eg. [S_high,S_30,S_...]
    private double [] S_state3;//eg. [S_unknown,S_unknown,S_...]

    public Node(int nodeID) {
        this.nodeID = nodeID;
        this.left = null;
        this.right = null;
        this.parent = null;
        this.name = "";
        this.characterStateWeight = "";
        this.characterStates = null;
    }

    public void addChild(Node node) {
        if (this.left == null) {
            this.left = node;
            node.parent = this;
        } else if (this.right == null) {
            this.right = node;
            node.parent = this;
        } else
            System.out.println("Already added left and right child");
    }

    public Node getLeftChild() {
        return left;
    }

    public void setLeftChild(Node leftChild) {
        this.left = leftChild;
    }

    public Node getRightChild() {
        return right;
    }

    public void setRightChild(Node rightChild) {
        this.right = rightChild;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    public String getCharacterStates(int ind){
        return this.characterStates[ind];
    }

    public void setCharacterStates(Dictionary<String,String[]> dictionary){
        this.characterStates = dictionary.get(this.name);
    }

    public String[] getCharacterStates(){
        return this.characterStates;
    }

    public String getCharacterStatesAtIndex(int index){
        return this.characterStates[index];
    }

    public void setParsimonyScores(ArrayList<double[]> scores){
       this.parsimonyScores= scores;

    }

    public ArrayList<double[]> getParsimonyScores(){
            return this.parsimonyScores;
    }
    public double[] getParsimonyScoresAtIndex(int index){
        return this.parsimonyScores.get(index);
    }

    public void setS_state1AtIndex(double score, int index){
        this.S_state1[index] = score;
    }

    public double[] getS_state1(){
        return this.S_state1;
    }

    public double getS_state1AtIndex(int index){
        return this.S_state1[index];
    }

    public void setS_state2AtIndex(double score, int index){
        this.S_state2[index] = score;
    }

    public double[] getS_state2() {
        return this.S_state2;
    }

    public double getS_state2AtIndex(double score, int index){
        return this.S_state2[index];
    }

    public void setS_state3AtIndex(double score, int index){
        this.S_state3[index] = score;
    }

    public double[] getS_state3() {return this.S_state3;}

    public double getS_state3AtIndex(double score, int index){return this.S_state3[index];}

    //public String getCharacterStateWeight() {return characterStateWeight;}

    //public void setCharacterStateWeight(String characterState) {this.characterStateWeight = characterState;}

    public double[] getScores() {return this.scores;}

    public double getScoresAtIndex(int index) {return this.scores[index];}

    public void setScores(double[] scores) {this.scores = scores;}

    public Node getParent() {return parent;}

    public void setParent(Node parent) {this.parent = parent;}

    public int getNodeID() {return nodeID;}

    public double getS_low() {return S_low;}

    public double getS_high() {return S_high;}

    public double getS_unknown() {return S_unknown;}



    public List<Node>  traversePreOrder(Node node, List<Node> list) {
        if (node != null) {
            list.add(node);
            //System.out.println(node.name);
            traversePreOrder(node.left,list);
            traversePreOrder(node.right,list);
        }
        return list;
    }






    //helper function to get minimum score in array
    public double minimum(double [] array) {
        double min = array[0];

        for (int i = 0; i < array.length; i++) {
            if (array[i] < min) {
                min = array[i];
            }
        }
        return min;
    }
    //helper function to get absolute value in array
    public double[] absolut(double[] array) {
        for (int i = 0; i < array.length; i++) {
            array[i] = Math.abs(array[i]);
        }
        return array;
    }
}