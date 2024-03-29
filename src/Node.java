import java.util.*;
/**
 * Project implementation for different uses of sankoff
 * Sequence Bioinformatics Group Project, WS 22/23
 * Node.java
 * Authors: Vincent Spath, Clarissa Auckenthaler
 */

public class Node {
    private Node left;
    private Node right;
    private Node parent;
    private String name;
    private String[] characterStates;
    private ArrayList<double[]> parsimonyScores;

    public Node(String name) {
        this.left = null;
        this.right = null;
        this.parent = null;
        this.name = name;
        this.characterStates = null;
    }

    public void addChild(Node node) {
        if (this.left == null) {
            this.left = node;
            node.parent = this;
        }
        else if (this.right == null) {
            this.right = node;
            node.parent = this;
        }
        else
            System.out.println("Already added left and right child");
    }

    public Node getLeftChild() {
        return left;
    }

    public Node getRightChild() {
        return right;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getCharacterStates(int ind){return this.characterStates[ind];}

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
    public Node getParent() {return parent;}

    public void setParent(Node parent) {this.parent = parent;}


    public List<Node>  traversePreOrder(Node node, List<Node> list) {
        if (node != null) {
            list.add(node);
            traversePreOrder(node.left,list);
            traversePreOrder(node.right,list);
        }
        return list;
    }


}