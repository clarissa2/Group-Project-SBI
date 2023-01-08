import java.util.Dictionary;

public class Node {
    private Node left;
    private Node right;
    private Node parent;
    private final int nodeID;
    private String name;
    private String characterStateWeight;
    private String characterStateDental;
    private int score;
    private double[] scores;
    private double S_low;
    private double S_high;
    private double S_unknown;

    public Node(int nodeID) {
        this.nodeID = nodeID;
        this.left = null;
        this.right = null;
        this.parent = null;
        this.name = "";
        this.characterStateWeight = "";
        this.characterStateDental ="";
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

    public String getCharacterStateWeight() {
        return characterStateWeight;
    }

    public void setCharacterStateWeight(String characterState) {
        this.characterStateWeight = characterState;
    }

    public String getCharacterStateDental() {
        return characterStateDental;
    }

    public void setCharacterStateDetal(String characterState) {
        this.characterStateDental = characterState;
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public double[] getScores() {
        return this.scores;
    }
    public double getScoresAtIndex(int index) {
        return this.scores[index];
    }

    public void setScores(double[] scores) {
        this.scores = scores;
    }

    public Node getParent() {
        return parent;
    }

    public void setParent(Node parent) {
        this.parent = parent;
    }

    public int getNodeID() {
        return nodeID;
    }

    public double getS_low() {
        return S_low;
    }
    public double getS_high() {
        return S_high;
    }
    public double getS_unknown() {
        return S_unknown;
    }

    public void initScores(String[] characterStates, Dictionary dictionary) {
        /**
         * G
         */
        Double inf = Double.POSITIVE_INFINITY;
        if(dictionary.get(this.name).equals(characterStates[0])){
            double[] l = new double[]{0, inf, inf};
            this.scores= l;
        } else if (dictionary.get(this.name).equals(characterStates[1])) {
            double[] h = new double[]{inf, 0, inf};
            this.scores = h;

        } else {
            double[] u = new double[]{inf, inf, 0};
            this.scores = u;

        }
    }

    public void updateScores(Node childnode1, Node childnode2, double [][]weightMatrix){
        /**
         * Update the scores, using Sankoffs small parsimoies algorithm
         */
        Double inf = Double.POSITIVE_INFINITY;

        double [] low_vectorChild1 =new double[]{(childnode1.getScoresAtIndex(0))+weightMatrix[0][0],childnode1.getScoresAtIndex(1)+weightMatrix[0][1],childnode1.getScoresAtIndex(2)+weightMatrix[0][2]} ;
        double [] low_vectorChild2 =new double[]{(childnode2.getScoresAtIndex(0))+weightMatrix[0][0],childnode2.getScoresAtIndex(1)+weightMatrix[0][1],childnode2.getScoresAtIndex(2)+weightMatrix[0][2]};
        this.S_low= (minimum(absolut(low_vectorChild1))+minimum(absolut(low_vectorChild2)));


        double [] high_vectorChild1 =new double[]{(childnode1.getScoresAtIndex(0))+weightMatrix[1][0],childnode1.getScoresAtIndex(1)+weightMatrix[1][1],childnode1.getScoresAtIndex(2)+weightMatrix[1][2]} ;
        double [] high_vectorChild2 =new double[]{(childnode2.getScoresAtIndex(0))+weightMatrix[1][0],childnode2.getScoresAtIndex(1)+weightMatrix[1][1],childnode2.getScoresAtIndex(2)+weightMatrix[1][2]};
        this.S_high = (minimum(absolut(high_vectorChild1))+minimum(absolut(high_vectorChild2)));

        this.S_unknown = inf;

        this.scores= new double[]{S_low,S_high,S_unknown};
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