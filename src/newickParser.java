import java.util.*;

public class newickParser {
    public static void main(String[] args) {
        Node test = parseNewick("(((Jaguarundi,Puma),Cheetah),Pallas)");
        System.out.println(test.getName());
    }

    public static Node parseNewick(String newick) {
        // Use a stack to keep track of the nodes as we build the tree
        Stack<Node> stack = new Stack<>();

        // Initialize some variables to keep track of the current node and its properties
        Node current = null;

        // Loop through the characters in the Newick string
        for (int i = 0; i < newick.length(); i++) {
            char c = newick.charAt(i);
            char c1 = newick.charAt(i+1);


            // If the character is an opening parenthesis, create a new node and push it onto the stack
             if (c == '(') {
                current = new Node("test");
                if (!stack.isEmpty()) {
                    stack.peek().addChild(current);
                }
                stack.push(current);

            }
            // If the character is a comma, set the current node's right child to the last node on the stack and push it back onto the stack
             else if (c == ',') {
                 current = new Node("test");
                 stack.peek().addChild(current);
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
                 current = new Node(name);
                 stack.peek().addChild(current);
             }
        }

        // The last node on the stack should be the root of the tree
        Node root = stack.pop();
        return root;
    }
    }

