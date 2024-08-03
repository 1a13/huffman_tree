#include "bits.h"
#include "treenode.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "priorityqueue.h"
#include "strlib.h"
#include "testing/SimpleTest.h"

using namespace std;

/**
 * Given a Queue<Bit> containing the compressed message bits and the encoding tree
 * used to encode those bits, decode the bits back to the original message text.
 *
 * You can assume that tree is a well-formed non-empty encoding tree and
 * messageBits queue contains a valid sequence of encoded bits.
 *
 * Your implementation may change the messageBits queue however you like. There
 * are no requirements about what it should look like after this function
 * returns. The encoding tree should be unchanged.
 *
 * messageBits guides a temporary pointer through the tree, and resets to the top once
 * it finds a leaf and adds the letter to the result string.
 */
string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits) {
    string result;
    EncodingTreeNode* temp = tree;
    while (!messageBits.isEmpty()) {
        Bit check = messageBits.dequeue();
        if (check == 0) {
            temp = temp->zero;
        } else if (check == 1) {
            temp = temp->one;
        }
        // check if at a leaf after every move
        if (temp->isLeaf()) {
            result += temp->getChar();
            temp = tree;
        }
    }
    return result;
}

/**
 * Reconstruct encoding tree from flattened form Queue<Bit> and Queue<char>.
 *
 * You can assume that the queues are well-formed and represent
 * a valid encoding tree.
 *
 * Your implementation may change the queue parameters however you like. There
 * are no requirements about what they should look like after this function
 * returns.
 *
 * This function takes in a queue of bits for the treeShape, and a queue of
 * characters for the leaves. The tree is created recursively, with the base
 * case being the creation of a leaf, and the recursive case the building of
 * the interior nodes.
 */
EncodingTreeNode* unflattenTree(Queue<Bit>& treeShape, Queue<char>& treeLeaves) {
    Bit build = treeShape.dequeue();
    if (build == 0) {
        EncodingTreeNode* leaf = new EncodingTreeNode(treeLeaves.dequeue());
        return leaf;
    } else {
        EncodingTreeNode* tree = new EncodingTreeNode(unflattenTree(treeShape, treeLeaves),
                                                      unflattenTree(treeShape, treeLeaves));
        return tree;
    }
}

/**
 * Decompress the given EncodedData and return the original text.
 *
 * You can assume the input data is well-formed and was created by a correct
 * implementation of compress.
 *
 * Your implementation may change the data parameter however you like. There
 * are no requirements about what it should look like after this function
 * returns.
 *
 * This function uses the tree shape and tree leaves data to build a tree
 * using the unflattenTree function and then inputs the result into the
 * decodeText function along with the message bits in data to return the
 * encoded message. The tree that's created is also deallocated from memory.
 */
string decompress(EncodedData& data) {
    EncodingTreeNode* tree = unflattenTree(data.treeShape, data.treeLeaves);
    string result = decodeText(tree, data.messageBits);
    deallocateTree(tree);
    return result;
}

/**
 * This function takes in a string and counts the number of times
 * a character appears in the string, tracking the occurences using
 * a map.
 */
Map<char, int> priority(string text) {
    // 106A style
    Map<char, int> counter;
    for (char letter : text) {
        if (counter.containsKey(letter)) {
            counter[letter] += 1;
        } else {
            counter[letter] = 1;
        }
    }
    return counter;
}

/**
 * Constructs an optimal Huffman coding tree for the given text, using
 * the algorithm described in lecture.
 *
 * Reports an error if the input text does not contain at least
 * two distinct characters.
 *
 * When assembling larger trees out of smaller ones, make sure to set the first
 * tree dequeued from the queue to be the zero subtree of the new tree and the
 * second tree as the one subtree.
 *
 * This function takes in a string of text, determines the priority of each
 * character in the text by the number of times it appears. Each character is
 * made into its own leaf node and inserted into the priority queue, and are
 * attached together by interior nodes to build the huffman tree by dequeueing
 * and enqueueing from the priority queue. The huffman tree is returned.
 */
EncodingTreeNode* buildHuffmanTree(string text) {
    if (text.size() < 2) {
        error("Not enough characters in the input string");
    }

    // establish priority values
    PriorityQueue<EncodingTreeNode*> order;
    Map<char, int> counter = priority(text);
    for (char key : counter) {
        EncodingTreeNode* leaf = new EncodingTreeNode(key);
        order.enqueue(leaf, counter[key]);
    }

    // create Huffman tree
    EncodingTreeNode* zero;
    EncodingTreeNode* one;
    EncodingTreeNode* huffman = nullptr;
    while (order.size() != 1) {
        double priority = order.peekPriority();
        zero = order.dequeue();
        priority += order.peekPriority();
        one = order.dequeue();
        EncodingTreeNode* interior = new EncodingTreeNode(zero, one);
        order.enqueue(interior, priority);
    }
    huffman = order.dequeue();
    return huffman;
}

/**
 * This helper function takes in an tree, tracker, and an empty map and adds
 * the characters on the leaves of the tree as keys in the map with their
 * corresponding values as the Vector of bits that lead to the character.
 */
void translation(EncodingTreeNode* tree, Vector<Bit> soFar, Map<char, Vector<Bit>>& key) {
    if (tree->isLeaf()) {
        key[tree->getChar()] = soFar;
    } else {
        translation(tree->zero, soFar + 0, key);
        translation(tree->one, soFar + 1, key);
    }
}

/**
 * Given a string and an encoding tree, encode the text using the tree
 * and return a Queue<Bit> of the encoded bit sequence.
 *
 * You can assume tree is a valid non-empty encoding tree and contains an
 * encoding for every character in the text.
 *
 * This function takes in a string of text to be encoded and an encoding tree,
 * creates a map with each character in the tree as a value with its corresponding
 * path as its value, and iteratively adds each encoded character of the string
 * to the resulting queue of bits.
 */
Queue<Bit> encodeText(EncodingTreeNode* tree, string text) {
    Vector<Bit> soFar;
    Map<char, Vector<Bit>> key;
    translation(tree, soFar, key);
    Queue<Bit> result;
    for (int i = 0; i < text.length(); i++) {
        char word = text[i];
        Vector<Bit> encoded = key[word];
        for (Bit bit : encoded) {
            result.enqueue(bit);
        }
    }
    return result;
}

/**
 * Flatten the given tree into a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment writeup.
 *
 * You can assume the input queues are empty on entry to this function.
 *
 * You can assume tree is a valid well-formed encoding tree.
 *
 * This function traverses a binary tree and adds a 0 to the treeShape queue
 * every time the node is a leaf and adds 1 if it's an interior node. Whenever
 * a leaf is found, the character is also added to the treeLeaves queue.
 */
void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeShape, Queue<char>& treeLeaves) {
    if (tree->isLeaf()) {
        treeShape.enqueue(0);
        treeLeaves.enqueue(tree->getChar());
    } else {
        treeShape.enqueue(1);
        flattenTree(tree->zero, treeShape, treeLeaves);
        flattenTree(tree->one, treeShape, treeLeaves);
    }
}

/**
 * Compress the input text using Huffman coding, producing as output
 * an EncodedData containing the encoded message and flattened
 * encoding tree used.
 *
 * Reports an error if the message text does not contain at least
 * two distinct characters.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 */
EncodedData compress(string messageText) {
    EncodedData result;

    // build huffman tree and flatten it, store in Encoded Data
    EncodingTreeNode* huffman = buildHuffmanTree(messageText);
    Queue<Bit> treeShape;
    Queue<char> treeLeaves;
    flattenTree(huffman, treeShape, treeLeaves);
    result.treeShape = treeShape;
    result.treeLeaves = treeLeaves;

    // encode message and store in EncodedData
    Queue<Bit> encoded = encodeText(huffman, messageText);
    result.messageBits = encoded;

    deallocateTree(huffman);
    return result;
}

/* * * * * * Testing Helper Functions Below This Point * * * * * */

/* This function creates the example tree and returns the pointer
 * to the root node.
 */
EncodingTreeNode* createExampleTree() {
    /* Example encoding tree used in multiple test cases:
     *                *
     *              /   \
     *             T     *
     *                  / \
     *                 *   E
     *                / \
     *               R   S
     */
    EncodingTreeNode* r = new EncodingTreeNode('R');
    EncodingTreeNode* s = new EncodingTreeNode('S');
    EncodingTreeNode* root3 = new EncodingTreeNode(r, s);

    EncodingTreeNode* e = new EncodingTreeNode('E');
    EncodingTreeNode* root2 = new EncodingTreeNode(root3, e);

    EncodingTreeNode* t = new EncodingTreeNode('T');
    EncodingTreeNode* root1 = new EncodingTreeNode(t, root2);

    return root1;
}

/* This function deallocates all the nodes in the tree
 * from memory.
 */
void deallocateTree(EncodingTreeNode* t) {
    // traverse through tree in post-order, as in lecture 22
    if (t == nullptr) {
        return;
    } else {
        deallocateTree(t->zero);
        deallocateTree(t->one);
        delete t;
    }
}

/* This function returns true if two trees have the same shape
 * and same characters in the leaves.
 */
bool areEqual(EncodingTreeNode* a, EncodingTreeNode* b) {
    // base case, both trees are empty
    if (a == nullptr && b == nullptr) {
        return true;
    }
    // recursive case, compare if interior nodes and leaves are in the same locations
    if (a != nullptr && b != nullptr) {
        if (a->isLeaf() && b->isLeaf()) {
            return (a->getChar() == b->getChar());
        } else if (!a->isLeaf() && !b->isLeaf()) {
            return (areEqual(a->zero, b->zero) && areEqual(a->one, b->one));
        }
    }
    return false;
}

/* * * * * * Test Cases Below This Point * * * * * */

STUDENT_TEST("deallocating example tree memory from heap") {
    EncodingTreeNode* tree = createExampleTree();
    deallocateTree(tree);
}

STUDENT_TEST("areEqual comparing two empty trees") {
    EncodingTreeNode* tree1 = nullptr;
    EncodingTreeNode* tree2 = nullptr;

    EXPECT(areEqual(tree1, tree2));

    deallocateTree(tree1);
    deallocateTree(tree2);
}

STUDENT_TEST("areEqual test with one empty tree and simple tree") {
    EncodingTreeNode* tree1 = nullptr;

    EncodingTreeNode* a = new EncodingTreeNode('A');
    EncodingTreeNode* b = new EncodingTreeNode('B');
    EncodingTreeNode* tree2 = new EncodingTreeNode(a, b);

    EXPECT(!areEqual(tree1, tree2));
    deallocateTree(tree1);
    deallocateTree(tree2);
}

STUDENT_TEST("areEqual comparing example tree to simple tree") {
    EncodingTreeNode* tree1 = createExampleTree();

    EncodingTreeNode* a = new EncodingTreeNode('A');
    EncodingTreeNode* b = new EncodingTreeNode('B');
    EncodingTreeNode* tree2 = new EncodingTreeNode(a, b);

    EXPECT(!areEqual(tree1, tree2));
    deallocateTree(tree1);
    deallocateTree(tree2);
}

STUDENT_TEST("areEqual comparing example trees") {
    EncodingTreeNode* tree1 = createExampleTree();
    EncodingTreeNode* tree2 = createExampleTree();

    EXPECT(areEqual(tree1, tree2));
    deallocateTree(tree1);
    deallocateTree(tree2);
}

STUDENT_TEST("areEqual comparing simple trees with different nodes") {
    EncodingTreeNode* a = new EncodingTreeNode('A');
    EncodingTreeNode* b = new EncodingTreeNode('B');
    EncodingTreeNode* tree1 = new EncodingTreeNode(a, b);

    EncodingTreeNode* c = new EncodingTreeNode('C');
    EncodingTreeNode* d = new EncodingTreeNode('D');
    EncodingTreeNode* tree2 = new EncodingTreeNode(c, d);

    EXPECT(!areEqual(tree1, tree2));
    deallocateTree(tree1);
    deallocateTree(tree2);
}

STUDENT_TEST("decodeText, empty message example") {
    EncodingTreeNode* tree = createExampleTree();
    EXPECT(tree != nullptr);

    Queue<Bit> messageBits = {  };
    EXPECT_EQUAL(decodeText(tree, messageBits), "");

    deallocateTree(tree);
}

STUDENT_TEST("unflattenTree, simple tree") {
    EncodingTreeNode* reference = createExampleTree();
    Queue<Bit>  treeShape  = { 1, 0, 0 };
    Queue<char> treeLeaves = { 'A', 'B' };
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    EXPECT(!areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

STUDENT_TEST("decompress using warmup example tree Q4") {
    EncodedData data = {
        { 1, 1, 0, 1, 0, 0, 0 }, // treeShape
        { 'N', 'M', 'S', 'O' },  // treeLeaves
        { 0, 1, 1, 1, 0, 0, 0, 1, 1 } // messageBits
    };

    EXPECT_EQUAL(decompress(data), "SONS");
}

STUDENT_TEST("decompress using warmup example tree Q5") {
    EncodedData data = {
        { 1, 1, 0, 1, 0, 0, 1, 0, 0 }, // treeShape
        { 'F', 'L', 'E', 'R', 'A' },  // treeLeaves
        { 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0 } // messageBits
    };

    EXPECT_EQUAL(decompress(data), "FEEL");
}

STUDENT_TEST("translation using warmup example Q2") {
    Queue<Bit> treeShape = { 1, 1, 0, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'N', 'M', 'S', 'O' };
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    Vector<Bit> soFar;
    Map<char, Vector<Bit>> key;
    translation(tree, soFar, key);

    cout << key << endl; // compare to key from warmup

    deallocateTree(tree);
}

STUDENT_TEST("encode using warmup example Q2") {
    Queue<Bit> treeShape = { 1, 1, 0, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'N', 'M', 'S', 'O' };
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    string word = "SONS";
    Queue<Bit> encode = encodeText(tree, word);

    EXPECT_EQUAL(decodeText(tree, encode), word); // decode the encoded text

    deallocateTree(tree);
}

STUDENT_TEST("flattenTree using warmup Q4") {
    Queue<Bit> treeShape = { 1, 1, 0, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'N', 'M', 'S', 'O' };
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    treeShape = { };
    treeLeaves = { };
    flattenTree(tree, treeShape, treeLeaves);

    EXPECT_EQUAL(treeShape, { 1, 1, 0, 1, 0, 0, 0 });
    EXPECT_EQUAL(treeLeaves, { 'N', 'M', 'S', 'O' });

    deallocateTree(tree);
}

STUDENT_TEST("buildHuffmanTree using warmup Q6") {
    // tree based on dequeue priority
    Queue<Bit> treeShape = { 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'O', 'K', 'B', 'R', 'P', 'E'};
    EncodingTreeNode* check = unflattenTree(treeShape, treeLeaves);

    string text = "BOOKKEEPER";
    EncodingTreeNode* huffman = buildHuffmanTree(text);

    EXPECT(areEqual(huffman, check));

    deallocateTree(huffman);
    deallocateTree(check);
}

STUDENT_TEST("compress using warmup example tree Q4") {
    Queue<Bit> treeShape = { 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'O', 'N', 'S' };
    Queue<Bit> message = { 1, 0, 0, 0, 1, 1 };

    EncodedData data = compress("SONS");

    EXPECT_EQUAL(data.treeShape, treeShape);
    EXPECT_EQUAL(data.treeLeaves, treeLeaves);
    EXPECT_EQUAL(data.messageBits, message);
}

STUDENT_TEST("compress using warmup example tree Q5") {
    Queue<Bit> treeShape = { 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'L', 'F', 'E' };
    Queue<Bit> message = { 0, 1, 1, 1, 0, 0 };

    EncodedData data = compress("FEEL");

    EXPECT_EQUAL(data.treeShape, treeShape);
    EXPECT_EQUAL(data.treeLeaves, treeLeaves);
    EXPECT_EQUAL(data.messageBits, message);
}

/* * * * * Provided Tests Below This Point * * * * */

PROVIDED_TEST("decodeText, small example encoding tree") {
    EncodingTreeNode* tree = createExampleTree(); // see diagram above
    EXPECT(tree != nullptr);

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(decodeText(tree, messageBits), "E");

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(decodeText(tree, messageBits), "SET");

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1}; // STREETS
    EXPECT_EQUAL(decodeText(tree, messageBits), "STREETS");

    deallocateTree(tree);
}

PROVIDED_TEST("unflattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  treeShape  = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'T', 'R', 'S', 'E' };
    EncodingTreeNode* tree = unflattenTree(treeShape, treeLeaves);

    EXPECT(areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

PROVIDED_TEST("decompress, small example input") {
    EncodedData data = {
        { 1, 0, 1, 1, 0, 0, 0 }, // treeShape
        { 'T', 'R', 'S', 'E' },  // treeLeaves
        { 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 } // messageBits
    };

    EXPECT_EQUAL(decompress(data), "TRESS");
}

PROVIDED_TEST("buildHuffmanTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
    EXPECT(areEqual(tree, reference));

    deallocateTree(reference);
    deallocateTree(tree);
}

PROVIDED_TEST("encodeText, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(encodeText(reference, "E"), messageBits);

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(encodeText(reference, "SET"), messageBits);

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1 }; // STREETS
    EXPECT_EQUAL(encodeText(reference, "STREETS"), messageBits);

    deallocateTree(reference);
}

PROVIDED_TEST("flattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  expectedShape  = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'T', 'R', 'S', 'E' };

    Queue<Bit>  treeShape;
    Queue<char> treeLeaves;
    flattenTree(reference, treeShape, treeLeaves);

    EXPECT_EQUAL(treeShape,  expectedShape);
    EXPECT_EQUAL(treeLeaves, expectedLeaves);

    deallocateTree(reference);
}

PROVIDED_TEST("compress, small example input") {
    EncodedData data = compress("STREETTEST");
    Queue<Bit>  treeShape   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeChars   = { 'T', 'R', 'S', 'E' };
    Queue<Bit>  messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0 };

    EXPECT_EQUAL(data.treeShape, treeShape);
    EXPECT_EQUAL(data.treeLeaves, treeChars);
    EXPECT_EQUAL(data.messageBits, messageBits);
}

PROVIDED_TEST("Test end-to-end compress -> decompress") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "Nana Nana Nana Nana Nana Nana Nana Nana Batman"
        "Research is formalized curiosity. It is poking and prying with a purpose. – Zora Neale Hurston",
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(input, output);
    }
}
