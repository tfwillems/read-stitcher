#ifndef LCA_H
#define LCA_H

#include "suffix_tree.h"

class LCA {
  int  max_num_nodes;
  int  num_bits;
  int* rightmost;
  int* leftmost;
  int* shifts;
  int* masks;
  int* I;
  int* run_heads;
  int* Av;

 protected:
  void createBitStructures();

 public:
  void processTree(SuffixTree& tree);

  LCA(int max_nnodes){
    max_num_nodes = max_nnodes;
    setNumBits();
    leftmost  = new int [1<<num_bits];
    rightmost = new int [1<<num_bits];
    shifts    = new int [num_bits+2];
    masks     = new int [num_bits+2];
    I         = new int [max_num_nodes+1];
    run_heads = new int [max_num_nodes+1];
    Av        = new int [max_num_nodes+1];
    createBitStructures();
  }

  ~LCA(){
    delete [] leftmost;
    delete [] rightmost;
    delete [] shifts;
    delete [] masks;
    delete [] I;
    delete [] run_heads;
    delete [] Av;
  }

  void setNumBits(){
    num_bits  = 0;
    int val   = max_num_nodes;
    while (val > 0){
      num_bits++;
      val >>= 1;
    }
  }

  Node* lca(SuffixTree& tree, Node* x, Node* y);
  Node* lca(SuffixTree& tree, int sfx_idx_1, int sfx_idx_2);
  int   longestPrefix(SuffixTree& tree, int sfx_idx_1, int sfx_idx_2);
};

#endif
