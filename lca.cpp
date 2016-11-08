#include <vector>
#include "lca.h"

void LCA::createBitStructures(){
  //Compute shifts
  int val = 1;
  for(int i = 1; i <= num_bits+1; i++){
    shifts[i] = val;
    val <<= 1;
  }

  // Compute masks
  val = 0;
  for(int i = num_bits+1; i >= 1; i--){
    val     |= shifts[i];
    masks[i] = val;
  }

  // Compute left-most and right-most 1 bits for each number
  int    prev_i = 0;
  int    n_new  = 1;
  int    index  = 1;
  leftmost[0]   = -1;
  rightmost[0]  = -1;
  while (prev_i <= (1 << num_bits)-1){
    int next;
    if(prev_i + n_new > (1 << num_bits)-1)
      next = (1<<num_bits)-1;
    else
      next = prev_i + n_new;

    for(int j = prev_i + 1; j <= next; j++){
      leftmost[j]  = index;
      rightmost[j] = leftmost[j & (-j)];
    }

    prev_i += n_new;
    n_new <<= 1;
    index++;
  }
}

void LCA::processTree(SuffixTree& tree){
  if (tree.getNodes().size() > max_num_nodes)
    printErrorAndDie("processTree()", "Number of nodes in tree exceeds the number allocated in LCA object");

  // Compute I for each node
  std::vector<Node*>& nodes = tree.getNodes();
  for(int i = nodes.size()-1; i >= 0; i--){
    int idx = i+1;
    int val = rightmost[i+1]; 
    std::vector<Node*> children = nodes[i]->getChildren();
    for(int j = children.size()-1; j >= 0; j--){
      if (rightmost[I[children[j]->getLabel()]] >= val){
	idx = I[children[j]->getLabel()];
	val = rightmost[idx];
      }
    }
    run_heads[idx] = i+1;
    I[i+1]         = idx;
  }

  // Compute Av for each node
  Av[1] = shifts[rightmost[I[1]]];
  for(int i=2; i <= nodes.size(); i++)
    Av[i] = (shifts[rightmost[I[i]]] | Av[nodes[i-1]->getParent()->getLabel()]);
}

Node* LCA::lca(SuffixTree& tree, Node* x, Node* y){
  int id1 = x->getLabel();
  int id2 = y->getLabel();

  // Check if lca is x or y
  if(id2 >= id1 && id2 <= id1 + x->getNumDescendants())
    return x;
  else if(id1 >= id2 && id1 <= id2 + y->getNumDescendants())
    return y;

  std::vector<Node*>& nodes = tree.getNodes();
  int idx, lca;
  if(I[id1] == I[id2])
    lca = I[id1];
  else {
    idx = leftmost[(I[id1] ^ I[id2])];
    if(rightmost[I[id1]] > idx || rightmost[I[id1]] == num_bits)
      lca = I[id1];
    else if(rightmost[I[id2]] > idx || rightmost[I[id2]] == num_bits)
      lca = I[id2];
    else
      lca = ((masks[idx+1]&I[id1]) | shifts[idx]);
  }
  int h     = rightmost[lca];
  int mask  = masks[h];
  int j     = rightmost[mask & Av[id1] & Av[id2]];

  int x_bar;
  int y_bar;
  int l;

  // Determine x_bar
  l = rightmost[Av[id1]];
  if(l == j)
    x_bar = x->getLabel();
  else{
    int k   = leftmost[Av[id1] & (~ masks[j])];
    int num = (I[id1] & masks[k+1]) | shifts[k];
    int w   = run_heads[num];
    x_bar   = nodes[w-1]->getParent()->getLabel();
  }

  // Determine y_bar
  l = rightmost[Av[id2]];
  if (l == j)
    y_bar = y->getLabel();
  else {
    int k   = leftmost[Av[id2] & (~ masks[j])];
    int num = (I[id2] & masks[k+1]) | shifts[k];
    int w   = run_heads[num];
    y_bar   = nodes[w-1]->getParent()->getLabel();
  }

  if (x_bar < y_bar)
    return nodes[x_bar-1];
  else
    return nodes[y_bar-1];  
}

Node* LCA::lca(SuffixTree& tree, int sfx_idx_1, int sfx_idx_2){
  return lca(tree, tree.getSuffix(sfx_idx_1), tree.getSuffix(sfx_idx_2));
}

int LCA::longestPrefix(SuffixTree& tree, int sfx_idx_1, int sfx_idx_2){
  return lca(tree, sfx_idx_1, sfx_idx_2)->getDepth();
}
