#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "suffix_tree.h"

void printErrorAndDie(std::string fname, std::string error) {
  std::cerr << "ERROR in function "<< fname << std::endl;
  std::cerr << "CAUSE: " << error << std::endl << "EXITING..." << std::endl;
  exit(1);
}

/* Beginning of Node-related functions */
void Node::extendEdge(int edge_idx){ 
  *stops[edge_idx] = *stops[edge_idx]+1;
}

void Node::splitEdge(int edge_idx, size_t pos_idx, size_t char_idx, std::vector<int>& ids, size_t*& global_stop, std::string& tokens, Node*& new_node){
  size_t old_char_idx = *starts[edge_idx] + pos_idx;
  Node*  new_child    = new Node(this, edge_idx);
  new_child->addEdge(ids[old_char_idx], old_char_idx, *stops[edge_idx]);
  new_child->addEdge(ids[char_idx], char_idx, char_idx);
  new_child->setGlobalStop(ids[char_idx], global_stop);
  new_node = new_child;

  if(children[edge_idx] != NULL){
    new_child->children[ids[old_char_idx]] = children[edge_idx];
    children[edge_idx]->edge_idx = ids[old_char_idx];
    children[edge_idx]->parent   = new_child;
  }
  children[edge_idx] = new_child;

  if (global[edge_idx]){
    new_child->setGlobalStop(ids[old_char_idx], global_stop);
    global[edge_idx] = false;
    stops[edge_idx]  = new size_t;
  }
  *stops[edge_idx]   = *starts[edge_idx]+pos_idx-1;
}

void Node::print(std::ostream& out, int depth, std::string& word){
  std::string spacing;
  for(int i = 0; i < depth; i++)
    spacing += " ";
  out << spacing << "*" << std::endl;
  out << spacing << label << std::endl;

  for (int i = 0; i < NUM_CHARS; i++){
    if (starts[i] != NULL){
      out << spacing << word.substr(*starts[i], getNumChars(i)) << std::endl;
      if (children[i] != NULL)
	children[i]->print(out, depth+getNumChars(i), word);
    }
  }
}

void Node::dfsProcess(size_t& cur_label, int cur_depth, std::vector< std::pair<int, Node*> >& stack, int max_depth, std::vector<Node*>& suffixes){
  depth = cur_depth;
  setLabel(cur_label);
  cur_label++;

  for(int i = NUM_CHARS-1; i >= 0; i--){
    if(starts[i] != NULL){
      if(children[i] != NULL)
	stack.push_back(std::pair<int, Node*>(cur_depth+getNumChars(i), children[i]));
      else{
	// Add leaf node
	children[i] = new Node(this, i);
	int new_depth = cur_depth + getNumChars(i);
	stack.push_back(std::pair<int, Node*>(new_depth, children[i]));
	suffixes[max_depth-new_depth] = children[i];
      }
    }
  }
}
/* End of Node-related functions */


/* Beginning of SuffixTree-related functions */
void SuffixTree::stringToIds(){
  ids.clear();
  if (tokens.size() == 0)
    printErrorAndDie("stringToIds()", "Empty string cannot be provided for suffix tree construction");

  for(int i = 0; i < tokens.size(); i++){
    int id = charID[tokens[i]];
    if(id == -1)
      printErrorAndDie("stringToIds()", "Invalid character ");
    ids.push_back(id);
  }
}

void SuffixTree::traverseChars(Node* start_node, size_t start_idx, size_t stop_idx, Node*& end_node, int& end_edge_idx, size_t& end_pos_idx){
  while(start_idx <= stop_idx && start_node != NULL && (stop_idx - start_idx + 1 > start_node->getNumChars(ids[start_idx]))){
    Node* new_node = start_node->getChildNode(ids[start_idx]);
    start_idx += start_node->getNumChars(ids[start_idx]);
    start_node = new_node;
  }

  end_node     = start_node;
  end_edge_idx = ids[start_idx];
  end_pos_idx  = stop_idx - start_idx + 1;
}

int SuffixTree::extendCharacter(Node* node, int edge_idx, size_t pos_index, int char_id, size_t char_idx, 
				Node*& ins_node, int& ins_edge_idx, size_t& ins_pos_idx, size_t*& global_stop){
  if (!node->edgeExists(edge_idx)){
    node->addEdge(edge_idx, char_idx, char_idx);
    node->setGlobalStop(edge_idx, global_stop);
    ins_node     = node;
    ins_edge_idx = edge_idx;
    ins_pos_idx  = 0;
    return 2;
  }

  size_t edge_len = node->getNumChars(edge_idx);
  if (pos_index == edge_len){
    if(node->isLeaf(edge_idx)){
      node->extendEdge(edge_idx);
      node->setGlobalStop(edge_idx, global_stop);
      ins_node     = node;
      ins_edge_idx = edge_idx;
      ins_pos_idx  = node->getNumChars(edge_idx);
      return 1;
    }
    else if(!node->getChildNode(edge_idx)->edgeExists(char_id)){
      node->getChildNode(edge_idx)->addEdge(char_id, char_idx, char_idx);
      node->getChildNode(edge_idx)->setGlobalStop(char_id, global_stop);
      ins_node     = node->getChildNode(edge_idx);
      ins_edge_idx = char_id;
      ins_pos_idx  = 0;
      return 2;
    }
    else {
      ins_node     = node->getChildNode(edge_idx);
      ins_edge_idx = char_id;
      ins_pos_idx  = 0;
      return 3;
    }
  }
  else if (pos_index <= edge_len-1){
    if (ids[node->getStart(edge_idx)+pos_index] != char_id){
      node->splitEdge(edge_idx, pos_index, char_idx, ids, global_stop, tokens, ins_node);
      ins_edge_idx = char_id;
      ins_pos_idx  = 0;
      return 2;
    }
    else {
      ins_node     = node;
      ins_edge_idx = edge_idx;
      ins_pos_idx  = pos_index;
      return 3;
    }
  }

  std::cerr << "ERROR: Invalid case reached. Exiting..." << std::endl;
  exit(1);
}

void SuffixTree::getLinkInfo(Node* prev_node, int edge_idx, size_t pos_idx, Node*& node, size_t& gamma_start, size_t& gamma_stop){
  if(pos_idx == 0){
    node        = prev_node->getParent();
    gamma_start = prev_node->getParent()->getStart(prev_node->getEdgeIndex());
    gamma_stop  = prev_node->getParent()->getStop(prev_node->getEdgeIndex());
  }
  else {
    node        = prev_node;
    gamma_start = prev_node->getStart(edge_idx);
    gamma_stop  = prev_node->getStop(edge_idx)-1;
  }
}

int SuffixTree::singleExtension(size_t i, size_t j, size_t*& global_stop, Node* prev_node, int edge_idx, 
				size_t pos_idx, Node*& new_node, int& new_edge_idx, size_t& new_pos_idx){
  Node* parent;
  size_t gamma_start, gamma_stop;
  Node*  trav_parent;
  int    trav_edge;
  size_t trav_pos;
  getLinkInfo(prev_node, edge_idx, pos_idx, parent, gamma_start, gamma_stop);

  if(!parent->isRoot())
    parent = parent->getSuffixLink();
  if(parent->isRoot())
    traverseChars(root, j, i-1, trav_parent, trav_edge, trav_pos);
  else
    traverseChars(parent, gamma_start, gamma_stop, trav_parent, trav_edge, trav_pos);
  return extendCharacter(trav_parent, trav_edge, trav_pos, ids[i], i, new_node, new_edge_idx, new_pos_idx, global_stop);
}

void SuffixTree::createTree(){
  global_end    = new size_t;
  *global_end   = 0;
  size_t j_star = 0;

  root = new Node(NULL, -1);
  root->addEdge(ids[0], 0, 0);
  root->setGlobalStop(ids[0], global_end);

  Node*  prev_node      = root;
  int    prev_edge_idx  = ids[0];
  size_t prev_pos_idx   = 0;
  Node*  new_node;
  int    new_edge_idx;
  size_t new_pos_idx;

  for(size_t i = 1; i < ids.size(); i++){
    *global_end = *global_end + 1;
    int prev_ext_type = extendCharacter(prev_node, prev_edge_idx, prev_pos_idx+1, ids[i], i, prev_node, prev_edge_idx, prev_pos_idx, global_end);
    if(prev_ext_type == 3)
      if(prev_node->isLeaf(prev_edge_idx) && prev_node->getNumChars(prev_edge_idx)-1 == prev_pos_idx)
	prev_ext_type = 1;

    if (prev_ext_type != 3){
      size_t j;
      for(j = j_star+1; j <= i; j++){
	int ext_type = singleExtension(i, j, global_end, prev_node, prev_edge_idx, prev_pos_idx, new_node, new_edge_idx, new_pos_idx);
	if(prev_ext_type == 2)
	  prev_node->setSuffixLink(new_node);

	prev_node     = new_node;
	prev_edge_idx = new_edge_idx;
	prev_pos_idx  = new_pos_idx;
	prev_ext_type = ext_type;

	if(ext_type == 3)
	  break;
      }
      if(j == i+1)
	j_star = j-1;
      else
	j_star = j;
    }
  }
}

void SuffixTree::dfsProcess(){
  nodes.clear();
  size_t label = 1;
  std::vector< std::pair<int,Node*> > stack;
  stack.push_back(std::pair<int,Node*>(0,root));
  suffixes = std::vector<Node*>(tokens.size(), NULL);
  while(stack.size() != 0){
    std::pair<int,Node*> info = stack.back(); stack.pop_back();
    nodes.push_back(info.second);
    info.second->dfsProcess(label, info.first, stack, tokens.size(), suffixes);
  }

  for(int i = nodes.size()-1; i > 0; i--)
    nodes[i]->getParent()->addDescendants(nodes[i]->getNumDescendants()+1);
}

SuffixTree::SuffixTree(std::string myString){
  charID['A'] = 0;
  charID['C'] = 1;
  charID['G'] = 2;
  charID['T'] = 3;
  charID['$'] = 4;
  charID['#'] = 5;
  tokens      = myString+"$";
  stringToIds();
  createTree();
  dfsProcess();
}

SuffixTree::~SuffixTree(){
  if (root != NULL)
    delete root;
  delete global_end;
}

std::string SuffixTree::getPathString(Node* end){
  return (end->isRoot() ? "" : tokens.substr(end->getParent()->getStop(end->getEdgeIndex())-end->getDepth() +1, end->getDepth()));
}

Node* SuffixTree::getSuffix(int idx){
  return suffixes[idx];
}

void SuffixTree::printTree(std::ostream& out){
  root->print(out, 1, tokens);
}
/* End of SuffixTree-related functions */
