#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#ifndef SUFFIX_TREE_H
#define SUFFIX_TREE_H
#define NUM_CHARS 6

void printErrorAndDie(std::string fname, std::string error);

class Node {
private:
  Node*    parent;
  size_t*  starts[NUM_CHARS];
  size_t*  stops[NUM_CHARS];
  bool     global[NUM_CHARS];
  Node*    children[NUM_CHARS];
  Node*    sfx_link;
  int      edge_idx;
  size_t   label;
  int      depth;
  int      num_descendants;

public:
  Node(Node* myParent, int my_edge_idx){
    parent = myParent;
    for(int i = 0; i < NUM_CHARS; i++){
      starts[i]   = NULL;
      stops[i]    = NULL;
      children[i] = NULL;
      global[i]   = false;
    }
    sfx_link = NULL;
    edge_idx = my_edge_idx;
    depth    = -1;
    num_descendants = 0;
  }

  ~Node(){
    for(int i = 0; i < NUM_CHARS; i++){
      delete starts[i];
      if(!global[i])
	delete stops[i];
      if(children[i] != NULL)
	delete children[i];
    }
  }

  void setGlobalStop(int edge_idx, size_t*& stop){
    delete stops[edge_idx];
    stops[edge_idx]  = stop;                      
    global[edge_idx] = true;
  }

  void unsetGlobalStop(int edge_idx){
    if(!global[edge_idx])
      printErrorAndDie("unsetGlobalStop()", "Cannot unset a non-global stop");
    global[edge_idx] = false;
    size_t val       = *stops[edge_idx];
    stops[edge_idx]  = new size_t;
    *stops[edge_idx] = val;
  }

  int    getEdgeIndex()              { return edge_idx;                             }
  Node*  getParent()                 { return parent;                               }
  int    getDepth()                  { return depth;                                }
  size_t getStart(int edge_idx)      { return *starts[edge_idx];                    }
  size_t getStop(int edge_idx)       { return *stops[edge_idx];                     }
  size_t getNumChars(int edge_idx)   { return *stops[edge_idx]-*starts[edge_idx]+1; } 
  void   setLabel(size_t my_label)   { label = my_label;  }
  size_t getLabel()                  { return label;      }
  void   addDescendants(int ndesc)   { num_descendants += ndesc;                    }
  int    getNumDescendants()         { return num_descendants;                      }
  bool   isLeaf(int edge_idx)        { return children[edge_idx] == NULL;           }
  bool   isRoot()                    { return parent == NULL;                       }
  bool   edgeExists(int edge_idx)    { return starts[edge_idx] != NULL;             }
  Node*  getChildNode(int edge_idx)  { return children[edge_idx];                   }
  void   setSuffixLink(Node* link)   { sfx_link = link;                             }
  Node*  getSuffixLink()             { return sfx_link;                             }

  std::vector<Node*> getChildren(){
    std::vector<Node*> my_children;
    for(int i = 0; i < NUM_CHARS; i++)
      if(children[i] != NULL)
	my_children.push_back(children[i]);
    return my_children;
  }

  void addEdge(int edge_idx, size_t start, size_t stop){
    starts[edge_idx]  = new size_t;
    stops[edge_idx]   = new size_t;
    *starts[edge_idx] = start;
    *stops[edge_idx]  = stop;
  }

  void extendEdge(int edge_idx);
  void splitEdge(int edge_idx, size_t pos_idx, size_t char_idx, std::vector<int>& ids, size_t*& global_stop, std::string& tokens, Node*& new_node);
  void print(std::ostream& out, int depth, std::string& word);
  void dfsProcess(size_t& cur_label, int cur_depth, std::vector< std::pair<int, Node*> >& stack, int max_depth, std::vector<Node*>& suffixes);
};

class SuffixTree {
private:
  std::string tokens;
  int charID[256];
  std::vector<int> ids;
  size_t* global_end;
  std::vector<Node*> nodes; // Indexed by label-1 
  std::vector<Node*> suffixes;

  void stringToIds();
  void createTree();
  void traverseChars(Node* start_node, size_t start_idx, size_t stop_idx, 
		     Node*& end_node, int& end_edge_idx, size_t& end_pos_idx);
  int  extendCharacter(Node* node, int edge_idx, size_t pos_index, int char_id, size_t char_idx, 
		       Node*& ins_node, int& ins_edge_idx, size_t& ins_pos_idx, size_t*& global_stop);
  void getLinkInfo(Node* prev_node, int edge_idx, size_t pos_idx, Node*& node, size_t& gamma_start, size_t& gamma_stop); 
  int  singleExtension(size_t i, size_t j, size_t*& global_stop, Node* prev_node, int edge_idx, size_t pos_idx, 
		       Node*& new_node, int& new_edge_idx, size_t& new_pos_idx);
  void dfsProcess();

public:
  Node* root;
  SuffixTree(std::string myString);
  ~SuffixTree();
  std::string getPathString(Node* node);
  Node* getSuffix(int idx);
  void printTree(std::ostream& out);
  std::vector<Node*>& getNodes(){ return nodes; }
};

#endif
