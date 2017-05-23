#ifndef HVC_P_H_
#define HVC_P_H_


typedef struct dlnode {
  double x[4];                    // The data vector              
  struct dlnode * closest[2]; // closest[0] == cx, closest[1] == cy
  struct dlnode * cnext[2]; //current next
  struct dlnode * head[2]; //lowest (0 - x, 1 - y)
  
  double area;
  double volume;
  double lastSlicez;
  
  struct dlnode * next[4]; //keeps the points sorted according to coordinates 2,3 and 4 (in the case of 2 and 3, only the points swept by 4 are kept)
  struct dlnode *prev[4]; //keeps the points sorted according to coordinates 2 and 3 (excepto o sentinela 3)
  
  //prev[0] and next[0] are used to keep a history of the order in which points were added to the list
  
  double hvolume; //4D
  int ndomr;    //number of dominators
  struct dlnode * domr; //dominator
  int id;
  
  double oldvolume; //HVC-ONLY 4D-U-ADD 
  
} dlnode_t;


int preprocessing(dlnode_t * list); //returns the number of dominated points found
double hvc3d(dlnode_t * list, int considerDominated);
int addToDataStructure(dlnode_t * list, dlnode_t * new, int determineInsertionPoints);
double updateContributions(dlnode_t * list, dlnode_t * newp, int adding);
// void saveContributions(dlnode_t * list, double * contribs, int d);
void removeFromDataStructure(dlnode_t * list, dlnode_t * old);

dlnode_t * point2Struct(dlnode_t * list, dlnode_t * p, double * v, int d);
void addToHistory(dlnode_t * list, dlnode_t * newp);
void removeFromHistory(dlnode_t * oldp);

double oneContribution3d(dlnode_t * list, dlnode_t * new); //private
double hv3dplus(dlnode_t * list); //private

dlnode_t *
setup_cdllist(double * data, int naloc, int n, int d, double *ref);
void free_cdllist(dlnode_t * list);



#endif
