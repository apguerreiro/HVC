 
/*************************************************************************

 hvc.c

 ---------------------------------------------------------------------

                        Copyright (c) 2013, 2016, 2017
                Andreia P. Guerreiro <apg@dei.uc.pt>
             

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 3 of the
 License.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------

 Reference:

 [1] A. P. Guerreiro, C. M. Fonseca, “Computing and Updating Hypervolume Contributions in Up to Four Dimensions”, CISUC Technical Report TR-2017-001, University of Coimbra, 2017

*************************************************************************/ 

#include "hvc.h"
#include "hvc-private.h"
#include "avl.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <string.h>



#if __GNUC__ >= 3
# define __hv_unused    __attribute__ ((unused))
#else
# define __hv_unused    /* no 'unused' attribute available */
#endif

#define AVL_DEPTH

/* ---------------------------------- Auxiliar Functions ---------------------------------------*/


static inline double max(double a, double b){
    return (a > b) ? a : b;
}


static inline double min(double a, double b){
    return (a < b) ? a : b;
}


static void join(double * p1, double * p2, double * pr){
    
    pr[0] = max(p1[0], p2[0]);
    pr[1] = max(p1[1], p2[1]);
    
}


//3D points
static inline int lexicographicLess(double * a, double * b){
    return (a[2] < b[2] || (a[2] == b[2] && (a[1] < b[1] || (a[1] == b[1] && a[0] <= b[0]))));
}




/* ---------------------------------- Data Structures Functions ---------------------------------------*/


/* next[0] and prev[0] are used to record the order in which points were added 
 * to the data structure. In the case of points given in the initialization phase,
 * then its order is the input order (next[0] and prev[0] are set in setup_cdllist).
 * Functions addToHistory and removeFromHistory are used to record this information
 * and to remove (permanently) the points from the data structure, respectively */
void addToHistory(dlnode_t * list, dlnode_t * newp){
    newp->next[0] = list;
    newp->prev[0] = list->prev[0];
    
    list->prev[0] = newp;
    newp->prev[0]->next[0] = newp;
}


/* Removes 'oldp' permanently from the data structure insertion history.
 * Assumes oldp is not in the hvc data structure (no 'closest' is pointing
 * to oldp nor any 'next[i]'/'prev[i]', for i > 2.)
 */
void removeFromHistory(dlnode_t * oldp){
    oldp->prev[0]->next[0] = oldp->next[0];
    oldp->next[0]->prev[0] = oldp->prev[0];
}




static dlnode_t * initSentinels(dlnode_t * list, const double * ref, int d){
 
    dlnode_t * s1 = list;
    dlnode_t * s2 = list + 1;
    dlnode_t * s3 = list + 2;
    
    
    s1->x[0] = -DBL_MAX;
    s1->x[1] = ref[1];
    s1->x[2] = -DBL_MAX;
    s1->x[3] = -DBL_MAX;
    s1->closest[0] = s2;
    s1->closest[1] = s1;  //to avoid errors in debug prints

    s1->area = 0; //HVC-ONLY
    s1->volume = 0; //HVC-ONLY
    s1->head[0] = s1->head[1] = s1; //HVC-ONLY
    s1->oldvolume = 0; //HVC-ONLY 4D-U-ADD
    
    s1->next[2] = s2;
    s1->next[3] = s2;
    
    s1->cnext[0] = NULL;
    
    s1->prev[2] = s3;
    s1->prev[3] = s3;
    s1->ndomr = 0;
    s1->domr = NULL; //HVC-ONLY
    s1->id = -1;  //HVC-ONLY
    
    
    s2->x[0] = ref[0];
    s2->x[1] = -DBL_MAX;
    s2->x[2] = -DBL_MAX;
    s2->x[3] = -DBL_MAX;
    s2->closest[0] = s2; 
    s2->closest[1] = s1; 
    s2->area = 0; //HVC-ONLY
    s2->volume = 0; //HVC-ONLY
    s2->oldvolume = 0; //HVC-ONLY 4D-U-ADD
    s2->head[0] = s2->head[1] = s2; //HVC-ONLY
    
    s2->next[2] = s3;
    s2->next[3] = s3;
    s2->cnext[1] = NULL;  
    s2->cnext[0] = NULL; 
    
    s2->prev[2] = s1;
    s2->prev[3] = s1;
    s2->ndomr = 0;
    s2->domr = NULL; //HVC-ONLY
    s2->id = -2; //HVC-ONLY
    
    
    s3->x[0] = -INT_MAX;
    s3->x[1] = -INT_MAX;
    s3->x[2] = ref[2];
    if(d == 4)
        s3->x[3] = ref[3];
    else
        s3->x[3] = - DBL_MAX;
    s3->closest[0] = s2;
    s3->closest[1] = s1;
    s3->area = 0; //HVC-ONLY
    s3->volume = 0; //HVC-ONLY
    s3->oldvolume = 0; //HVC-ONLY 4D-U-ADD
    s3->head[0] = s3->head[1] = s3; //HVC-ONLY
    
    s3->next[2] = s1;
    s3->next[3] = NULL;
    s3->cnext[1] = NULL;
    s3->cnext[0] = NULL;
    
    s3->prev[2] = s2;
    s3->prev[3] = s2;
    s3->ndomr = 0;
    s3->domr = NULL; //HVC-ONLY
    s3->id = -3; //HVC-ONLY
    
    
    return s1;
    
}





static void clearPoint(dlnode_t * list, dlnode_t * p){
    
    p->closest[1] = list;
    p->closest[0] = list->next[2]; 
    
    /* because of printfs */
    p->cnext[1] = list;
    p->cnext[0] = list->next[2];
    
    p->head[0] = p->cnext[0]; //HVC-ONLY
    p->head[1] = p->cnext[1]; //HVC-ONLY
    /* until here */
    
    p->area = 0; //HVC-ONLY
    p->volume = 0; //HVC-ONLY
    p->oldvolume = 0; //HVC-ONLY 4D-U-ADD
    p->hvolume = 0; //HVC-ONLY 4D
    p->lastSlicez = p->x[2]; //HVC-ONLY
    p->ndomr = 0;
    p->domr = NULL; //HVC-ONLY
    
}


dlnode_t * point2Struct(dlnode_t * list, dlnode_t * p, double * v, int d){
    
    int i;
    for(i = 0; i < d; i++)
        p->x[i] = v[i];

    clearPoint(list, p);
    
    return p;
    
}

/* ---------------------------------- Update data structure ---------------------------------------*/




static void addToZ(dlnode_t * new){
    
    new->next[2] = new->prev[2]->next[2]; //in case new->next[2] was removed for being dominated
    
    new->next[2]->prev[2] = new;
    new->prev[2]->next[2] = new;
}


static void removeFromz(dlnode_t * old){
    
    old->prev[2]->next[2] = old->next[2];
    old->next[2]->prev[2] = old->prev[2];
}



/* check if 'new' is dominated, find cx and cy ('closest[0]' and 'closest[1]')
 * of the 'new' point and find where to insert 'new' in the list sorted
 * lexicographically by z, y, x (prev[2]/next[2]).
 */
static void setupZandClosest(dlnode_t * list, dlnode_t * new){
    
            
    double * closest1 = (double *) (list);
    double * closest0 = (double *) (list->next[2]);

    dlnode_t * q = (list->next[2]->next[2]);
    double * newx = new->x;
    while(lexicographicLess(q->x, newx)){
        if(q->x[0] <= newx[0] && q->x[1] <= newx[1]){
                
            new->ndomr += 1;
            new->domr = q; //HVC-ONLY 4D
//                 return new;
                
        }else if(q->x[1] < newx[1] && (q->x[0] < closest0[0] || (q->x[0] == closest0[0] && q->x[1] < closest0[1]))){
            closest0 = (double *) q;
        }else if(q->x[0] < newx[0] && (q->x[1] < closest1[1] || (q->x[1] == closest1[1] && q->x[0] < closest1[0]))){
            closest1 = (double *) q;
        }
        
        q = q->next[2];
    }
    
    new->closest[0] = new->cnext[0] = (dlnode_t *) closest0;
    new->closest[1] = new->cnext[1] = (dlnode_t *) closest1;
    
    new->prev[2] = q->prev[2];
    new->next[2] = q;
    
}



/**
 * Inserts 'new' in the data structure.
 * (Adds 'new' to the list sorted lexicographically by z, determines its outer
 * delimiters cx and cy (closest[0] and closest[1], respectively) if needed,
 * and updates the outer delimiters of the points previously in the list and
 * that are lexicographically above 'new' in z.
 * Points (that become) dominated by more than one point are removed from list z.
 * 
 * param list - list of points
 * param new - the point to be added
 * determineInsertionPoints - boolean indicating that: 
 *      0 - new->closest[0] and new->closest[1] are known as well as where to
 *          insert in list z (new->prev[2], new->next[2])
 * 
 * Note: In doubt, set determineInsertionPoints to 1.
 *      
 * returns: the number of points dominated by 'new'
 */
int addToDataStructure(dlnode_t * list, dlnode_t * new, int determineInsertionPoints){

    if(determineInsertionPoints) setupZandClosest(list, new);
    addToZ(new);
    
    dlnode_t * p = new->next[2];
    dlnode_t * stop = list->prev[2];
    int ndom = 0;
    int allDelmtrsVisited = 0;
    while(p != stop && allDelmtrsVisited < 2){ //HVC-ONLY 4D (allDelmtrsVisited < 2 instead of !allDelmtrsVisited)
        
        if(p->x[0] <= new->x[0] && p->x[1] <= new->x[1] && (p->x[0] < new->x[0] || p->x[1] < new->x[1])){
            allDelmtrsVisited += 1;
        }else {
            if(allDelmtrsVisited == 0 || p->ndomr > 0){ //HVC-ONLY 4D (if condition)
                if(new->x[0] <= p->x[0]){
                    //new <= p
                    if(new->x[1] <= p->x[1]){
                        p->ndomr++;
                        p->domr = new; //HVC-ONLY 4D
                        ndom += 1;
    //                     removeFromz(p); //HV-ONLY (does not need dominated to compute HV)
                        
                    }else if(new->x[0] < p->x[0] && (new->x[1] < p->closest[1]->x[1] || (new->x[1] == p->closest[1]->x[1] && (new->x[0] < p->closest[1]->x[0] || (new->x[0] == p->closest[1]->x[0] && new->x[2] < p->closest[1]->x[2]))))){
                        p->closest[1] = new;
                    }
                }else if(new->x[1] < p->x[1] && (new->x[0] < p->closest[0]->x[0] || (new->x[0] == p->closest[0]->x[0] && (new->x[1] < p->closest[0]->x[1] || (new->x[1] == p->closest[0]->x[1] && new->x[2] < p->closest[0]->x[2]))))){
                    p->closest[0] = new;
                }
                
                
                if(p->ndomr > 1){ //HVC-ONLY 4D
                    removeFromz(p); //HVC-ONLY 4D
                }
            }
        }
        p = p->next[2];
    }
    
    return ndom;
}




/* ---------------------------------- Sort ---------------------------------------*/
//lexicographic order of coordinates (z,y,x)
static int compare_point3d(const void *p1, const void* p2)
{
    int i;
    for(i = 2; i >= 0; i--){
        double x1 = (*((const double **)p1))[i];
        double x2 = (*((const double **)p2))[i];

        if(x1 < x2)
            return -1;
        if(x1 > x2)
            return 1;
    }
    return 0;
}



//lexicographic order of coordinates (w,z,y,x)
static int compare_point4d(const void *p1, const void* p2)
{
    int i;
    for(i = 3; i >= 0; i--){
        double x1 = (*((const double **)p1))[i];
        double x2 = (*((const double **)p2))[i];

        if(x1 < x2)
            return -1;
        if(x1 > x2)
            return 1;
    }
    return 0;
}






/*
 * Setup circular double-linked list in each dimension
 */
dlnode_t *
setup_cdllist(double * data, int naloc, int n, int d, double *ref)
{
    int di = d-1;
    
    dlnode_t * head = (dlnode_t *) malloc((naloc+3) * sizeof(dlnode_t));
    int i;
    dlnode_t * list = head;
    initSentinels(head, ref, d);
    
    if(n > 0){
        double **scratchd;
        int j;
        scratchd = malloc(n * sizeof(double*));
        double * data2 = (double *) malloc(d * n * sizeof(double));
        
        for (i = 0; i < n; i++) {
            scratchd[i] = &data[d*i];   
        }
        
        if(d == 3)
            qsort(scratchd, n, sizeof(double*), compare_point3d);
        else if(d == 4)
            qsort(scratchd, n, sizeof(double*), compare_point4d);
        
        for(i = 0; i < n; i++){
            for(j = 0; j < d; j++){
                data2[d * i + j] = scratchd[i][j];
            }
        }
        
        dlnode_t ** scratch = (dlnode_t **) malloc(n * sizeof(dlnode_t *));

        dlnode_t ** originalOrder = (dlnode_t **) malloc((n+1) * sizeof(dlnode_t*)); //to keep track where the i-th point is in scratch
        originalOrder[0] = list;
        
        for (i = 0; i < n; i++) {
            scratch[i] = point2Struct(list, head+i+3, &data2[i*d], d);
            scratch[i]->id = (scratchd[i]-data)/d;
            originalOrder[scratch[i]->id+1] = scratch[i];
        }
        
        free(scratchd);
        
        
        dlnode_t * s = head->next[di];
        s->next[di] = scratch[0];
        scratch[0]->prev[di] = s;

                
        for(i = 0; i < n-1; i++){
            scratch[i]->next[di] = scratch[i+1];
            scratch[i+1]->prev[di] = scratch[i];
            
            originalOrder[i+1]->prev[0] = originalOrder[i];
            originalOrder[i+1]->next[0] = originalOrder[i+2];
        }
        
        originalOrder[n]->prev[0] = originalOrder[n-1];
        originalOrder[n]->next[0] = originalOrder[0];
            
        list->next[0] = originalOrder[1];
        list->prev[0] = originalOrder[n];
        
        
        s = head->prev[di];
        s->prev[di] = scratch[n-1];
        scratch[n-1]->next[di] = s;
        
        free(scratch);
        free(data2);
        
        free(originalOrder);
    }else{
        
        list->next[0] = list->prev[0] = list;
    }
    
    return head;
}



void free_cdllist(dlnode_t * list)
{
    free(list);
}




/* ---------------------------------- Preprocessing ---------------------------------------*/


static int compare_tree_asc_y( const void *p1, const void *p2)
{
    const double x1= *((const double *)p1+1);
    const double x2= *((const double *)p2+1);

    if (x1 < x2)
        return -1;
    else if (x1 > x2)
        return 1;
    else return 0;
}




static inline double *node_point(const avl_node_t *node)
{
    return (double*) node->item;
}


int preprocessing(dlnode_t * list){


//     restartListy(list);
    
    avl_tree_t * avltree = avl_alloc_tree ((avl_compare_t) compare_tree_asc_y, NULL);
    
    dlnode_t * p = list;
    int ndom = 0;
    

    avl_node_t * node = malloc(sizeof(avl_node_t));
    node = avl_init_node(node, p->x);
    avl_insert_top(avltree, node);
    p = p->next[2];
    
    
    avl_node_t * nodeaux = malloc(sizeof(avl_node_t));
    nodeaux = avl_init_node(nodeaux, p->x);
    
    avl_insert_before(avltree, node, nodeaux);
    
    p = p->next[2];
    

    dlnode_t * stop = list->prev[2];
    double * prev;
    double * point;
    
    while(p != stop){
        
        node = malloc(sizeof(avl_node_t));
    
        node = avl_init_node(node, p->x);
            
        
        if(avl_search_closest(avltree, p->x, &nodeaux) == 1)
            nodeaux = nodeaux->next;
        point = node_point(nodeaux);

        if(point[1] == p->x[1] && point[0] <= p->x[0]){
            nodeaux = nodeaux->next;
        }
        
        prev = node_point(nodeaux->prev);
        
        if(prev[0] <= p->x[0] && prev[1] <= p->x[1] && prev[2] <= p->x[2]){
//             p->ndomr = 2; //to know that p is not in the hvc data structure
            p->ndomr = 1;
            p->domr = list; //HVC-ONLY //to avoid errors in debug prints
            free(node);
            ndom++;
        }else{
            while(node_point(nodeaux)[0] >= p->x[0]){
                
                nodeaux = nodeaux->next;
                avl_delete_node(avltree, nodeaux->prev);            
            }
            
            avl_insert_before(avltree, nodeaux, node);
            p->closest[0] = node->prev->item;
            p->closest[1] = node->next->item;
        }
        p = p->next[2];
    }
    
    
    avl_free_tree(avltree);
    return ndom;
    
}




/* ----------------------Hypervolume Indicator Algorithms ---------------------------------------*/



static void restartListy(dlnode_t * list){
    
    list->next[2]->cnext[1] = list; //link sentinels (p: -inf ref[1] -inf e p->next[2]: ref[0] -inf -inf)
    list->cnext[0] = list->next[2];
    
    
}


/*
 * Returns the area dominated by 'p', by sweeping points in asceding order of
 * coordinate 'di' (which is either 0 or 1, i.e., x or y), starting from
 * point 's' and stopping when a point nondominated by 'p' and with coordinate
 * di higher than that of 'p' on the (x,y)-plane is reached.
*       p  - the point whose contributions in 2D is to be computed
*       di - dimension used for sweeping points (in asceding order)
*       s  - outer delimiter of p (with lower di-coordinate value than p) from
*            which to start the sweep. 
*       u  - The delimiter of p with lowest di-coordinate which is not s. If p has
*            inner delimiters, then u is the inner delimiter of p with lowest
*            di-coordinate, otherwise, u is the outer delimiter with higher
*            di-coordinate than p. (Note: u is given because of the cases for
*            (which p has inner delimiter(s). When p does not have any inner delimiters
*            then u=s->cnext[di], but this is not true when inner delimiters exist)
*/
static double computeAreaSimple(double * p, int di, dlnode_t * s, dlnode_t * u){

    int dj = 1 - di;
    double area = 0;
    dlnode_t * q = s;
    area += (q->x[dj] - p[dj]) * (u->x[di] - p[di]);
    
    while(p[dj] < u->x[dj]){

        q = u;
        u = u->cnext[di];
        area += (q->x[dj] - p[dj]) * (u->x[di] - q->x[di]);

    }
 
    
    return area;
    
}




static void setupNDPoint(dlnode_t * p){
    
    
    p->cnext[0] = p->closest[0];
    p->cnext[1] = p->closest[1];
    
    p->head[1] = p->cnext[0]->cnext[1];
    p->head[0] = p->cnext[1]->cnext[0];
    
    
}



static void setupDomPoint(dlnode_t * p){
    
    p->cnext[0] = p->closest[0];
    p->cnext[1] = p->closest[1];
    
    dlnode_t * domr = p->domr;
    
    if(p->cnext[0]->cnext[1] == domr)
        p->head[1] = domr->head[1];
    else
        p->head[1] = p->cnext[0]->cnext[1];
        
    if(p->cnext[1]->cnext[0] == domr)
        p->head[0] = domr->head[0];
    else
        p->head[0] = p->cnext[1]->cnext[0];
    
}





static void updateVolume(dlnode_t * p, double z){

    p->volume += p->area * (z - p->lastSlicez);
    p->lastSlicez = z;
    
}



static void updateVolumeSimple(double * x, int di, dlnode_t * u){

    int dj = 1 - di;
    dlnode_t * q = u->cnext[dj];
    
    updateVolume(q, x[2]);
    while(x[dj] < u->x[dj]){
        q = u;
        u = u->cnext[di];
        updateVolume(q, x[2]);
    }
    updateVolume(u, x[2]);
    
}



static void addNDPoint(dlnode_t * p){
    
    
    //update 'head's of neighbour of 'p' only if 'p' is nondominated
    if(p->ndomr == 0){
        
        if(p->cnext[0]->head[1]->x[1] >= p->x[1]){
            p->cnext[0]->head[1] = p;
            p->cnext[0]->head[0] = p->cnext[0]->cnext[0];
        }else{
            dlnode_t * q = p->cnext[0]->head[0];
            while(q->x[1] >= p->x[1]){
                q = q->cnext[0];
            }
            
            p->cnext[0]->head[0] = q;
            q->cnext[1] = p;
        }
        
        
        if(p->cnext[1]->head[0]->x[0] >= p->x[0]){
            p->cnext[1]->head[0] = p;
            p->cnext[1]->head[1] = p->cnext[1]->cnext[1];
        }else{
            dlnode_t * q = p->cnext[1]->head[1];
            while(q->x[0] >= p->x[0]){
                q = q->cnext[1];
            }
            
            p->cnext[1]->head[1] = q;
            q->cnext[0] = p;
        }
    }

    if(p->cnext[0]->cnext[1]->x[1] > p->x[1] || (p->cnext[0]->cnext[1]->x[1] == p->x[1] && p->cnext[0]->cnext[1]->x[0] > p->x[0]))
        p->cnext[0]->cnext[1] = p;

    if(p->cnext[1]->cnext[0]->x[0] > p->x[0] || (p->cnext[1]->cnext[0]->x[0] == p->x[0] && p->cnext[1]->cnext[0]->x[1] > p->x[1]))
        p->cnext[1]->cnext[0] = p;   
    
     
}





static void addDomPoint(dlnode_t * p){
    
    p->head[0] = p->cnext[0];
    p->head[1] = p->cnext[1];

    if(p->cnext[0]->cnext[1] == p->domr)
        p->domr->head[1] = p;
    else
        p->cnext[0]->cnext[1] = p;

    if(p->cnext[1]->cnext[0] == p->domr)
        p->domr->head[0] = p;
    else
        p->cnext[1]->cnext[0] = p;   
    
     
}



/*
 * Compute all contributions.
 *   list - list of points
 *   considerDominated - 1 indicates whether dominated points are admitted, i.e, they decrease
 *                      the contribution of the only point dominating them (domr).
 *                      0 indicates that dominated points should be ignored.
 * 
 * This code corresponds to HVC3D algorithm as described in the paper. 
 * The main difference is that, althought each p maintains a list of its delimiters in 2D (p.L)
 * through 'head', this list does not contain the outer delimiters (these are accessible
 * through 'cnext') nor does it copy points to its list. Only one copy of each point exist in
 * the whole program.
 */
double hvc3d(dlnode_t * list, int considerDominated){
    
    dlnode_t * p;
    dlnode_t * q;

    double area = 0;
    double volume = 0;
    double x[3];
    
    restartListy(list);
    p = list->next[2]->next[2];
    
    dlnode_t * stop = list->prev[2];
    while(p != stop){
        p->area = 0;
        p->volume = 0;
        //p->oldvolume = 0;
        p->lastSlicez = p->x[2];
        //volume += area * (p->x[2]- p->prev[2]->x[2]);
    
        if(p->ndomr < 1){

            setupNDPoint(p);
            
            updateVolumeSimple(p->x, 1, p->head[1]);
            p->area = computeAreaSimple(p->x, 1, p->cnext[0], p->head[1]);
            area += p->area;

            q = p->cnext[0]; 
            x[0] = q->x[0]; x[1] = p->x[1]; x[2] = p->x[2]; // join(p,q) - (x[2] is not important)
            q->area -= computeAreaSimple(x, 0, p->head[1], q->head[0]);
    
            q = p->cnext[1]; 
            x[0] = p->x[0]; x[1] = q->x[1];
            q->area -= computeAreaSimple(x, 1, p->head[0], q->head[1]);
    
            addNDPoint(p);
    
        }else if(considerDominated && p->ndomr == 1){
            //if dominated points must be taken into account, then remove
            //the area of their single dominating point that they dominate
            
            updateVolume(p->domr, p->x[2]);
            setupDomPoint(p);
            p->domr->area -= computeAreaSimple(p->x, 1, p->cnext[0], p->head[1]);
            
            addDomPoint(p);
        }/*else{
            //too dominated
        }*/
        
        volume += area * (p->next[2]->x[2]- p->x[2]);
        p = p->next[2];
        
    }
    
    setupNDPoint(p);    
    updateVolumeSimple(p->x, 1, p->head[1]);
    return volume;
    
}




/* Compute all hypervolume contributions in d=4 by iteratively
 * computing all hypervolume contributions in d=3 (using hvc3d)
 */
double hvc4dR(dlnode_t * list)
{
    double height = 0, volume = 0, hv = 0;
    
    dlnode_t * last = list->prev[3];
    dlnode_t * new = list->next[3]->next[3];
    dlnode_t * stop = list->prev[2];
    dlnode_t * p;
    
    int considerDominated = 1;
    
    while(new != last){
       
        addToDataStructure(list, new, 1);
        
        volume = hvc3d(list, considerDominated);  // compute hv indicator and contributions in d=3, in linear time 
        p = list->next[2]->next[2];
        while(p != stop){
            p->hvolume += p->volume * (new->next[3]->x[3] - new->x[3]);
            p = p->next[2];
        }
        
        height = new->next[3]->x[3] - new->x[3];
        hv += volume * height;                // update hypervolume in d=4
        new = new->next[3];
    }
        
    return hv;
}



static void setupPoint(dlnode_t * p){
    
    p->cnext[0] = p->closest[0];
    p->cnext[1] = p->closest[1];
    
    dlnode_t * head1 = p->cnext[0]->cnext[1];
    dlnode_t * head0 = p->cnext[1]->cnext[0];
    
    if(p->ndomr > 0 && head1 == p->domr){
        head1 = head1->head[1];
    }
    if(p->ndomr > 0 && head0 == p->domr){
        head0 = head0->head[0];
    }
    
    p->head[0] = head0;
    p->head[1] = head1;
    
}


static void clearAreaVol(dlnode_t * p){
    
    p->oldvolume = p->volume;
    p->volume = 0;
    p->area = 0;
    
}





static double computeArea(double * p, int di, dlnode_t * l1, dlnode_t * l2, dlnode_t * u1){

    int dj = 1 - di;
    dlnode_t * q;
    double area = 0;
    double lastx = l1->x[dj];
    
    q = l2;    

    while(q->x[dj] > p[dj]){
        area += (lastx - q->x[dj]) * (q->x[di] - p[di]);
        
        lastx = q->x[dj];
        q = q->cnext[di];
        
    }
    
    area += (lastx - p[dj]) * (min(q->x[di], u1->x[di]) - p[di]);
    
    return area;
    
}


static void updateDominatorAreaVol(dlnode_t * p, dlnode_t * domr, dlnode_t * new, int di){
    
    int dj = 1 - di;
    dlnode_t * l2 = p->head[di];
    
    double pj[3]; pj[2] = 0;
    join(p->x, new->x, pj);
    
    double area = computeArea(pj, di, p->cnext[dj], l2, p->cnext[di]);
    p->area = area;
    updateVolume(domr, p->x[2]);
    domr->area = domr->area - area;

}





static void updateDominatorHeads(dlnode_t * p, dlnode_t * domr){
    
        
        if(p->x[0] < domr->head[0]->x[0] || (p->x[0]== domr->head[0]->x[0] && p->x[1] < domr->head[0]->x[1])){
            domr->head[0] = p;
        }
        if(p->x[1] < domr->head[1]->x[1] || (p->x[1] == domr->head[1]->x[1] && p->x[0] < domr->head[1]->x[0])){
            domr->head[1] = p;
        }
}





/*
 * Update the (3d) contribution of 'q' and subtract from the 2d contribution
 * of 'q' the area dominated by 'p'. Use coordinate di to sweep points.
 * p - new point that partially dominated the area exclusively dominated by 'q'
 * q - point in the hvc data strucutre
 * di - coordinate used to sweep the points in ascending order in the computation
 *      of the joint 2d contribution (dominated by p and q)
 */
static double updateVolPartialArea(dlnode_t * p, dlnode_t * q, int di){
    
    double area = 0;
    
    if(q->x[2] != -DBL_MAX){
        
        int dj = 1 - di;
        double pj[3]; pj[2] = 0;
        dlnode_t * l2;
        double z = max(p->x[2], q->x[2]);
            
        updateVolume(q, z);
        
        join(p->x, q->x, pj);
        l2 = q->head[di];
        area = computeArea(pj, di, q->cnext[dj], l2, q->cnext[di]);
        q->area -= area;
    }
    
    return area;
    
}


//Example: p = (15 8 7), q = (12 10 2) e new = (9 14 1)
/* Similar to updateVolPartialArea but computes the area exclusively dominated by p, q and new
 */
static double updateVolSubPartialArea(dlnode_t * p, dlnode_t * q, dlnode_t * new, int di){
    double area = 0;
    
    if(q->x[2] != -DBL_MAX){
        
        int dj = 1 - di;
        double pj[2];
        dlnode_t * l2;
        dlnode_t * l1;
        double z = max(p->x[2], q->x[2]);
        
        updateVolume(q, z);
        
        join(p->x, q->x, pj);
        join(pj, new->x, pj);

        l2 = q->head[di];
        if(l2->x[di] < pj[di]){
            while(l2->x[di] < pj[di]){
                l1 = l2;
                l2 = l2->cnext[di];
            }
        }else{
            l2 = q->head[di];
            l1 = q->cnext[dj];
            
        }
        area = computeArea(pj, di, l1, l2, q->cnext[di]);
        q->area -= area;

    }
    
    return area;
    
}


static void updateExistingAreas(dlnode_t * new, dlnode_t * p, int di, int adding){
    int dj = 1 - di;
    dlnode_t * q = p->cnext[dj];
    int factor = (adding) ? -1 : 1;
    
    if(p->ndomr == 1){
//         dlnode_t * l2;
        dlnode_t * domr = p->domr;
        updateDominatorAreaVol(p, domr, new, di);
        updateDominatorHeads(p, domr);
        
    }else{
    
    //     updateVolPartialArea(p, p->cnext[dj], dj);
        updateVolSubPartialArea(p, p->cnext[dj], new, dj);
        
        if(p->cnext[dj]->x[dj] >= new->x[dj]){
            q = p->head[di];
            
    //         if(q != NULL){
            if(q != NULL && q->x[2] != -DBL_MAX){
                double maxv = max(new->x[dj], p->x[dj]); //if p is dominated by new, it stops when an outer delimiter of p is reached
//                 while(new->x[dj] <= q->x[dj]){ // && p->x[di] <= q->x[di]
                while(maxv <= q->x[dj]){ // && p->x[di] <= q->x[di]
                    updateVolume(q, p->x[2]);
                    q->volume = q->oldvolume + factor * q->volume;
                    q = q->cnext[di];
                }
            }
                
            else{
                q = p->cnext[di];
            }
            
            if(p->x[dj] >= new->x[dj]){
                updateVolSubPartialArea(p, q, new, di);
            }else{
                updateVolume(q, p->x[2]);
                q->volume = q->oldvolume + factor * q->volume;
            }
        }
    }
}


//assumes new->cnext[1-di] = p
//updates new->head[di] (and new->head[1-di] and new->head[di]->cnext[dj] if it is an inner delimiter of p)
static void updateHeads(dlnode_t * new, dlnode_t * p, int di){

    int dj = 1 - di;
    
    if(new->head[dj]->x[dj] >= p->x[dj]){
        new->head[dj] = p;
        new->head[di] = new->cnext[di];
    }else{
            
        dlnode_t * q = new->head[di];
        while(q->x[dj] >= p->x[dj]){
            q = q->cnext[di];
        }

        new->head[di] = q;
    //         q->cnext[dj] = new->cnext[dj];
        q->cnext[dj] = p; //to assure that the links of inner delimiters is correct
        
    }
    
}





static int interactsIndirectly(dlnode_t * p, dlnode_t * new, int di){

    int dj = 1 - di;
    dlnode_t * next = new->cnext[dj];
    
    if(p->x[dj] >= next->cnext[dj]->x[dj]) return 0; //if p definitely does not interact with new 
                                                    //(there is, at least, two points between them)
    
    dlnode_t * q = next->head[di];
    dlnode_t * h = q;
    
    while(q->x[dj] >= p->x[dj] && q->x[di] <= new->x[di]){
        h = q;
        q = q->cnext[di];
    }
    
    if(q->x[di] > new->x[di]){
        h->cnext[dj] = next->head[di]->cnext[dj];
        next->head[di] = h;
        return 1;
    }
    return 0;
    
    
}


//reconstructing L at z = stop->x[2] (only with points lexicographically less than 'stop')
static void restartBase(dlnode_t * list, dlnode_t * stop){

    dlnode_t * p = list->next[2]->next[2];
    restartListy(list);
    
    while(p != stop){
            //reconstruct
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            /*this is needed in updateContributions when removing a point
             * but it is not necessary in the reconstruction of the base when
             * removing a point from the data structure */
            //p->head[0] = p->cnext[1]->cnext[0];
            //p->head[1] = p->cnext[0]->cnext[1];
            
            p->cnext[0]->cnext[1] = p;
            p->cnext[1]->cnext[0] = p;
            
            
        p = p->next[2];
    }
}


/* reconstructing L at z = stop->x[2] (only with points lexicographically less than 'stop'),
 * including inner delimiters of each point in L */
static void restartContributionsBase(dlnode_t * list, dlnode_t * stop){
//     printf("restartContributionsBase\n");

    dlnode_t * p = list->next[2]->next[2];
    restartListy(list);
    
    while(p != stop){
            //reconstruct
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            if(p->ndomr == 0){
                /*this is needed in updateContributions when removing a point
                * but it is not necessary in the reconstruction of the base when
                * removing a point from the data structure */
                p->head[0] = p->cnext[1]->cnext[0];
                p->head[1] = p->cnext[0]->cnext[1];
                
                p->cnext[0]->cnext[1] = p;
                p->cnext[1]->cnext[0] = p;
                
                updateHeads(p->cnext[0], p, 0);
                updateHeads(p->cnext[1], p, 1);
            }else{
                p->head[0] = p->cnext[0];
                p->head[1] = p->cnext[1];
                
                if(p->domr->cnext[1] == p->cnext[1]){
                    p->domr->head[0] = p;
                }else{
                    p->cnext[1]->cnext[0] = p;
                }
                if(p->domr->cnext[0] == p->cnext[0]){
                    p->domr->head[1] = p;
                }else{
                    p->cnext[0]->cnext[1] = p;
                }
                
            }
            
        p = p->next[2];
    }
}





/* does what setupZandClosest does while reconstructing L at z = new->x[2]
 * (visiting only the points that are lexicographically less than new)*/
static void restartBaseSetupZandClosest(dlnode_t * list, dlnode_t * new){
    
//     printf("restartBaseSetupZandClosest\n");
    dlnode_t * p = list->next[2]->next[2];
    double * closest1 = (double *) (list);
    double * closest0 = (double *) (list->next[2]);

    double * newx = new->x;
    
    restartListy(list);
    while(lexicographicLess(p->x, newx)){
        if(p->ndomr <= 1){ //HVC-ONLY 4D-U-ADD
        
            //reconstruct
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            if(p->ndomr == 0){
                p->head[0] = p->cnext[1]->cnext[0]; // HVC-ONLY UPDATE hvc4dU
                p->head[1] = p->cnext[0]->cnext[1]; // HVC-ONLY UPDATE hvc4dU
                
                p->cnext[0]->cnext[1] = p;
                p->cnext[1]->cnext[0] = p;
                
                updateHeads(p->cnext[0], p, 0); // HVC-ONLY UPDATE hvc4dU
                updateHeads(p->cnext[1], p, 1); // HVC-ONLY UPDATE hvc4dU
            }else{
                p->head[0] = p->cnext[0];
                p->head[1] = p->cnext[1];
                
                if(p->domr->cnext[1] == p->cnext[1]){
                    p->domr->head[0] = p;
                }else{
                    p->cnext[1]->cnext[0] = p;
                }
                if(p->domr->cnext[0] == p->cnext[0]){
                    p->domr->head[1] = p;
                }else{
                    p->cnext[0]->cnext[1] = p;
                }
                
//                 updateDominatorHeads(p, p->domr); // HVC-ONLY UPDATE hvc4dU
                /*
                if(p->domr->head[1] != p){
                    p->cnext[0]->cnext[1] = p;
                }
                if(p->domir->head[0] != p){
                    p->cnext[1]->cnext[0] = p;
                }
                */
                
            }
            
            
            //setupZandClosest
            if(p->x[0] <= newx[0] && p->x[1] <= newx[1]){
                
                new->ndomr += 1;
                new->domr = p; //HVC-ONLY 4D-U-ADD
    //                 return new;
                    
            }else if(p->x[1] < newx[1] && (p->x[0] < closest0[0] || (p->x[0] == closest0[0] && p->x[1] < closest0[1]))){
                closest0 = (double *) p;
            }else if(p->x[0] < newx[0] && (p->x[1] < closest1[1] || (p->x[1] == closest1[1] && p->x[0] < closest1[0]))){
                closest1 = (double *) p;
            }
        }else{ //HVC-ONLY 4D-U-ADD
            p->prev[2]->next[2] = p->next[2]; //HVC-ONLY 4D-U-ADD
            p->next[2]->prev[2] = p->prev[2]; //HVC-ONLY 4D-U-ADD
        } //HVC-ONLY 4D-U-ADD
        
        p = p->next[2];
    }
    
    new->closest[0] = (dlnode_t *) closest0;
    new->closest[1] = (dlnode_t *) closest1;
    
    new->prev[2] = p->prev[2];
    new->next[2] = p;

}




static void addPoint(dlnode_t * p){
    
    //update heads of neighbour (only do this for nondominated points)
    if(p->ndomr == 0){
        dlnode_t * q = p->cnext[0]->head[0];
        while(q->x[1] >= p->x[1]){
            q = q->cnext[0];
        }
        
        p->cnext[0]->head[0] = q;
        if(q->x[1] >= p->cnext[0]->x[1]){
            q->cnext[1] = p;
        }else{
            p->cnext[0]->head[1] = p;
        }

        q = p->cnext[1]->head[1];
        while(q->x[0] >= p->x[0]){
            q = q->cnext[1];
        }
        
        p->cnext[1]->head[1] = q;
        if(q->x[0] >= p->cnext[1]->x[0]){
            q->cnext[0] = p;
        }else{
            p->cnext[1]->head[0] = p;
        }
    }

    if(p->cnext[0]->cnext[1]->x[1] > p->x[1] || (p->cnext[0]->cnext[1]->x[1] == p->x[1] && p->cnext[0]->cnext[1]->x[0] > p->x[0]))
        p->cnext[0]->cnext[1] = p;

    if(p->cnext[1]->cnext[0]->x[0] > p->x[0] || (p->cnext[1]->cnext[0]->x[0] == p->x[0] && p->cnext[1]->cnext[0]->x[1] > p->x[1]))
        p->cnext[1]->cnext[0] = p;   
     
}


/* compute the area dominated by every delimiter of new's contribution in L
 * and by new and assign it to the corresponding delimiter */
static void restartCoveredAreas(dlnode_t * new){
    dlnode_t * q = new->cnext[0];
    dlnode_t * l2;
    
    //right outer delimiter
    if(q->x[2] != -DBL_MAX){
        clearAreaVol(q);
        double pj[3] = {q->x[0], new->x[1], new->x[2]};
        
        l2 = q->head[0];
        q->area = computeArea(pj, 0, q->cnext[1], l2, q->cnext[0]);
        q->lastSlicez = new->x[2];
    }
    q = q->cnext[1];
    
    //inner delimiters
    while(new->x[0] <= q->x[0] && new->x[1] <= q->x[1]){
        clearAreaVol(q);
        l2 = q->head[1];
        q->area = computeArea(q->x, 1, q->cnext[0], l2, q->cnext[1]);
        q->lastSlicez = new->x[2];    
        q->domr = new;
        
        q = q->cnext[1];
    }
    
    //left outer delimiter
    if(q->x[2] != -DBL_MAX){
        clearAreaVol(q);
        double pj[3] = {new->x[0], q->x[1], new->x[2]};
        l2 = q->head[1];
        q->area = computeArea(pj, 1, q->cnext[0], l2, q->cnext[1]);
        q->lastSlicez = new->x[2];
    }
    
    
}



/*
 * Updates the list L (of nondominated points in 2D) only. Assumes that there is no need
 * to keep track of the delimiters of the points in L (For example, when solving the
 * OneContribution problem regarding 'new' and X, where only the delimiters of 'new' are needed).
 * Therefore, the 'head' of the outer delimiters of p do not have to be updated.
 * Assumes p is nondominated (p->ndomr == 0).
 */
static void addPointNoHeads(dlnode_t * p){
    
    p->head[0] = p->cnext[0];
    p->head[1] = p->cnext[1];
    
    p->cnext[0]->cnext[1] = p;

    p->cnext[1]->cnext[0] = p;   
    
}


/* Assumes the base of new is already setup. The sweep begins in p (which is assumed
 * to be the first point above new in z). Updates the contribution of each delimiter
 * of new above it in z. Updates new's contribution.
 */
static void incrementSlicing(dlnode_t * new, dlnode_t * p){
    
     while(p->x[0] > new->x[0] || p->x[1] > new->x[1]){
        
        if(new->x[0] <= p->x[0] && new->x[1] <= p->x[1]){
            
            if(p->ndomr == 0 || (p->domr->x[1] <= new->x[1] && p->domr->x[0] <= new->x[0])){ //
                    
                setupPoint(p);
                dlnode_t * l2 = p->head[1];
                double area = computeArea(p->x, 1, p->cnext[0], l2, p->cnext[1]);
                p->area = area; //this is not really needed, although it is helpful for debug

                updateVolume(new, p->x[2]);
                new->area -= area;
                addPointNoHeads(p);
                
                updateDominatorHeads(p, new);
            }
            
       
        }else if(p->x[1] < new->x[1] && p->x[0] > new->x[0] && (p->x[0] < new->cnext[0]->x[0] || (p->x[0] == new->cnext[0]->x[0] && p->x[1] < new->cnext[0]->x[1]))){
            p->cnext[1] = p->closest[1];
            //p->cnext[1] = new->cnext[1];
            
            updateVolPartialArea(p, new, 1);
            new->cnext[0] = p;
            updateHeads(new, p, 1);

            p->cnext[1] = new->head[1];
            p->cnext[1]->cnext[0] = p; //to avoid any problems when a dominated point is added
            
        }else if(p->x[0] < new->x[0] && p->x[1] > new->x[1] && (p->x[1] < new->cnext[1]->x[1] || (p->x[1] == new->cnext[1]->x[1] && p->x[0] < new->cnext[1]->x[0]))){
            p->cnext[0] = p->closest[0];

            updateVolPartialArea(p, new, 0);
            new->cnext[1] = p;
            updateHeads(new, p, 0);

            p->cnext[0] = new->head[0];
            p->cnext[0]->cnext[1] = p; //to avoid any problems when a dominated point is added
            
        }
        p = p->next[2];
        
    }
    updateVolume(new, p->x[2]);
    
}


//restartBaseSetupZandClosest from hv-plus (it performs slightly less operations than restartBaseSetupZandClosest for contributions)
static void HVrestartBaseSetupZandClosest(dlnode_t * list, dlnode_t * new){
    

    dlnode_t * p = list->next[2]->next[2];
    double * closest1 = (double *) (list);
    double * closest0 = (double *) (list->next[2]);

    double * newx = new->x;
    
    restartListy(list);
    while(lexicographicLess(p->x, newx)){
            
        if(p->ndomr == 0){
            //reconstruct
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            p->cnext[0]->cnext[1] = p;
            p->cnext[1]->cnext[0] = p;
            
            //setupZandClosest
            if(p->x[0] <= newx[0] && p->x[1] <= newx[1]){
                
                new->ndomr += 1;
                //new->domr = p;
                    //return new;
                    
            }else if(p->x[1] < newx[1] && (p->x[0] < closest0[0] || (p->x[0] == closest0[0] && p->x[1] < closest0[1]))){
                closest0 = (double *) p;
            }else if(p->x[0] < newx[0] && (p->x[1] < closest1[1] || (p->x[1] == closest1[1] && p->x[0] < closest1[0]))){
                closest1 = (double *) p;
            }
        }
        
        p = p->next[2];
    }
    
    new->closest[0] = (dlnode_t *) closest0;
    new->closest[1] = (dlnode_t *) closest1;
    
    new->prev[2] = p->prev[2];
    new->next[2] = p;
    
    
}

/* 
 * Computes the contribution of the new point.
 * Note: assumes that all points are mutually nondominated.
 */
double oneContribution3d(dlnode_t * list, dlnode_t * new){
    
    dlnode_t * p = list;
    double area = 0;
    double volume = 0;
    double x[3];
    
    HVrestartBaseSetupZandClosest(list, new);
    if (new->ndomr > 0)
        return 0;
    
    new->cnext[0] = new->closest[0];
    new->cnext[1] = new->closest[1];

    area = computeAreaSimple(new->x, 1, new->cnext[0], new->cnext[0]->cnext[1]);
    
    
    p = new->next[2];
    double lastz = new->x[2];
    
    while(p->x[0] > new->x[0] || p->x[1] > new->x[1]){
        if(p->ndomr < 1){ //skip dominated points

            volume += area * (p->x[2]- lastz);
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            
            if(p->x[0] >= new->x[0] && p->x[1] >= new->x[1]){
                area -= computeAreaSimple(p->x, 1, p->cnext[0], p->cnext[0]->cnext[1]);
                p->cnext[1]->cnext[0] = p;
                p->cnext[0]->cnext[1] = p;
                
            }else if(p->x[0] >= new->x[0]){
                if(p->x[0] <= new->cnext[0]->x[0]){
                    x[0] = p->x[0]; x[1] = new->x[1]; x[2] = p->x[2]; // join(p,q) - (x[2] is not important)
                    area -= computeAreaSimple(x, 1, new->cnext[0], new->cnext[0]->cnext[1]);
                    p->cnext[0] = new->cnext[0];
                    p->cnext[1]->cnext[0] = p;
                    new->cnext[0] = p;
                }
            }else{
                if(p->x[1] <= new->cnext[1]->x[1]){
                    x[0] = new->x[0]; x[1] = p->x[1]; x[2] = p->x[2]; // join(p,q) - (x[2] is not important)
                    area -= computeAreaSimple(x, 0, new->cnext[1], new->cnext[1]->cnext[0]);
                    p->cnext[1] = new->cnext[1];
                    p->cnext[0]->cnext[1] = p;
                    new->cnext[1] = p;
                }
            }
            lastz = p->x[2];
        }
        p = p->next[2];
        
    }
    volume += area * (p->x[2]- lastz);
    return volume;
    
}



double hv3dplus(dlnode_t * list){
    
    dlnode_t * p = list;
    double area = 0;
    double volume = 0;
    
    restartListy(list);
    p = p->next[2]->next[2];
    
    dlnode_t * stop = list->prev[2];
    
    while(p != stop){
        if(p->ndomr < 1){
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            area += computeAreaSimple(p->x, 1, p->cnext[0], p->cnext[0]->cnext[1]);
            
            p->cnext[0]->cnext[1] = p;
            p->cnext[1]->cnext[0] = p;
        }/*else{
            removeFromz(p);
            p->prev[2]->next[2] = p->next[2];
            p->next[2]->prev[2] = p->prev[2];
        }*/
        
        volume += area * (p->next[2]->x[2]- p->x[2]);
        
        p = p->next[2];
    }
    
    
    return volume;
    
}




double updateContributions(dlnode_t * list, dlnode_t * new, int adding){
    dlnode_t * p;
    
    if(adding){
        restartBaseSetupZandClosest(list, new); //new->next[2] and new->closest[?] only need to be set if
                                                // the point is being processed for the first time (adding)
    }else{
        restartContributionsBase(list, new->prev[2]->next[2]); //so it works even when new is not in list z
    }
    
    new->cnext[0] = new->closest[0];
    new->cnext[1] = new->closest[1];
    
    new->head[1] = new->closest[0]->cnext[1];
    new->head[0] = new->closest[1]->cnext[0];
    new->lastSlicez = new->x[2];
    
    restartCoveredAreas(new);
    //new->oldvolume = new->volume; //it is not really necessary
    new->volume = 0;
    new->area = computeArea(new->x, 1, new->closest[0], new->closest[0]->cnext[1], new->closest[1]);
    p = new->next[2];
    

    
    double factor = (adding) ? -1 : 1;
    dlnode_t * domr;
    
    
    while(p->x[0] > new->x[0] || p->x[1] > new->x[1]){
        if(p->ndomr <= 1){
            
            setupPoint(p);
            updateVolume(new, p->x[2]);
                
            if(new->x[0] <= p->x[0] && new->x[1] <= p->x[1]){
                clearAreaVol(p);
                p->lastSlicez = p->x[2];
                if(p->ndomr == 0){

                    domr = new;
                    updateDominatorAreaVol(p, domr, new, 1);
                    updateDominatorHeads(p, domr);
                    updateExistingAreas(new, p, 0, adding);
                    
                }else{ //p->ndomr == 1
                    
                    domr = p->domr;
                    
                    if(p->x[1] < domr->cnext[1]->x[1] && p->x[0] < domr->cnext[0]->x[0]){ //interacts
                        updateDominatorAreaVol(p, domr, new, 1);
                        updateDominatorHeads(p, domr);
                    }
                }
                
            }else if(p->x[0] > new->x[0]){
                if(p->x[0] < new->cnext[0]->x[0] || (p->x[0] == new->cnext[0]->x[0] && p->x[1] < new->cnext[0]->x[1])){
                    clearAreaVol(p);
                    p->lastSlicez = p->x[2];
                    p->area = updateVolPartialArea(p, new, 1);
                    
                    updateExistingAreas(new, p, 0, adding);
                    updateHeads(new, p, 1);
                    //p->closest[1] = new;

                    new->cnext[0] = p;
                    
                }else if(interactsIndirectly(p, new, 1)){
                //update the areas and volumes of the affected points
                    updateExistingAreas(new, p, 0, adding);
                }
                
            }else if(p->x[1] > new->x[1]){
                if(p->x[1] < new->cnext[1]->x[1] || (p->x[1] == new->cnext[1]->x[1] && p->x[0] < new->cnext[1]->x[0])){
                    clearAreaVol(p);
                    p->lastSlicez = p->x[2];
                    p->area = updateVolPartialArea(p, new, 0);
                    updateExistingAreas(new, p, 1, adding);
                    updateHeads(new, p, 0);
                    //p->closest[0] = new;
                      
                    new->cnext[1] = p;
                    
                }else if(interactsIndirectly(p, new, 0)){
                //update the areas and volumes of the affected points
                    updateExistingAreas(new, p, 1, adding);
                }
            }
            //update L (list of nondominated points in 2D)
            addPoint(p);
        }
        p = p->next[2];
    }
    
    updateVolume(new, p->x[2]);
    
    dlnode_t * q = new->cnext[0];
    while(q->x[0] >= new->x[0]){
        updateVolume(q, p->x[2]);
        q->volume = q->oldvolume + factor * q->volume;
        q = q->cnext[1];
    }
    updateVolume(q, p->x[2]);
    if(q->x[2] != DBL_MAX){
        q->volume = q->oldvolume + factor * q->volume;
    }
    
    if(p->next[2] != list){
        double volume = new->volume;
        new->volume = 0;
        
        new->lastSlicez = p->x[2];
        incrementSlicing(new, p->next[2]);
        p->volume = p->volume + factor * new->volume;
        new->volume = volume;
    }

    return new->volume;
}
    
    

    

/* Compute all hypervolume contributions in d=4 by iteratively
 * computing the contributions that change in d=3
 */
double hvc4dU(dlnode_t * list){
    double height = 0, volume = 0, hv = 0;
    
    dlnode_t * last = list->prev[3];
    dlnode_t * new = list->next[3]->next[3];
    dlnode_t * stop = list->prev[2];
    dlnode_t * p;
    
    int adding = 1;
    int knowInsertionPoints;
    
    while(new != last){
        
        updateContributions(list, new, adding);
        knowInsertionPoints = 1; //insertion points (prev[2]/next[2] and closest[?]) are setup up in updateContributions
        addToDataStructure(list, new, !knowInsertionPoints);
        
        //if(new->ndomr == 1) new->ndomr = 2;
        p = list->next[2]->next[2];
        while(p != stop){
            p->hvolume += p->volume * (new->next[3]->x[3] - new->x[3]);
            p = p->next[2];
        }
        
        height = new->next[3]->x[3] - new->x[3];
        hv += volume * height;                // update hypervolume in d=4
        new = new->next[3];
    }
        
    return hv;
}




/* Contributions are stored in "contribs". The order in which they are stored in based
 * on the history of (addition of) points (stores in the order in which points were added)*/
static void saveContributions(dlnode_t * list, double * contribs, int d){
    
    int i;
    dlnode_t * p;
    p = list->next[0];
    if(d == 3){
        for(i = 0; p != list; i++, p = p->next[0]){
            contribs[i] = p->volume;
        }
    }else if(d == 4){
        for(i = 0; p != list; i++, p = p->next[0]){
            contribs[i] = p->hvolume;
        }
    }
}

/*
static void saveContributions(dlnode_t * list, double * contribs, int d){
    
    dlnode_t * last = list->prev[d-1];
    dlnode_t * p;

    p = list->next[d-1]->next[d-1];
    while(p != last){
        if(d == 3)
            contribs[p->id] = p->volume;
        else if(d == 4)
            contribs[p->id] = p->hvolume;
        p = p->next[d-1];
    }
}*/


/* Compute the hypervolume contributions in d=3,4
 *      data    - array with all points
 *      d       - number of dimensions (either 3 or 4)
 *      n       - number of points
 *      ref     - reference point
 *      recompute - available for d=4 only
 *                  if 0: recompute the hypervolume in d-1 (HVC4D+-R)
 *                  otherwise: compute one contribution in d-1 (HVC4D+-U)
 * TODO: Implement 2D case
 */
double hvc(double *data, int d, int n, double *ref, double * contribs, int recompute){
    double hv = 0;
    dlnode_t * list = setup_cdllist(data, n, n, d, ref);
    
    if(d == 2){
        //hv = hvc2d(data, n, ref, contribs, stdout);
    }else if(d == 3){
        
        preprocessing(list);
        int considerDominated = 0;
        hv = hvc3d(list, considerDominated);
        
    }else if(d == 4){
        if(recompute)
            hv = hvc4dR(list);
        else
            hv = hvc4dU(list);
    }
    
    saveContributions(list, contribs, d);
    
    free_cdllist(list);
    
    return hv;
}



static dlnode_t * leastContributor(dlnode_t * list, int d){
    dlnode_t * p = list->next[d-1]->next[d-1];
    dlnode_t * stop = list->prev[d-1];
    
    double leastcontribution = DBL_MAX;
    dlnode_t * lcontributor = NULL;
    
    while(p != stop){
        if(p->ndomr > 0 || p->volume < leastcontribution){ // ->volume assume d == 3
            leastcontribution = p->volume;
            lcontributor = p;
        }
     
        p = p->next[d-1];
    }
    
    return lcontributor;
}




static int updateInsertingPoints(dlnode_t * old, dlnode_t * p, int di){
    
    if(p->closest[di] == old){
        int dj = 1 - di;
        if(p->x[dj] > old->head[dj]->x[dj]){
            dlnode_t * q = old->head[di];
            while(q->x[dj] >= p->x[dj]){
                q = q->cnext[di];
            }
            p->closest[di] = q;
        }else{
            p->closest[di] = old->cnext[di];
        }
        return 1;
    }
    
    return 0;
    
}

void removeFromDataStructure(dlnode_t * list, dlnode_t * old){
    
    dlnode_t * stop = list->prev[2];
    dlnode_t * p;
    
    restartBase(list, old);
    
    old->cnext[0] = old->closest[0];
    old->cnext[1] = old->closest[1];
    
    old->head[1] = old->closest[0]->cnext[1];
    old->head[0] = old->closest[1]->cnext[0];
    
    p = old->next[2];
    
    
    while(p != stop){

        if(old->x[0] <= p->x[0] && old->x[1] <= p->x[1]){
            p->ndomr--; //this should never happen... we assumed that a nondominated point is removed
                        //only if there is not any dominated points in the data structure
            
        }else if(updateInsertingPoints(old, p, 0)){
            old->cnext[1] = p;
            if(p->closest[0]->x[1] >= old->x[1]){ //if updateInsertingPoints == True then p->closest[0]->x[0] >= old->x[0]
                old->head[0] = p->closest[0];
            }else{
                old->head[0] = old->cnext[0];
                old->head[1] = old->cnext[1];
            }
        }else{
            if(updateInsertingPoints(old, p, 1)){

                old->cnext[0] = p;
                if(p->closest[1]->x[0] >= old->x[0]){  //if updateInsertingPoints == True then p->closest[1]->x[1] >= old->x[1]
                    old->head[1] = p->closest[1];
                }else{
                    old->head[0] = old->cnext[0];
                    old->head[1] = old->cnext[1];
                }
            }
        }
        
        p = p->next[2];
    }
    
    removeFromz(old);
    
}



// Decreasing version of the maximization of the HV of k points
/* Decremental greedy hypervolume subset selection in 3D */
static double gHSSD3D(dlnode_t * list, int n, int k, double * contribs, int * selected, int recompute){
    // double gHSSD3DTimes(double *data, int d, int n, int k, double *ref, double * contribs, int * selected, int recompute, double * times){
    //Timer_start();
    int ki;
    double hv = 0;
    int considerDominated = 0;    
    int adding = 0;
    dlnode_t * lcontributor;
    
    preprocessing(list);
    hv = hvc3d(list, considerDominated);

    for(ki = n-1; ki >= k; ki--){
        lcontributor = leastContributor(list, 3);

        selected[ki] = lcontributor->id;
        contribs[ki] = lcontributor->volume; // assumes d == 3

        hv -= lcontributor->volume;
    
        removeFromDataStructure(list, lcontributor);
        
        if(recompute){
            hv = hvc3d(list, considerDominated);
        }else{
            updateContributions(list, lcontributor, adding);
        }
        
        removeFromHistory(lcontributor);
        //times[ki] = Timer_elapsed_virtual();
    }
    
    dlnode_t * p = list->prev[0];
    for(ki = k-1; ki >= 0; ki--){
        selected[ki] = p->id;
        contribs[ki] = p->volume;
        p = p->prev[0];
        //times[ki] = 0;
    }
    
    return hv;
    
}


/* Compute the hypervolume indicator in d=3,4
 *      data    - array with all points
 *      d       - number of dimensions (either 3 or 4)
 *      n       - number of points
 *      ref     - reference point
 *      recompute - available for d=3 (and )
 *                  if 0: recompute the hypervolume contributions in 3 (uses hvc3d+-R)
 *                  otherwise: update only the hypervolume contributions that change in 3 (uses hvc3d+-R)
 * 
 * TODO: Implement gHSSD4D
 */
double gHSSD(double *data, int d, int n, int k, double *ref, double * contribs, int * selected, int recompute)
{
    //printf("gHSSD\n");
    dlnode_t * list = setup_cdllist(data, n, n, d, ref);
    double hv = 0;
    
    if(d == 3){
        
        hv = gHSSD3D(list, n, k, contribs, selected, recompute);
        
    }else{
        printf("Not supported yet\n");
    }
    free_cdllist(list);
    
    return hv;
}






