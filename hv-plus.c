 
/*************************************************************************

 hv-plus.c

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



#include "avl.h"
#include "timer.h"
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


/* ---------------------------------- Auxiliar Functions ----------------------------------*/
//3D points
static inline int lexicographicLess(double * a, double * b){
    return (a[2] < b[2] || (a[2] == b[2] && (a[1] < b[1] || (a[1] == b[1] && a[0] <= b[0]))));
}




/* ---------------------------------- Data Structure ---------------------------------------*/

typedef struct dlnode {
  double x[4];                    /* The data vector              */
  struct dlnode * closest[2]; // closest[0] == cx, closest[1] == cy
  struct dlnode * cnext[2]; //current next

  struct dlnode * next[4]; //keeps the points sorted according to coordinates 2,3 and 4
                           // (in the case of 2 and 3, only the points swept by 4 are kept)
  struct dlnode *prev[4]; //keeps the points sorted according to coordinates 2 and 3 (except the sentinel 3)
  
  int ndomr;    //number of dominators
} dlnode_t;




/* ---------------------------------- Data Structures Functions ---------------------------------------*/


static dlnode_t * initSentinels(dlnode_t * list, const double * ref, int d){
 
    dlnode_t * s1 = list;
    dlnode_t * s2 = list + 1;
    dlnode_t * s3 = list + 2;
    
    s1->x[0] = -DBL_MAX;
    s1->x[1] = ref[1];
    s1->x[2] = -DBL_MAX;
    s1->x[3] = -DBL_MAX;
    s1->closest[0] = s2;
    s1->closest[1] = s1;  

    s1->next[2] = s2;
    s1->next[3] = s2;
    s1->cnext[1] = NULL;  
    s1->cnext[0] = NULL; 
    
    s1->prev[2] = s3;
    s1->prev[3] = s3;
    s1->ndomr = 0;

    
    s2->x[0] = ref[0];
    s2->x[1] = -DBL_MAX;
    s2->x[2] = -DBL_MAX;
    s2->x[3] = -DBL_MAX;
    s2->closest[0] = s2; 
    s2->closest[1] = s1; 

    s2->next[2] = s3;
    s2->next[3] = s3;
    s2->cnext[1] = NULL;  
    s2->cnext[0] = NULL;  
    
    s2->prev[2] = s1;
    s2->prev[3] = s1;
    s2->ndomr = 0;

    
    
    s3->x[0] = -INT_MAX; 
    s3->x[1] = -INT_MAX; 
    s3->x[2] = ref[2];
    if(d == 4)
        s3->x[3] = ref[3];
    else
        s3->x[3] = - DBL_MAX;
    s3->closest[0] = s2;
    s3->closest[1] = s1;
    
    s3->next[2] = s1;
    s3->next[3] = NULL;
    s3->cnext[1] = NULL;  
    s3->cnext[0] = NULL;  
    
    s3->prev[2] = s2;
    s3->prev[3] = s2;
    s3->ndomr = 0;

    
    return s1;
    
}





static void clearPoint(dlnode_t * list, dlnode_t * p){
    
    p->closest[1] = list;
    p->closest[0] = list->next[2]; 
    
    /* because of printfs */
    p->cnext[1] = list;
    p->cnext[0] = list->next[2];
    
    
    p->ndomr = 0;
    
}



static dlnode_t * point2Struct(dlnode_t * list, dlnode_t * p, double * v, int d){
    
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


/* check if new is dominated, find cx and cy of the 'new' point and find where to insert 'new' in the
 * list sorted by z
 */
static void setupZandClosest(dlnode_t * list, dlnode_t * new){
    
            
    double * closest1 = (double *) (list);
    double * closest0 = (double *) (list->next[2]);

    dlnode_t * q = (list->next[2]->next[2]);
    
    double * newx = new->x;
    
    
    while(lexicographicLess(q->x, newx)){
        if(q->x[0] <= newx[0] && q->x[1] <= newx[1]){
                
            new->ndomr += 1;
            //new->domr = q;
            //return new;
                
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




static int updateLinks(dlnode_t * list, dlnode_t * new, dlnode_t * p){
    
    dlnode_t * stop = list->prev[2];

    int ndom = 0;
    int allDelmtrsVisited = 0;
//     while(p != stop){
    while(p != stop && !allDelmtrsVisited){
//         q = p->next[2];
        
        
        if(p->x[0] <= new->x[0] && p->x[1] <= new->x[1] && (p->x[0] < new->x[0] || p->x[1] < new->x[1])){
            
            allDelmtrsVisited = 1;
            
        }else {
            
            if(new->x[0] <= p->x[0]){
                //new <= p
                if(new->x[1] <= p->x[1]){
                    p->ndomr++;
//                     p->domr = new;
                    ndom += 1;
                    removeFromz(p); //HV-ONLY (does not need dominated to compute HV)
                    
                }else if(new->x[0] < p->x[0] && (new->x[1] < p->closest[1]->x[1] || (new->x[1] == p->closest[1]->x[1] && (new->x[0] < p->closest[1]->x[0] || (new->x[0] == p->closest[1]->x[0] && new->x[2] < p->closest[1]->x[2]))))){ // new->x[1] > p->x[1]
                    p->closest[1] = new;
                }
            }else if(new->x[1] < p->x[1] && (new->x[0] < p->closest[0]->x[0] || (new->x[0] == p->closest[0]->x[0] && (new->x[1] < p->closest[0]->x[1] || (new->x[1] == p->closest[0]->x[1] && new->x[2] < p->closest[0]->x[2]))))){//new->x[0] > p->x[0]
                p->closest[0] = new;
            }
            
            
//             if(p->ndomr > 1){
//                 printf("remove\n");
//                 removeFromz(p);
//             }
        }
//         p = q;
    p = p->next[2];
    }
    
    return ndom;
}




/* ---------------------------------- Sort ---------------------------------------*/

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
static dlnode_t *
setup_cdllist(double * data, int naloc, int n, int d, const double *ref)
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
    //             printf("%f\n", scratchd[i][2]);
            for(j = 0; j < d; j++){
                data2[d * i + j] = scratchd[i][j];
            }
        }
        
        
    //     int d = 3;
        dlnode_t ** scratch = (dlnode_t **) malloc(n * sizeof(dlnode_t *));

        for (i = 0; i < n; i++) {
            scratch[i] = point2Struct(list, head+i+3, &data2[i*d], d);
    //             scratch[i] = newPoint(head, &data[i*d], d);   
//             scratch[i]->id = order[i];
//             scratch[i]->id = (scratchd[i]-data)/3;
        }

        
        free(scratchd);
        
        
        dlnode_t * s = head->next[di];
        s->next[di] = scratch[0];
        scratch[0]->prev[di] = s;

                
        for(i = 0; i < n-1; i++){
            scratch[i]->next[di] = scratch[i+1];
            scratch[i+1]->prev[di] = scratch[i];
        }
        
        s = head->prev[di];
        s->prev[di] = scratch[n-1];
        scratch[n-1]->next[di] = s;
        
        free(scratch);
        free(data2);
    }
    
    return head;
}



static void free_cdllist(dlnode_t * list)
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



static void preprocessing(dlnode_t * list){


//     restartListy(list);
    
    avl_tree_t * avltree = avl_alloc_tree ((avl_compare_t) compare_tree_asc_y, NULL);
    
    dlnode_t * p = list;
    

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
            p->ndomr = 1;
            free(node);
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
    
}




/* ----------------------Hypervolume Indicator Algorithms ---------------------------------------*/



static void restartListy(dlnode_t * list){
    
    list->next[2]->cnext[1] = list; //link sentinels sentinels ((-inf ref[1] -inf) and (ref[0] -inf -inf))
    list->cnext[0] = list->next[2];
    
}

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

//does what setupZandClosest does while reconstructing L at z = new->x[2]
static void restartBaseSetupZandClosest(dlnode_t * list, dlnode_t * new){
    

    dlnode_t * p = list->next[2]->next[2];
    double * closest1 = (double *) (list);
    double * closest0 = (double *) (list->next[2]);

    double * newx = new->x;
    
    restartListy(list);
    
    while(lexicographicLess(p->x, newx)){
        
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
    
        p = p->next[2];
    }
    
    
    new->closest[0] = (dlnode_t *) closest0;
    new->closest[1] = (dlnode_t *) closest1;
    
    new->prev[2] = p->prev[2];
    new->next[2] = p;
    
    
}

static double oneContribution3d(dlnode_t * list, dlnode_t * new){
    
//     int considerDominated = 1;
    
    dlnode_t * p = list;
    double area = 0;
    double volume = 0;
    double x[3];
    
    restartBaseSetupZandClosest(list, new);
    if (new->ndomr > 0)
        return 0;
    
    new->cnext[0] = new->closest[0];
    new->cnext[1] = new->closest[1];
    area = computeAreaSimple(new->x, 1, new->cnext[0], new->cnext[0]->cnext[1]);
    
    p = new->next[2];
    double lastz = new->x[2];
    
    while(p->x[0] > new->x[0] || p->x[1] > new->x[1]){
        volume += area * (p->x[2]- lastz);
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            
            if(p->x[0] >= new->x[0] && p->x[1] >= new->x[1]){
                area -= computeAreaSimple(p->x, 1, p->cnext[0], p->cnext[0]->cnext[1]);
                p->cnext[1]->cnext[0] = p;
                p->cnext[0]->cnext[1] = p;
                
            }else if(p->x[0] >= new->x[0]){
                if(p->x[0] <= new->cnext[0]->x[0]){
                    x[0] = p->x[0]; x[1] = new->x[1]; x[2] = p->x[2]; 
                    area -= computeAreaSimple(x, 1, new->cnext[0], new->cnext[0]->cnext[1]);
                    p->cnext[0] = new->cnext[0];
                    p->cnext[1]->cnext[0] = p;
                    new->cnext[0] = p;
                }
            }else{
                if(p->x[1] <= new->cnext[1]->x[1]){
                    x[0] = new->x[0]; x[1] = p->x[1]; x[2] = p->x[2]; 
                    area -= computeAreaSimple(x, 0, new->cnext[1], new->cnext[1]->cnext[0]);
                    p->cnext[1] = new->cnext[1];
                    p->cnext[0]->cnext[1] = p;
                    new->cnext[1] = p;
                }
                
            }
        lastz = p->x[2];
        p = p->next[2];
        
    }
    
    volume += area * (p->x[2]- lastz);
    return volume;
    
}



static double hv3dplus(dlnode_t * list){
    
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
        }else{
            removeFromz(p);
        }
        
        volume += area * (p->next[2]->x[2]- p->x[2]);
        
        p = p->next[2];
    }
    
    
    return volume;
    
}


/* Compute the hypervolume indicator in d=4 by iteratively
 * computing the hypervolume indicator in d=3 (using hv3d+)
 */
double hv4dplusR(dlnode_t * list)
{
    double height = 0, volume = 0, hv = 0;
    
    dlnode_t * last = list->prev[3];
    dlnode_t * new = list->next[3]->next[3];
    
    while(new != last){
        
        setupZandClosest(list, new);          // compute cx and cy of 'new' and determine next and prev in z
        addToZ(new);                          // add to list sorted by z
        updateLinks(list, new, new->next[2]); // update update cx and cy of the points above 'new' in z
                                              // and removes dominated points

        volume = hv3dplus(list);              // compute hv indicator in d=3 in linear time 

        height = new->next[3]->x[3] - new->x[3];
        hv += volume * height;                // update hypervolume in d=4
    
        new = new->next[3];
    }
        
    return hv;
}




/* Compute the hypervolume indicator in d=4 by iteratively
 * computing the one contribution problem in d=3
 */
double hv4dplusU(dlnode_t * list)
{

    double height = 0, volume = 0, hv = 0;
    
    dlnode_t * last = list->prev[3];
    dlnode_t * new = list->next[3]->next[3];
    
    
    while(new != last){
    
        volume += oneContribution3d(list, new);
        addToZ(new);
        updateLinks(list, new, new->next[2]);

        height = new->next[3]->x[3] - new->x[3];
        hv += volume * height;
    
        new = new->next[3];
    }
        

    return hv;
    
}


/* Compute the hypervolume indicator in d=3,4
 *      data    - array with all points
 *      d       - number of dimensions (either 3 or 4)
 *      n       - number of points
 *      ref     - reference point
 *      recompute - available for d=4 only
 *                  if 0: recompute the hypervolume in d-1 (hv4d+-R)
 *                  otherwise: compute one contribution in d-1 (hv4d+-U)
 */
double hvplus(double *data, int d, int n, double *ref, int recompute)
{
    double hv = 0;
    dlnode_t * list = setup_cdllist(data, n, n, d, ref);
    if(d == 3){
        
        preprocessing(list);
        hv = hv3dplus(list);
        
    }else{
        if(recompute)
            hv = hv4dplusR(list);
        else
            hv = hv4dplusU(list);
    }
    
    free_cdllist(list);
    
    return hv;
}

