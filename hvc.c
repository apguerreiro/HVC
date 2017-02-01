 
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




/* ---------------------------------- Data Structure ---------------------------------------*/

//TODO: colocar volume e hvolume num array
//TODO: corrigir estrutura (nao e preciso a 4a dimensao)
typedef struct dlnode {
  double x[4];                    /* The data vector              */
  struct dlnode * closest[2]; // closest[0] == cx, closest[1] == cy
  struct dlnode * cnext[2]; //current next
  struct dlnode * head[2]; //lowest (0 - x, 1 - y)
  
  double area;
  double volume;
  double lastSlicez;
  
  struct dlnode * next[4]; //keeps the points sorted according to coordinates 2,3 and 4 (in the case of 2 and 3, only the points swept by 4 are kept)
  struct dlnode *prev[4]; //keeps the points sorted according to coordinates 2 and 3 (excepto o sentinela 3)
  
  double hvolume; //4D
  int ndomr;    //number of dominators
  struct dlnode * domr; //dominator
  int id;
  
  double oldvolume; //HVC-ONLY 4D-U-ADD
//   double volumecp;
//   double xcp[3];
//   struct dlnode * closestcp[2];
//   int ndomrcp;
//   struct dlnode * domrcp;
  
  
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
//     s1->head[0] = s1->head[1] = s1;
    s1->area = 0; //HVC-ONLY
    s1->volume = 0; //HVC-ONLY
    s1->head[0] = s1->head[1] = s1; //HVC-ONLY
    s1->oldvolume = 0; //HVC-ONLY 4D-U-ADD
    
//     s1->prev[1] = s3;
//     s1->next[1] = NULL;
    s1->next[2] = s2;
    s1->next[3] = s2;
    s1->cnext[1] = NULL;  //para nao dar erro ao imprimir
    s1->cnext[0] = NULL;  //para nao dar erro ao imprimir
    
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
// // //     s2->head[0] = s2->head[1] = NULL;
    s2->head[0] = s2->head[1] = s2; //HVC-ONLY
    
//     s2->prev[1] = NULL;
// //     s2->next[1] = s1;
//     s2->next[1] = s3;
    s2->next[2] = s3;
    s2->next[3] = s3;
    s2->cnext[1] = NULL;  
    s2->cnext[0] = NULL;  
    
    s2->prev[2] = s1;
    s2->prev[3] = s1;
    s2->ndomr = 0;
    s2->domr = NULL; //HVC-ONLY
    s2->id = -2; //HVC-ONLY
    
// //     s2->cnext[1] = s1; //ligar as duas sentinelas (p: -inf ref[1] -inf e p->next[2]: ref[0] -inf -inf)
// //     s1->cnext[0] = s2;
    
    
    s3->x[0] = -INT_MAX; //minv - 10;
    s3->x[1] = -INT_MAX; //minv - 10;
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
// // //     s3->head[0] = s3->head[1] = NULL;
    s3->head[0] = s3->head[1] = s3; //HVC-ONLY
    
//     s3->next[1] = s1;
//     s3->prev[1] = s2;
// //     s3->next[1] = NULL;
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
    
    /* Esta parte e so por causa dos printfs */
    p->cnext[1] = list;
    p->cnext[0] = list->next[2];
    
// //     p->head[0] = p->head[1] = NULL;   
    p->head[0] = p->cnext[0]; //HVC-ONLY
    p->head[1] = p->cnext[1]; //HVC-ONLY
    
    p->area = 0; //HVC-ONLY
    p->volume = 0; //HVC-ONLY
    p->oldvolume = 0; //HVC-ONLY 4D-U-ADD
    p->hvolume = 0; //HVC-ONLY 4D
    p->lastSlicez = p->x[2]; //HVC-ONLY
    p->ndomr = 0;
    p->domr = NULL; //HVC-ONLY
    
}



static dlnode_t * point2Struct(dlnode_t * list, dlnode_t * p, double * v, int d){
    
//     int d = 3;
    
    int i;
    for(i = 0; i < d; i++)
        p->x[i] = v[i];
    
    
//     p->next[1] = NULL;
//     p->prev[1] = NULL;
    
    
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
    
            
//     dlnode_t * closest0 = list->next[2];
//     dlnode_t * closest1 = list;
//     new->closest[0] = closest0;
//     new->closest[1] = closest1;

//     double * closest1 = (double *) (new->closest[1]);
//     double * closest0 = (double *) (new->closest[0]);
    
    double * closest1 = (double *) (list);
    double * closest0 = (double *) (list->next[2]);

    dlnode_t * q = (list->next[2]->next[2]);
    
    double * newx = new->x;
    
    
    while(q->x[2] <= newx[2]){
//         printf("%f %f %f\n", q->x[0], q->x[1], q->x[2]);
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




/* Update cx and cy of the points above 'new' in z and remove dominated points
 */
static int updateLinks(dlnode_t * list, dlnode_t * new, dlnode_t * p){
    
    dlnode_t * stop = list->prev[2];

    int ndom = 0;
    int allDelmtrsVisited = 0;
//     while(p != stop){
    while(p != stop && allDelmtrsVisited < 2){ //HVC-ONLY 4D
//         q = p->next[2];
        
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
    //                     removeFromz(p); //HV-ONLY 
                        
                    }else if(new->x[0] < p->x[0] && (new->x[1] < p->closest[1]->x[1] || (new->x[1] == p->closest[1]->x[1] && (new->x[0] < p->closest[1]->x[0] || (new->x[0] == p->closest[1]->x[0] && new->x[2] < p->closest[1]->x[2]))))){ // new->x[1] > p->x[1] //testeDE.6
                        p->closest[1] = new;
                    }
                }else if(new->x[1] < p->x[1] && (new->x[0] < p->closest[0]->x[0] || (new->x[0] == p->closest[0]->x[0] && (new->x[1] < p->closest[0]->x[1] || (new->x[1] == p->closest[0]->x[1] && new->x[2] < p->closest[0]->x[2]))))){//new->x[0] > p->x[0] //testeDE.6
                    p->closest[0] = new;
                }
                
                
                if(p->ndomr > 1){ //HVC-ONLY 4D
    //                 printf("remove\n"); //HVC-ONLY 4D
                    removeFromz(p); //HVC-ONLY 4D
                } //HVC-ONLY 4D
            }
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
            scratch[i]->id = (scratchd[i]-data)/d;
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

//     if (x1 != x2)
//         return (x1 < x2) ? -1 : 1;
//     else
//         return 0;
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
//             p->ndomr = 2;
            p->ndomr = 1;
            p->domr = list; 
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
    
    list->next[2]->cnext[1] = list; 
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






static void setupNDPoint(dlnode_t * p){
    
    
    p->cnext[0] = p->closest[0];
    p->cnext[1] = p->closest[1];
    
    p->head[1] = p->cnext[0]->cnext[1];
    p->head[0] = p->cnext[1]->cnext[0];
    
    
}



static void setupDomPoint(dlnode_t * p){
    
#if VARIANT == 2
    printf("setupDomPoint id: %d\n", p->id);
#endif
    
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
    
#if VARIANT == 2
    printf("SP0 - head1: %d\n", head1->id);
    printf("SP0 - head0: %d\n", head0->id);
#endif
    
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
#if VARIANT == 2
    printf("addPoint: %f %f %f | %d\n", p->x[0], p->x[1], p->x[2], p->id);
#endif
    
    
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

//     if(p->cnext[0]->cnext[1]->x[1] > p->x[1] || (p->cnext[0]->cnext[1]->x[1] == p->x[1] && p->cnext[0]->cnext[1]->x[0] > p->x[0]))
    if(p->cnext[0]->cnext[1] == p->domr)
        p->domr->head[1] = p;
    else
        p->cnext[0]->cnext[1] = p;

//     if(p->cnext[1]->cnext[0]->x[0] > p->x[0] || (p->cnext[1]->cnext[0]->x[0] == p->x[0] && p->cnext[1]->cnext[0]->x[1] > p->x[1]))
    if(p->cnext[1]->cnext[0] == p->domr)
        p->domr->head[0] = p;
    else
        p->cnext[1]->cnext[0] = p;   
    
     
}




static double hvc3d(dlnode_t * list, int considerDominated){
    
//     int considerDominated = 1;
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
//         p->oldvolume = 0; 
        p->lastSlicez = p->x[2];
//         volume += area * (p->x[2]- p->prev[2]->x[2]);
//         printf("hvc3d p: %f %f %f | %d\n", p->x[0], p->x[1], p->x[2], p->ndomr);
    
        if(p->ndomr < 1){
//         if(p->ndomr < 2){

            setupNDPoint(p);
            
            updateVolumeSimple(p->x, 1, p->head[1]);
            p->area = computeAreaSimple(p->x, 1, p->cnext[0], p->head[1]);
            area += p->area;

            q = p->cnext[0]; 
            x[0] = q->x[0]; x[1] = p->x[1]; x[2] = p->x[2]; //o x[2] nao tem importancia (join(p,q))
            q->area -= computeAreaSimple(x, 0, p->head[1], q->head[0]);
    
            q = p->cnext[1]; 
            x[0] = p->x[0]; x[1] = q->x[1];
            q->area -= computeAreaSimple(x, 1, p->head[0], q->head[1]);
    
            addNDPoint(p);
    
        }else if(considerDominated && p->ndomr == 1){
//             p->volume = 0;
            
            updateVolume(p->domr, p->x[2]);
            setupDomPoint(p);
            p->domr->area -= computeAreaSimple(p->x, 1, p->cnext[0], p->head[1]);
//             addNDPoint(p);
            addDomPoint(p);
        }else{
            //too dominated
        }
        
        volume += area * (p->next[2]->x[2]- p->x[2]);
//         printf("volume: %f\n", volume);
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
        
        setupZandClosest(list, new);          // compute cx and cy of 'new' and determine next and prev in z
        addToZ(new);                          // add to list sorted by z
        updateLinks(list, new, new->next[2]); // update update cx and cy of the points above 'new' in z
                                              // and removes dominated points
        
        volume = hvc3d(list, considerDominated);  // compute hv indicator in d=3 in linear time 
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






static void saveContribs(dlnode_t * list, double * contribs, int d){
    
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
double hvc(double *data, int d, int n, double *ref, double * contribs, int recompute){
//     printf("hvc\n");
    double hv = 0;
    dlnode_t * list = setup_cdllist(data, n, n, d, ref);
    
    if(d == 2){
//         hv = hvc2d(data, n, ref, contribs, stdout);
    }else if(d == 3){
        
        preprocessing(list);
        int considerDominated = 0;
        hv = hvc3d(list, considerDominated);
        
        
    }else if(d == 4){
//         hv = fpli_hvc4d2(data, n, ref, contribs, stdout);
        if(recompute)
            hv = hvc4dR(list);
//         else
//             hv = hvc4dU(list);
    }
    
    saveContribs(list, contribs, d);
    
    free_cdllist(list);
    
    return hv;
}



static dlnode_t * leastContributor(dlnode_t * list, int d){
    dlnode_t * p = list->next[d-1]->next[d-1];
    dlnode_t * stop = list->prev[d-1];
    
    double leastcontribution = DBL_MAX;
    dlnode_t * lcontributor = NULL;
    
    while(p != stop){
        if(p->volume < leastcontribution){ // ->volume assume d == 3
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

//reconstructing L at z = old->x[2]
static void restartBase(dlnode_t * list, dlnode_t * old){
    

    dlnode_t * p = list->next[2]->next[2];
    restartListy(list);
    
    while(p != old){
            
            //reconstruct
            p->cnext[0] = p->closest[0];
            p->cnext[1] = p->closest[1];
            
            p->cnext[0]->cnext[1] = p;
            p->cnext[1]->cnext[0] = p;
            
        p = p->next[2];
    }
}



static void removeFromDataStructure(dlnode_t * list, dlnode_t * old){
    
    dlnode_t * stop = list->prev[2];
    dlnode_t * p;
    
    restartBase(list, old);
    p = old->next[2];
    
    
    while(p != stop){
        
        if(p->x[2] == old->x[2]){
            old->next[2] = p->next[2];
            old->prev[2] = p;
        }
        if(old->x[0] <= p->x[0] && old->x[1] <= p->x[1]){
            p->ndomr--;
            
        }else if(updateInsertingPoints(old, p, 0)){
            old->cnext[1] = p;
            if(p->closest[0]->x[1] >= old->x[1]){ 
                old->head[0] = p->closest[0];
//                 old->head[0]->cnext[1] = p;
            }else{
//                 old->head[0] = old->head[1] = NULL;
                old->head[0] = old->cnext[0];
                old->head[1] = old->cnext[1];
            }
        }else{
            if(updateInsertingPoints(old, p, 1)){

                old->cnext[0] = p;
                if(p->closest[1]->x[0] >= old->x[0]){ 
                    old->head[1] = p->closest[1];
//                     old->head[1]->cnext[0] = p;
                }else{
//                     old->head[0] = old->head[1] = NULL;
                    old->head[0] = old->cnext[0];
                    old->head[1] = old->cnext[1];
                }
            }
            
        }
    
//     old->head[0] = oldhead[0];
//     old->head[1] = oldhead[1];
//     old->cnext[0] = oldcnext[0];
//     old->cnext[1] = oldcnext[1];

        p = p->next[2];
    
    }
    
}

// Decreasing version of the maximization of the HV of k points
static double gHSSD3D(dlnode_t * list, int n, int k, double * contribs, int * selected, int recompute){
    
//     printf("gHSSD3D(%d,%d)\n", n, k);
    int ki;
    double hv = 0;
    dlnode_t * lcontributor;
    preprocessing(list);
    int considerDominated = 0;
    hv = hvc3d(list, considerDominated);
//     printf("hv: %f\n", hv);
    for(ki = n-1; ki >= k; ki--){
        
        lcontributor = leastContributor(list, 3);
        
        selected[ki] = lcontributor->id;
        contribs[ki] = lcontributor->volume; // assume d == 3
//         printf("leastContributor(%d): %f\n", lcontributor->id, lcontributor->volume);
        
        hv -= lcontributor->volume;
        
        if(recompute){
            removeFromDataStructure(list, lcontributor);
            removeFromz(lcontributor); 
            hv = hvc3d(list, considerDominated);
        }
    }
    dlnode_t * p = list->next[2]->next[2];
    for(ki = k-1; ki >= 0; ki--){
        selected[ki] = p->id;
        contribs[ki] = p->volume;
        p = p->next[2];
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
double gHSSD(double *data, int d, int n, int k, double *ref, double * contribs, int * selected, int recompute)
{
    if(d != 3){
        printf("Not supported yet\n");
        return 0;
    }
    
//     printf("gHSSD\n");
    dlnode_t * list = setup_cdllist(data, n, n, d, ref);
    double hv = 0;
    
    if(d == 3){
        
        hv = gHSSD3D(list, n, k, contribs, selected, recompute);
        
    }
    free_cdllist(list);
    
    return hv;
}


