 
/*************************************************************************

 hv-class.c

 ---------------------------------------------------------------------

                        Copyright (c) 2017
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



#include "hvc-class.h"
#include "hvc-private.h"
#include "io.h"

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>



struct hvcstruct {
    dlnode_t * list;
    double * contribs;
    double hv;

    double * points;
    double * ref;
    int n;
    int naloc;
    int d;
    
    dlnode_t ** i2struct; //to be removed
    
    int * freeSpaces;
    int * freeIds;
    
    int leastContributorix;
    
    int updated; //are contributions updated?
    int hvupdated;

    int saved; //are the contributions saved in contribs vector?

    int ndom; // not fully implemented yet
  
}; 

/************** static functions **********************************/

//TODO: Manage memory - double naloc
static void allocateMoreSpace(){
    warnprintf("NOT IMPLEMENTED YET!\n");
}



static dlnode_t * newPoint(hvc_s * hvcs, double * point){
    
    dlnode_t * p = hvcs->list + hvcs->freeSpaces[hvcs->n];
    p->id = hvcs->freeIds[hvcs->n];
        
    p = point2Struct(hvcs->list, p, point, hvcs->d);

    return p;
}

//assumes hvcs->updated == 1
static void saveContributions(hvc_s * hvcs){
    if(!hvcs->updated) updateAllContributions(hvcs);
    
    int i, j, d = hvcs->d;
    dlnode_t * list = hvcs->list;
    dlnode_t * p;
    
    double * points = hvcs->points;
    double * contribs = hvcs->contribs;
    dlnode_t ** i2struct = hvcs->i2struct;
    
    double minContr = DBL_MAX;
    int leastContrix = -1;
    
    p = list->next[0];
    for(i = 0; p != list; i++, p = p->next[0]){
        if(hvcs->d == 3){
            contribs[i] = p->volume;
        }else{ // if(hvcs->d == 4){
            contribs[i] = p->hvolume;
        }
        i2struct[i] = p;
        for(j = 0; j < d; j++) points[d * i +j] = p->x[j];
        if(contribs[i] < minContr){
            minContr = contribs[i];
            leastContrix = i;
        }
    }
    hvcs->leastContributorix = leastContrix;
    hvcs->saved = 1;

}




static int removePointNode(hvc_s * hvcs, dlnode_t * point, int updateContr){
    if(hvcs->ndom > 0) return -1;
    
    dlnode_t * list = hvcs->list;
    
    if(point->ndomr < 2){
        removeFromDataStructure(list, point);
        if(updateContr){
            if(!hvcs->updated){
                int considerDominated = 1;
                hvcs->hv = hvc3d(hvcs->list, considerDominated);
            }else{
//                 printf("fast update - remove\n");
                hvcs->hv -= updateContributions(hvcs->list, point, 0);
            }
            hvcs->updated = 1;
        }else{
            hvcs->updated = 0;
            //hvcs->saved = 0;
        }
        if(point->ndomr == 1) hvcs->ndom -= 1;
    }else{
        //in this case, the contribution should be zero...
        if(updateContr && !hvcs->updated){
            int considerDominated = 1;
            hvcs->hv = hvc3d(hvcs->list, considerDominated);
        }
        hvcs->ndom -= 1;
    }
    removeFromHistory(point);
    
    int id = point->id;
    hvcs->n -= 1;
    if(hvcs->n == 0) hvcs->hv = 0; //force zero, otherwise, it could be different from zero due to rounding problems 
    hvcs->freeSpaces[hvcs->n] = (point-list);
    hvcs->freeIds[hvcs->n] = id;
    point->id = -1; 
    hvcs->saved = 0;
    hvcs->leastContributorix = -1;
    return 1;
}

//as this function is called only by removeLeastContributor, then it is not necessary to save info in hvcs->leastContributorix
static dlnode_t * getLeastContributorNode(hvc_s * hvcs){
//     printf("updated: %d, saved: %d\n", hvcs->updated, hvcs->saved);
    if(!hvcs->updated) updateAllContributions(hvcs);
    if(hvcs->saved)
        return hvcs->i2struct[hvcs->leastContributorix];
        
//     int i;
    dlnode_t * list = hvcs->list;
    dlnode_t * p;
    
    
    double contr = -1, minContr = DBL_MAX;
    dlnode_t * leastContr = NULL;
    
    for(p = list->next[0]; p != list; p = p->next[0]){
        if(hvcs->d == 3){
            contr = p->volume;
        }else{
            contr = p->hvolume;
        }
        if(contr < minContr){
            minContr = contr;
            leastContr = p;
        }
    }
//     printf("getLeastContributorNode id: %d (%f)\n", leastContr->id, minContr);
    return leastContr;
    
}



static int isequal(double * p, double * q, int d){
    int i;
    for(i = 0; i < d; i++){
        if(p[i] != q[i])
            return 0;
    }
    return 1;
}



/************** init and free functions ***************************/


hvc_s * init(double * data, int d, int n, int naloc, double *ref){

    int i;
    
    if(d != 3){
        warnprintf("ERROR!! d = %d not supported! Only d=3 is supported!\n", d);
        return NULL;
    }
    
    naloc = (n > naloc) ? n : naloc;
    
    double * newdata = (double *) malloc(d * naloc * sizeof(double));
    for(i = 0; i < d*n; i++) newdata[i] = data[i];
    double * contribs = (double *) malloc(naloc * sizeof(double));
    for(i = 0; i < n; i++){
        contribs[i] = 0;
    }
    double * newRef = (double *) malloc(d * sizeof(double));
    for(i = 0; i < d; i++) newRef[i] = ref[i];
    
    
    hvc_s * hvcs = (hvc_s *) malloc(sizeof(hvc_s));
    hvcs->list = setup_cdllist(newdata, naloc, n, d, newRef);

    hvcs->points = newdata; 
    hvcs->ref = newRef;
    hvcs->contribs = contribs;
    hvcs->hv = 0;
    hvcs->naloc = naloc;
    hvcs->n = n;
    hvcs->ndom = 0;
    
    
    hvcs->d = d;
    hvcs->updated = (n > 0) ? 0 : 1;
    hvcs->hvupdated = 0; 
    hvcs->saved = (n > 0) ? 0 : 1;
    hvcs->leastContributorix = -1;
    
    hvcs->i2struct = (dlnode_t **) malloc((naloc) * sizeof(dlnode_t *)); 
    
    int * freeSpaces = (int *) malloc((naloc) * sizeof(int));
    int * freeIds = (int *) malloc((naloc) * sizeof(int));
    for(i = 0; i < naloc; i++){
        freeSpaces[i] = i + 3; //skip sentinels
        freeIds[i] = i;
    }
    hvcs->freeSpaces = freeSpaces;
    hvcs->freeIds = freeIds;
    
    //setup data strucuture to compute contributions - O(n log n)
    if(n > 0){
        int ndom = preprocessing(hvcs->list);
        hvcs->ndom = ndom;
    }
    
    
    
    return hvcs;
}



double dealloc(hvc_s * hvcs){
    free_cdllist(hvcs->list);
    free(hvcs->contribs);    
    
    free(hvcs->i2struct);
    free(hvcs->freeSpaces);
    free(hvcs->freeIds);
    free(hvcs->points);
    free(hvcs->ref);
    
    free(hvcs);
    return 0;
}


double totalHV(hvc_s * hvcs){
    if(!hvcs->updated)
        updateAllContributions(hvcs);
    return hvcs->hv;
    
}

double * getContributions(hvc_s * hvcs){
    //printf("getContributions (saved: %d, updated: %d)\n", hvcs->saved, hvcs->updated);
    if(hvcs->saved) return hvcs->contribs;
    if(!hvcs->updated) updateAllContributions(hvcs);
    
    //save contributions in hvcs->contribs
    saveContributions(hvcs);
    
    
    return hvcs->contribs;
}

void updateAllContributions(hvc_s * hvcs){
    //printf("updateAllContributions\n");
    if(!hvcs->updated){
        int considerDominated = 1; 
        hvcs->hv = hvc3d(hvcs->list, considerDominated);
        
        hvcs->updated = 1;
        hvcs->saved = 0;
    }
}



//getters
int getSize(hvc_s * hvcs){
    return hvcs->n;
}

int getAlocsize(hvc_s * hvcs){
    return hvcs->naloc;
}

double * getPoints(hvc_s * hvcs){
    if(!hvcs->saved) saveContributions(hvcs);
    return hvcs->points;
}


double addPoint(hvc_s * hvcs, double * point, int updateContribs){
    if(hvcs->n == hvcs->naloc){
        allocateMoreSpace(hvcs);
    }
    
    dlnode_t * newp = newPoint(hvcs, point);
    
    if(!updateContribs){
        addToDataStructure(hvcs->list, newp, 1);
        hvcs->updated = 0;
        
    }else if(hvcs->updated){
        hvcs->hv += updateContributions(hvcs->list, newp, 1);
        addToDataStructure(hvcs->list, newp, 0);
        hvcs->updated = 1;
        
    }else{ //contributions are not updated and do not update them
        addToDataStructure(hvcs->list, newp, 1);

        updateAllContributions(hvcs);
        
    }
    
    hvcs->saved = 0;
    addToHistory(hvcs->list, newp);
    hvcs->n++;
    
    hvcs->leastContributorix = -1;
    
    return newp->volume;
}

int removePointAt(hvc_s * hvcs, int i, int updateContribs){
    if(hvcs->ndom > 0) return -1;
    if(!hvcs->saved) saveContributions(hvcs);
    removePointNode(hvcs, hvcs->i2struct[i], updateContribs);
    return 1;
}

//if there are repeated points, remove the oldest of the most dominated ones
int removePoint(hvc_s * hvcs, double * point, int updateContribs){
    
    dlnode_t * epoint = NULL;
    int epointndom = -1;
    int d = hvcs->d;
    dlnode_t * list = hvcs->list;
    dlnode_t * p = list->next[0];
    int ndom = hvcs->ndom;
    
    while(p != list && (epoint == NULL || ndom > 0)){                
        if(isequal(point, p->x, d) && p->ndomr > epointndom){
            epoint = p;
            epointndom = p->ndomr;
        }
        p = p->next[0];
    }
    
    if(epoint == NULL) return 0;
        
    removePointNode(hvcs, epoint, updateContribs);
    return 1;
    
}

void removeLeastContributor(hvc_s * hvcs, int updateContribs){
    dlnode_t * leastContributor = getLeastContributorNode(hvcs);
    removePointNode(hvcs, leastContributor, updateContribs);
}


int getLeastContributorIndex(hvc_s * hvcs){
    if(!hvcs->saved) saveContributions(hvcs);
    return hvcs->leastContributorix;
}


//private - just for gHSSD, to get the original index
int getLeastContributorId(hvc_s * hvcs){
    if(!hvcs->saved) saveContributions(hvcs);
    return hvcs->i2struct[hvcs->leastContributorix]->id;
}



double getLeastContributor(hvc_s * hvcs, double * point){
    if(!hvcs->saved) saveContributions(hvcs);
    int i, d = hvcs->d, ix = hvcs->leastContributorix;
    dlnode_t * p = hvcs->i2struct[ix];
    for(i = 0; i < d; i++){
        point[i] = p->x[i];
    }

    return hvcs->contribs[ix];
}

double getLeastContribution(hvc_s * hvcs){
    if(!hvcs->saved) saveContributions(hvcs);
    return hvcs->contribs[hvcs->leastContributorix];
}



void forceRecomputationNext(hvc_s * hvcs){
    hvcs->updated = 0;
}

int isUpToDate(hvc_s * hvcs){
    return hvcs->updated;
}


void printIds(hvc_s * hvcs){
    dlnode_t * p = hvcs->list->next[0];
    
    while(p != hvcs->list){
        printf("%d\n", p->id);
        p = p->next[0];
    }
    
}

//computes the contribution of point
double oneContribution(hvc_s * hvcs, double * point){
     if(hvcs->n == hvcs->naloc){
        allocateMoreSpace(hvcs);
    }
    
    dlnode_t * newp = newPoint(hvcs, point); //temporary point (it is not added to the data structure)
    double contr = oneContribution3d(hvcs->list, newp);
    return contr;
}


//adds the new point and computes only its contribution
double addOneContribution(hvc_s * hvcs, double * point){
     if(hvcs->n == hvcs->naloc){
        allocateMoreSpace(hvcs);
    }
    
    dlnode_t * newp = newPoint(hvcs, point);
    double contr = oneContribution3d(hvcs->list, newp);
    addToDataStructure(hvcs->list, newp, 0);
    
    hvcs->saved = 0;
    addToHistory(hvcs->list, newp);
    hvcs->n++;
    hvcs->leastContributorix = -1;
    
    return contr;
}




//TODO: keep track of hvupdated, set it correctly (look where hvcs->updated is being updated. Also, change totalHV() to check whether hvupdated == 1)
//updates the hypervolume only and not the contributions
double updateHypervolume(hvc_s * hvcs){
//     if(hvcs->updated || hvcs->hvupdated) return hvcs->hv;
    
    double hv = hv3dplus(hvcs->list);
    hvcs->hv = hv;
    hvcs->hvupdated = 1;
    return hv;
    
}


//TODO: check that the reference point was properly set
void changeReferencePoint(hvc_s * hvcs, double * ref){
    double * newRef = (double *) malloc(hvcs->d * sizeof(double));
    int i;
    for(i = 0; i < hvcs->d; i++) newRef[i] = ref[i];
    free(hvcs->ref);
    hvcs->ref = newRef;
    hvcs->updated = 0;
    hvcs->hvupdated = 0;
    hvcs->saved = 0;
    
    
}
