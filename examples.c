 
/*************************************************************************

 examples.c

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

#include "timer.h"
#include "examples.h"
#include "hvc-class.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))


int _quiet = 0; //1 - quiet (just prints the result), 2 - extra quiet (do not print anything)
int _checkCorrectness;
double _time_elapsed_cpu;




void printPoints(double * points, int d, int n){
    printf("POINTS:\n");
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < d; j++){
            printf("%f ", points[i*d + j]);
        }
        printf("\n");
    }
}


void printContributions(double * contribs, int n){
    printf("CONTRIBUTIONS:\n");
    int i;
    for(i = 0; i < n; i++){
        printf("%f\n", contribs[i]);
    }
}

void printContributionsFromHVCS(hvc_s * hvcs){
    int sz = getSize(hvcs);
    
    // O(n) because contributions are copied to an array, although if called a
    // second time without any changes in the data structure, then O(1)
    double * contribs = getContributions(hvcs); 
    
    
    if(_quiet < 1) printf("\t\tcontributions (size: %d):\n", sz);
    int j;
    for(j = 0; j < sz; j++){
        if(_quiet < 1) printf("\t\t%d: ", j);
        printf("%-16.15g\n", contribs[j]);
    }
    if(_quiet < 1) printf("\t\tHV: %f\n", totalHV(hvcs));
    
}





void testInitClose(double * data, int d, int n, double * ref){
    if(_quiet < 1) printf("---- init and free data structures ----\n");
    
    if(_quiet < 1) printf("\tinit: %d points in %d dimensions, ref = (%f, %f, %f)\n", n, d, ref[0], ref[1], ref[2]);
    hvc_s * hvcs = init(data, d, n, n, ref);
    
    if(_quiet < 1) printf("\tfree\n");
    dealloc(hvcs);
}



//TODO: test with dominated points!
double computeHypervolumeOnce(double * data, int d, int n, double * ref){
    if(_quiet < 1) printf("--- compute hypervolume once ----\n");
    
    
    
    if(_quiet < 1) printf("\tinit: %d points in %d dimensions, ref = (%f, %f, %f)\n", n, d, ref[0], ref[1], ref[2]);
    hvc_s * hvcs = init(data, d, n, n, ref);
    
    
    if(_quiet < 1) printf("\tcompute hypervolume\n");
    
    /* fastest way */
    double totalhv = updateHypervolume(hvcs); //uses HV3D+
    
    /*alternative that uses HVC3D (slightly less eficient) */
    //double * contribs = getContributions(hvcs); //compute contributions (and hypervolume)
    //double totalhv = totalHV(hvcs); //O(1) because everything is up-to-date

    
    dealloc(hvcs);
    
    
    if(_quiet < 1) printf("Hypervolume: ");
    if(_quiet < 2) printf("%-16.15g\n", totalhv);
    
    

    return totalhv;
}


double computeContributionsOnce(double * data, int d, int n, double * ref){
    if(_quiet < 1) printf("--- compute contributions once ----\n");
    double totalhv = 0;
    
    
    if(_quiet < 1) printf("\tinit: %d points in %d dimensions, ref = (%f, %f, %f)\n", n, d, ref[0], ref[1], ref[2]);
    hvc_s * hvcs = init(data, d, n, n, ref);
    
    
    if(_quiet < 1) printf("\tcompute contributions\n");
    double * contribs = getContributions(hvcs); //although contributions were not explicitly computed, that is checked here and in that case, the contributions are computed before returning
    
    
    if(_quiet < 2){
        int j, sz = getSize(hvcs);
        contribs = getContributions(hvcs); //O(1), because no point was added nor removed after the first call to this function
        if(_quiet < 1) printf("\t\tcontributions (size: %d):\n", sz);
        for(j = 0; j < sz; j++){
            printf("%-16.15g\n", contribs[j]);
        }
        
        totalhv = totalHV(hvcs); //O(1) because everything is up-to-date
    }
    
    dealloc(hvcs);
    
    
    return totalhv;
}


//outputFormat - 0 (default), prints the hv after each point is added, 1 - prints the contributions at the end
//update contributions by recomputing only those that change (HVC3D-U)
double unboundedArchive(double * data, int d, int n, double * ref, int outputFormat){
    if(_quiet < 1) printf("--- compute unbounded archive ----\n");
    
    
    if(_quiet < 1) printf("\tinit empty archive (%d dimensions) and set ref = (%f, %f, %f)\n", d, ref[0], ref[1], ref[2]);
    
    hvc_s * hvcs = init(data, d, 0, n, ref);
    int i;
    double * point;
    int updateContributionsRightAway = 1;
    
    
    for(i = 0; i < n; i++){
        if(_quiet < 1 && outputFormat == 0) printf("hypervolume of archive: ");
        if(_quiet < 2 && outputFormat == 0) printf("%-16.15g\n", totalHV(hvcs));
        point = &data[d*i];
        if(_quiet < 1) printf("\tadd Point: %f %f %f\n", point[0], point[1], point[2]);
        
        addPoint(hvcs, point, updateContributionsRightAway); //O(n)
        if(_quiet < 1) printf("\tcontributions are up-to-date? %d\n", isUpToDate(hvcs));
        
    }
    
    double totalhv = totalHV(hvcs);
    if(_quiet < 2 && outputFormat == 1) printContributionsFromHVCS(hvcs);
    if(_quiet < 2 && outputFormat == 0) printf("%-16.15g\n", totalhv);
    
    dealloc(hvcs);
    
    return totalhv;
}


//outputFormat - 0 (default), prints the hv after each point is added, 1 - prints the contributions at the end
//update contributions by recomputing all contributions (HVC3D-R)
double unboundedArchiveSlower(double * data, int d, int n, double * ref, int outputFormat){
    if(_quiet < 1) printf("--- compute unbounded archive ----\n");
    
    
    if(_quiet < 1) printf("\tinit empty archive (%d dimensions) and set ref = (%f, %f, %f)\n", d, ref[0], ref[1], ref[2]);
    
    hvc_s * hvcs = init(data, d, 0, n, ref);
    int i;
    double * point;
    int updateContributionsRightAway = 0;
    
    
    for(i = 0; i < n; i++){
        if(_quiet < 1 && outputFormat == 0) printf("hypervolume of archive: ");
        if(_quiet < 2 && outputFormat == 0) printf("%-16.15g\n", totalHV(hvcs));
        point = &data[d*i];
        if(_quiet < 1) printf("\tadd Point: %f %f %f\n", point[0], point[1], point[2]);
        
        addPoint(hvcs, point, updateContributionsRightAway); //O(n) - just add to the data structure
        
        if(_quiet < 1) printf("\tcontributions are up-to-date? %d\n", isUpToDate(hvcs));
        updateAllContributions(hvcs); //O(i) - compute all contributions from scratch (HVC3D) if not up-to-date, otherwise, O(1)
        
    }
    if(_quiet < 2 && outputFormat == 1) printContributionsFromHVCS(hvcs);
    if(_quiet < 2 && outputFormat == 0) printf("%-16.15g\n", totalHV(hvcs));
    
    double totalhv = totalHV(hvcs);

    
    dealloc(hvcs);
    
    return totalhv;
}




//update contributions by recomputing only those that change (HVC3D-U)
double replicateGHSSD3D(double * data, int d, int n, double * ref, int k, int updateContributionsRightAway, int outputFormat){
    if(_quiet < 1) printf("--- compute gHSSD3D (size: %d, k: %d)----\n", n, k);
    
    int removedId[n];
    double removedContr[n];
    
    if(_quiet < 1) printf("\tcompute contribs (%d dimensions) and set ref = (%f, %f, %f)\n", d, ref[0], ref[1], ref[2]);
    
    hvc_s * hvcs = init(data, d, n, n, ref);
    int i, j, sz;
    double * contribs;
    
    
    for(i = n; i > k; i--){
        //point = &data[d*i];
        //printf("i: %d (k: %d)\n", i, k);
        
        
        if(_quiet < 1) printf("\tremove least contribution: %f\n", getLeastContribution(hvcs));
        
        removedId[i-1] = getLeastContributorId(hvcs);
        removedContr[i-1] = getLeastContribution(hvcs);

        //printf("least contribution: %f\n", getLeastContribution(hvcs));
        removeLeastContributor(hvcs, updateContributionsRightAway);
        
        if(_quiet < 1) printf("\tcontributions are up-to-date? %d\n", isUpToDate(hvcs));
        //if(_quiet < 1) printContributionsFromHVCS(hvcs);
        
        if(!updateContributionsRightAway) updateAllContributions(hvcs);
        
    }
    double totalhv = totalHV(hvcs);
    
    
    if(_quiet < 1) printf("Contributions\n");
    sz = getSize(hvcs);
    contribs = getContributions(hvcs); //O(n) because contributions are copied to an array, although if called a second time without any changes in the data structure, then O(1)
    
    if(_quiet < 2){
        
        switch(outputFormat){
            
            case 0:
                printIds(hvcs);
                printf("%-16.15g\n", totalHV(hvcs));
                break;
                
            case 1:
                printIds(hvcs);
                break;
            case 2:
                printf("%-16.15g\n", totalHV(hvcs));
                break;
                
            
            case 3:
                contribs = getContributions(hvcs);
                if(_quiet < 1) printf("\t\tcontributions (size: %d):\n", sz);
                for(j = 0; j < sz; j++){
                    printf("%-16.15g\n", contribs[j]);
                }
                printf("%-16.15g\n", totalHV(hvcs));
                break;
            case 4:
                //index and contribution at time of removal (for the (n-k) points removed)
                //note that, for k=0, it prints the results for all points
                for(j = k; j < n; j++){
                    printf("%d %-16.15g\n", removedId[j], removedContr[j]);
                }
                break;
                
            

        }
    }
    
    if(_quiet < 1) printf("\tfree\n");
    dealloc(hvcs);
    
    return totalhv;
}



/* The fist 'popsize' points in 'data' represent the initial population.
 * This function adds the remaining points in data (n-popsize points) one at a time.
 * After a point is added, the currently least contributor is removed.
 * Flag 'updateContributionsRightAway' indicates whether the contributions should be
 * updated immediately after applying changes to the population (flag is set to 1, and 
 * applies the HVC3D-U update algorithm). If the flag is set to 0, then the contributions
 * will be computed (using HVC3D-R that recomputes all contributions from scratch)
 * only when need, i.e., only the least contributor has to be removed).
 * 
 * mode (0|1|2) - is used only for testing purposes. It also shows different ways of
 *              how the least contributor can be removed.
 */
double boundedArchiveSteadyState(double * data, int d, int n, double * ref, int popsize, int updateContributionsRightAway, int mode){
    
    if(_quiet < 1) printf("--- compute bounded archive steady-state (popsize: %d) ----\n", popsize);
    
    
    int i, ix, removed;
    double * point;
    double contr = 0;
//     int updateContributionsRightAway = 1;
    
    hvc_s * hvcs = init(data, d, popsize, n, ref);
    
    //computes all contributions
    updateAllContributions(hvcs); //should also work without this line

    
    for(i = popsize; i < n; i++){

        if(_quiet < 1) printf("hypervolume of the archive: ");
        if(_quiet < 2) printf("%-16.15g\n", totalHV(hvcs));
        ix = n+1;
        point = &data[d*i];
        
        addPoint(hvcs, point, updateContributionsRightAway); //O(n)
        
        
        double * p = (double *) malloc(d * sizeof(double));
        ix = -1, removed = 1;
        if(_quiet < 1){
            //hv = totalHV(hvcs);
            
            switch(mode){
                case 0:
                    contr = getLeastContribution(hvcs); //mode 1
                    ix = getLeastContributorIndex(hvcs);
                    break;
                case 1:
                    //mode 2
                    ix = getLeastContributorIndex(hvcs);
                    contr = getContributions(hvcs)[ix];
                    break;
                case 2:
                    //mode 3
                    contr = getLeastContributor(hvcs, p);
                    ix = getLeastContributorIndex(hvcs);
                    break;
            }
            printf("Least contribution (ix: %d): %f\n", ix, contr);
//             checkLeastContributor(hvcs, contr);
            
        }

            /* how to get the index and/or the contribution of the least contributor */
            /* //mode 1
                contr = getLeastContribution(hvcs); //mode 1
                ix = getLeastContributorIndex(hvcs);
            */
            /* //mode 2
                ix = getLeastContributorIndex(hvcs);
                contr = getContributions(hvcs)[ix];
            */
            /* //mode 3
                    contr = getLeastContributor(hvcs, p);
                    ix = getLeastContributorIndex(hvcs);
            */

        /* different ways to remove the least contributor */
        switch(mode){
            case 0:
                removeLeastContributor(hvcs, updateContributionsRightAway); //mode 1
                break;
            case 1:
                ////mode 2
                ix = getLeastContributorIndex(hvcs);
                removePointAt(hvcs, ix, updateContributionsRightAway);
                break;
            case 2:
                //mode 3
                contr = getLeastContributor(hvcs, p);
                removed = removePoint(hvcs, p, updateContributionsRightAway);
                
        }
        free(p);
        
    }
    
    double totalhv = totalHV(hvcs);
    

    if(_quiet < 2){
        printf("%-16.15g\n", totalHV(hvcs));
    }
    
    dealloc(hvcs);
    
    return totalhv;
}






//update hv sequentially, starting from the empy set (HV3D+-U)
double sequentialIncrementalHV(double * data, int d, int n, double * ref){
    if(_quiet < 1) printf("--- compute hv sequentially ----\n");
    
    if(_quiet < 1) printf("\tinit empty set (%d dimensions) and set ref = (%f, %f, %f)\n", d, ref[0], ref[1], ref[2]);
    
    hvc_s * hvcs = init(data, d, 0, n, ref);
    int i;
    double * point;
    double totalhv = 0;
    
    
    for(i = 0; i < n; i++){
        point = &data[d*i];
        if(_quiet < 1) printf("\tadd Point: %f %f %f\n", point[0], point[1], point[2]);
        
        //mode 1
        totalhv += addOneContribution(hvcs, point);
        
        /* //mode 2 (slower than the previous one)
        double contr = 0;
        contr = oneContribution(hvcs, point);
        totalhv += contr;
        addPoint(hvcs, point, 0); //O(n)
        */
    }
    
    if(_quiet < 2){
        printf("%-16.15g\n", totalhv);
    }
    dealloc(hvcs);
    
    return totalhv;
}




//update hv sequentially, starting from the empy set (HV3D+-R), but recomputes the
//hypervolume for each new point
double sequentialIncrementalHVSlow(double * data, int d, int n, double * ref){
    if(_quiet < 1) printf("--- compute hv sequentially ----\n");
    
    
    if(_quiet < 1) printf("\tinit empty set (%d dimensions) and set ref = (%f, %f, %f)\n", d, ref[0], ref[1], ref[2]);
    
    hvc_s * hvcs = init(data, d, 0, n, ref);
    int i;
    double * point;
    double totalhv = 0;
    
    for(i = 0; i < n; i++){
        point = &data[d*i];
        
        if(_quiet < 1) printf("\tadd Point: %f %f %f\n", point[0], point[1], point[2]);
        
        addPoint(hvcs, point, 0); //O(n)
        totalhv = updateHypervolume(hvcs);
//         printf("totalhv: %f\n", totalhv);
        
    }
    
    
    if(_quiet < 2){
        printf("%-16.15g\n", totalhv);
    }
    dealloc(hvcs);
    
    return totalhv;
}




//update contributions by recomputing only those that change (HVC3D-U)
//replicate gHSS (incremental greedy), but in O(nk^2)
double replicateGHSS3D(double * data, int d, int n, double * ref, int k, int outputFormat){
    if(_quiet < 1) printf("--- compute gHSS3D (size: %d, k: %d)----\n", n, k);
    
    int maxix;
    double maxc;
    int i, j;
    
    double * contribs = (double *) malloc(n * sizeof(double));
    int * selected = (int *) malloc(n * sizeof(int));
    double * selcontribs = (double *) malloc(n * sizeof(double));
    int * outp = (int *) malloc(n * sizeof(int));
    for(i = 0; i < n ; i++) outp[i] = i;
    
    hvc_s * hvcs = init(data, d, 0, n, ref); //start with an empty set
    
    for(i = 0; i < k; i++){
        maxix = -1;
        maxc = 0;
        for(j = 0; j < n-i; j++){
            contribs[j] = oneContribution(hvcs, &data[outp[j]*d]);
            if(contribs[j] > maxc){
                maxc = contribs[j];
                maxix = j;
            }
        }
        
        addPoint(hvcs, &data[outp[maxix]*d], 0); //No need to update the contribution of the selected points (points in hvcs)
        selected[i] = outp[maxix];
        selcontribs[i] = maxc;
        outp[maxix] = outp[n-i-1]; //the first n-k indexes in outp indicate which points in 'data' were not yet (selected) added to the data structure hvcs
        
    }
    
    free(contribs);
    free(outp);
    double totalhv = totalHV(hvcs);
    
    
//     if(_quiet < 1) printf("Contributions\n");
//     int sz = getSize(hvcs);
//     contribs = getContributions(hvcs); //O(n) because contributions are copied to an array, although if called a second time without any changes in the data structure, then O(1)
    
    if(_quiet < 2){
        double hv = 0;
        switch(outputFormat){
            
            case 0:
                for(j = 0; j < k; j++) printf("%d\n", selected[j]);
                printf("%-16.15g\n", totalHV(hvcs));
                break;
                
            case 1:
                for(j = 0; j < k; j++) printf("%d\n", selected[j]);
                break;
            case 2:
                printf("%-16.15g\n", totalHV(hvcs));
                break;
                
            
            case 3:
//                 contribs = getContributions(hvcs);
                for(j = 0; j < k; j++){
                    printf("%d\t%-16.15g\n", selected[j], selcontribs[j]);
                }
                break;
            case 4:
                //index and contribution at time of removal (for the (n-k) points removed)
                //note that, for k=0, it prints the results for all points
                for(j = 0; j < k; j++){
                    hv += selcontribs[j];
                    printf("%d\t%-16.15g\n", selected[j], hv);
                }
                break;
        }
    }
    
    free(selected);
    free(selcontribs);
    
    if(_quiet < 1) printf("\tfree\n");
    dealloc(hvcs);
    
    return totalhv;
}




double runExample(double * data, int d, int n, double *ref, int exampleID, int mu, int lambda, int quiet, int checkCorrectness, int mode, int outputFormat){
//     _quiet = 1;
    _quiet = quiet;
    _checkCorrectness = checkCorrectness;
    _time_elapsed_cpu = 0;
    double res = 1;
    int minpop = max(1, n/10);
    int popsize = (mu <= 0 || mu >= n) ? minpop : mu; //to assure a popsize smaller than n (just for testing)
    int k = (mu >= 0) ? min(mu, n) : n/2; //default
    lambda = (lambda > 0 && lambda < (n-k)) ? lambda : max((n-k)/10, 1);
//     printf("quiet: %d\n", _quiet);
//     printf("(mu+lambda) = (%d + %d)\n", popsize, lambda);
    
    Timer_start ();
    
    switch(exampleID){
        case 0: testInitClose(data, d, n, ref); break;
        case 1: res = computeHypervolumeOnce(data, d, n, ref); break;
        case 2: res = computeContributionsOnce(data, d, n, ref); break;
        
        case 3: res = replicateGHSSD3D(data, d, n, ref, k, 1, outputFormat); break; //fast version
        case 103: res = replicateGHSSD3D(data, d, n, ref, k, 0, outputFormat); break; //slower version
       
        case 4: res = unboundedArchive(data, d, n, ref, outputFormat); break;
        case 104: res = unboundedArchiveSlower(data, d, n, ref, outputFormat); break;
        
        case 5: res = boundedArchiveSteadyState(data, d, n, ref, popsize, 1, mode); break; //fast version
        case 105: res = boundedArchiveSteadyState(data, d, n, ref, popsize, 0, mode); break; //slower version
        case 8: res = sequentialIncrementalHV(data, d, n, ref); break;
        case 108: res = sequentialIncrementalHVSlow(data, d, n, ref); break;
        
        case 9: res = replicateGHSS3D(data, d, n, ref, k, outputFormat); break;
        
    }
    if(res < 0){
        printf("%f\n", res);
    }
    
    _time_elapsed_cpu = Timer_elapsed_virtual ();
    
    return _time_elapsed_cpu;
        
    
}
