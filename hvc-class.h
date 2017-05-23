#ifndef HVC_CLASS_H_
#define HVC_CLASS_H_

//These functions are available only for 3D


/* -------- data strucutre -------- */
typedef struct hvcstruct hvc_s;


/* -------- init function --------*/
/**
 * Pre: 'data' is expected to represent a nondominated point set of 'n' points in 'd' dimensions. 
 * TODO: Any dominated points detected at this phase are immediately discarded!
 */
hvc_s * init(double * data, int d, int n, int naloc, double *ref);

/* -------- compute and update the data strucutre --------*/
//operation functions
//contributions/hv is (re)computed only if points were added/removed (without contributions/hv update) since the last time these were computed
void updateAllContributions(hvc_s * hvcs);

double addPoint(hvc_s * hvcs, double * point, int updateContribs);

int removePoint(hvc_s * hvcs, double * point, int updateContribs);
int removePointAt(hvc_s * hvcs, int i, int updateContribs);
void removeLeastContributor(hvc_s * hvcs, int updateContribs);


//one contribution and hypervolume computation
double oneContribution(hvc_s * hvcs, double * point); //3D only - O(n)
double addOneContribution(hvc_s * hvcs, double * point);
double updateHypervolume(hvc_s * hvcs);



void changeReferencePoint(hvc_s * hvcs, double * ref);

/* -------- get Information -------- */
double * getContributions(hvc_s * hvcs);
double totalHV(hvc_s * hvcs);
double * getPoints(hvc_s * hvc_s); //the order is the same as contributions - contributions are updated after this call (even if they were not before the call)
int getSize(hvc_s * hvcs);
int getAllocSize(hvc_s * hvcs);

int isUpToDate(hvc_s * hvcs);


int getLeastContributorIndex(hvc_s * hvcs);
double getLeastContributor(hvc_s * hvcs, double * point);
double getLeastContribution(hvc_s * hvcs);


/* -------- setters -------- */
void forceRecomputationNext(hvc_s * hvcs); //O(1)


/* -------- close function -------- */
double dealloc(hvc_s * hvcs);


/* -------- others --------*/
// TODO
// hvc_s * copy(hvc_s * hvcs);
void printIds(hvc_s * hvcs); //private - just for testing gHSSD
int getLeastContributorId(hvc_s * hvcs); //private - just for testing gHSSD

#endif
