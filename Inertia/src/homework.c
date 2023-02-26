#include "fem.h"


#ifndef NORHOSTEEL
double inertiaSteelRho()
{

    //
    // Modifier la valeur pour avoir la masse volumique de l'acier [kg/m3]
    // Une tolerance de 10% sur la valeur est admise
    //

    double rho = 7800.0;
    return rho;
}
#endif
double interpolate(double u[3], double xsi, double eta) {
    return u[0] * xsi + u[1] * eta + u[2] * (1.0 - xsi - eta);
}

#ifndef NOINERTIA
double inertiaIntegrate(femMesh* theMesh, femIntegration* theRule, double rho)
{
    double I = 0, Ipart = 0;
    double jac = 0;
    double xLoc[3];
    double yLoc[3];
    double xInteg = 0, yInteg = 0;

    for (int i = 0; i < theMesh->nElem; i++) {
        Ipart = 0;
        for (int j = 0; j < 3; j++) {
            xLoc[j] = theMesh->X[theMesh->elem[3 * i + j]];
            yLoc[j] = theMesh->Y[theMesh->elem[3 * i + j]];
        }
        jac = fabs((xLoc[0] - xLoc[1]) * (yLoc[0] - yLoc[2]) - (xLoc[0] - xLoc[2]) * (yLoc[0] - yLoc[1]));

        for (int j = 0; j < theRule->n; j++) {
            xInteg = interpolate(xLoc, theRule->xsi[j], theRule->eta[j]);
            yInteg = interpolate(yLoc, theRule->xsi[j], theRule->eta[j]);
            Ipart += (xInteg * xInteg + yInteg * yInteg) *theRule->weight[j];
        }
        I += Ipart * jac;
        return I * rho * 1e-7;
    }
}
#endif

#ifndef NOREAD
femMesh* inertiaMeshRead(const char* filename) {
    femMesh *theMesh = malloc(sizeof(femMesh));
    

        int i, trash;
        int trois = 3;
        theMesh->nLocalNode = trois;

        FILE* file = fopen(filename, "r");
        if (file == NULL) Error("No mesh file !");

        ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
        theMesh->X = malloc(sizeof(double)*theMesh->nNode);
        theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
        if (theMesh->X == NULL)
            return NULL;
        if (theMesh->Y == NULL) 
            return NULL;
            for (i = 0; i < theMesh->nNode; ++i) {
                ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theMesh->X[i], &theMesh->Y[i]));
            }
        ErrorScan(fscanf(file, "Number of triangles %d \n", &theMesh->nElem));
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        if (theMesh->elem == NULL)
            return NULL;
        for (i = 0; i < theMesh->nElem; ++i) {
            ErrorScan(fscanf(file, "%d : %d %d %d \n", &trash, &theMesh -> elem[3 * i], &theMesh ->elem[3 * i + 1], &theMesh -> elem[3 * i + 2]));

        }

        fclose(file);
    return theMesh;
}
#endif

#ifndef NOFREE
void inertiaMeshFree(femMesh* theMesh)
{
    if (theMesh) {
        free(theMesh->X);
        free(theMesh->Y);
        free(theMesh->elem);
        free(theMesh);
    }
}
#endif