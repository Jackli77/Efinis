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

#ifndef NOINERTIA
double inertiaIntegrate(femMesh *theMesh, femIntegration *theRule, double rho)
{
    double I = 0;
 
//
// A completer :-)
//
 
    return I;
}
#endif

#ifndef NOREAD
femMesh *inertiaMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));


    int i,trash;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }
    
    ErrorScan(fscanf(file, "Number of triangles %d \n", &theMesh->nElem));
    theMesh->elem = malloc(sizeof(int) * theMesh->nElem);
    for (i = 0; i < theMesh->nNode; i++) {
        ErrorScan(fscanf(file, "%d : %d %d %d \n", &trash, &theMesh->elem[3*i], &theMesh->elem[3*i+1], &theMesh->elem[3 * i + 2]));
    }

//
// A completer :-)
//     Allocation dynamique du tableau d'appartenance
//     Lecture du tableau d'appartenance
// 
// Le premiere partie fournie de la fonction peut largement servir d'inspiration....
// N'oubliez pas : Google est votre ami : par exemple google C malloc peut etre tres utile !
//     
  


    fclose(file);
    return theMesh;
}
#endif

#ifndef NOFREE
void inertiaMeshFree(femMesh *theMesh)
{

//
// A completer :-)
//
//     Libérer la mémoire dynamique allouée 
// 
    if (theMesh) {
        free(theMesh->nElem);
        free(theMesh->nNode);
        free(theMesh->X);
        free(theMesh->Y);
        free(theMesh->elem);
    }
}
#endif

