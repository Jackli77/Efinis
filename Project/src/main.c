/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"
#include <math.h>

int main(void)
{  
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("C:/Users/jackl/Documents/Efinis/Project/data/mesh.txt");
    femProblem* theProblem = femElasticityRead(theGeometry,"C:/Users/jackl/Documents/Efinis/Project/data/problem.txt");
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"C:/Users/jackl/Documents/Efinis/Project/data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"C:/Users/jackl/Documents/Efinis/Project/data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();
    return 0;  
}

 
