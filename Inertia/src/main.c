/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Integration of Baltic sea :-)
 *
 *  Copyright (C) 2022 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include <stdio.h>
#include <math.h>


#include "fem.h"
#include "glfem.h"

femMesh    *inertiaMeshRead(const char *filename);
double      inertiaSteelRho();
double      inertiaIntegrate(femMesh *theMesh, femIntegration *theRule, double rho);
void        inertiaMeshFree(femMesh *theMesh);


int main(void)
{   
    femMesh *theMesh = inertiaMeshRead("../data/gear60.txt");
    femIntegration *theRule = femIntegrationCreate(3,FEM_TRI);
    double rho = inertiaSteelRho();
    double   I = inertiaIntegrate(theMesh,theRule,rho); 
    
    double *theField = malloc(sizeof(double)*theMesh->nNode); int i;
    for (i=0; i < theMesh->nNode; ++i) {
        double xLoc = theMesh->X[i];
        double yLoc = theMesh->Y[i];
        theField[i]  = xLoc*xLoc + yLoc*yLoc; }
        
    char theMessage[256];       
    sprintf(theMessage,"Inertia = %14.7e kg m^2",I);
    printf("%s\n",theMessage);


//
//  Partie graphique :-)
//  Peut être commentée si vous n'arrivez pas à compiler la partie graphique
//  Il faut aussi retirer l'include de glfem.h et retirer la partie graphique de
//  de la compilation.  Il suffit alors de compiler main.c fem.c homework.c
//
//  So easy !
// 
 
    glfemWindowCreate("EPL1110 : Inertia",400,400);

    do
    {
        int w,h;
        glfemReshape(theMesh);
        glfemPlotField(theMesh,theField); 
        glfemDrawMessage(theMessage,(double[2]){20,380});
        glfemWindowUpdate();
         
    } while(!glfemWindowShouldClose());

    glfemWindowFree();
    femIntegrationFree(theRule);
    inertiaMeshFree(theMesh);
    exit(EXIT_SUCCESS);
    return 0;  

//
//  Fin de la partie graphique du code !
//    
    
  
}

 
