/*
 *  main.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat    
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>

int main(void)
{  

    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    F-B-I : Full solver - Band solver - Iterative solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    printf("    < >   : Unrefine - refine the mesh \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_NO;
    int testConvergence;
    
 
    femMesh *theMesh = femMeshCreateBasicSquare(2);
  //  femMeshWrite(theMesh,"../data/msh2.txt");
    femDiffusionProblem* theProblem = femDiffusionCreate(theMesh,solverType,renumType);
    femDiffusionSetSource(theProblem,convergenceSource);
     
    clock_t tic = clock();
    do {
        femDiffusionCompute(theProblem);  
        femSolverPrintInfos(theProblem->solver); 
        testConvergence = femSolverConverged(theProblem->solver); }
    while ( testConvergence == 0);
    if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    
    int iRef=0, nRef= 2, maxRef = 8;
    double errors[maxRef], rate = 0;
    femDiffusionComputeError(theProblem,convergenceSoluce);
    printf("    Error Solution L2      : %10.4e \n",theProblem->errorSoluceL2);
    printf("    Error Interpolation L2 : %10.4e \n",theProblem->errorInterpolationL2);
    printf("    Error Solution H1      : %10.4e \n",theProblem->errorSoluceH1);
    printf("    Error Interpolation H1 : %10.4e \n",theProblem->errorInterpolationH1);
    errors[iRef] = theProblem->errorSoluceL2; 
    fflush(stdout);  
    
    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("EPL1110 : Manufactured solutions ");
    glfwMakeContextCurrent(window);
     

    do 
    {
       
        
        int testConvergence,w,h;
        char theMessage[256];
        sprintf(theMessage, "Estimated rate : %.2f ",rate);
        
        glfwGetFramebufferSize(window,&w,&h);
          
        if (option == 1) {
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotField(theProblem->mesh,theProblem->soluce);   }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
        
        
        int action = glfemGetAction();        
        if (action != 0) {
            if (action == -1  && nRef !=  2)        {nRef = nRef/2; iRef = iRef - 1;}
            if (action == 1   && iRef != maxRef-1 ) {nRef = nRef*2; iRef = iRef + 1;}
            femMeshFree(theMesh);
            theMesh =  femMeshCreateBasicSquare(nRef);
            solverType = UNDEFINED;  }
        
        
        
        if (solverType != newSolverType || renumType != newRenumType ) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(theMesh,solverType,renumType);
            femDiffusionSetSource(theProblem,convergenceSource);
   
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            
            femDiffusionComputeError(theProblem,convergenceSoluce);
            printf("    Error Solution L2      : %10.4e \n",theProblem->errorSoluceL2);
            printf("    Error Interpolation L2 : %10.4e \n",theProblem->errorInterpolationL2);
            printf("    Error Solution H1      : %10.4e \n",theProblem->errorSoluceH1);
            printf("    Error Interpolation H1 : %10.4e \n",theProblem->errorInterpolationH1);
            errors[iRef] = theProblem->errorSoluceL2;   
            if (iRef != 0){
              rate = convergenceEstimateRate(&errors[0],iRef+1,2.0);
              printf("    Mean estimated L2 rate of convergence  : %10.4e \n",rate);
              rate = convergenceEstimateRate(&errors[iRef-1],2,2.0);
              printf("    Last estimated L2 rate of convergence  : %10.4e \n",rate); }
            
            fflush(stdout); }           
          
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
       
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed


 
   
   
               
    glfwTerminate(); 
  //  femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

