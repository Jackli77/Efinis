#include "glfem.h"

int main(void)
{  
    femMesh  *theMesh  = femMeshRead("C:/Users/jackl/Documents/Efinis/Edges/data/gear8.txt");    
    femEdges *theEdges = femEdgesCreate(theMesh);    
    
    edgesExpand(theEdges);               //   femEdgesPrint(theEdges);
    edgesSort(theEdges);                 //   femEdgesPrint(theEdges);
    edgesShrink(theEdges);               //   femEdgesPrint(theEdges);
    printf("Boundary edges  : %i \n", theEdges->nBoundary);
    printf("Boundary length : %14.7e \n", edgesBoundaryLength(theEdges));

    char theMessage[256];
    sprintf(theMessage, "Boundary edges : %i", theEdges->nBoundary);
      
//
//  On superpose le maillage (en bleu), 
//  tous les segments frontieres (en noir),
//  et la frontiere (en rouge)
//
//  Au depart de votre travail, vous devriez obtenir un maillage bleu....
//  et a la fin de l'exercice un maillage noir avec bord rouge :-)
//

    glfemWindowCreate("EPL1110 : Edges",480,480,theMesh->nNode,theMesh->X,theMesh->Y);
    do
    {
        glfemReshape(theMesh->X,theMesh->Y,theMesh->nNode);  
        
        glfemSetLineWidth(0.0001); 
        glfemSetColor(GLFEM_BLUE); 	  glfemPlotMesh(theMesh);  
        glfemSetColor(GLFEM_BLACK); 	glfemPlotEdges(theEdges);        
        glfemSetLineWidth(0.0015);
        glfemSetColor(GLFEM_RED); 	  glfemPlotBnd(theEdges);
       
        glfemDrawMessage(theMessage,(double[2]){16.0, 30.0}); 
        glfemWindowUpdate();  
    } while(!glfemWindowShouldClose());

    glfemWindowFree();
    femEdgesFree(theEdges);
    femMeshFree(theMesh);
    exit(EXIT_SUCCESS);
    return 0;
    
 
}





