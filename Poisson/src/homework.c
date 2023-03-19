#include"fem.h"


# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename) {
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh = femMeshRead(filename);
    theProblem->edges = femEdgesCreate(theProblem->mesh);
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4, FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
    } else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
    }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem) {
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}


# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y) {
    int nbrnoeuds = theMesh->nLocalNode;

    for (int j = 0; j < nbrnoeuds; ++j) { // Lister les nœuds locaux
        int noeudcorrespondant = theMesh->elem[nbrnoeuds * i + j]; // *i pour atteindre l'élément recherché
        map[j] = noeudcorrespondant;
        x[j] = theMesh->X[map[j]];
        y[j] = theMesh->Y[map[j]];
        // Le tableau map : pour un seul élément à la fois (la fonction sera appelée pour chaque élément)
    }
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem) {
    if (theProblem->space->n > 4) {
        Error("Yeeek piyou piyou, c tro gran");
    }
    for (int j = 0; j < theProblem->mesh->nElem; ++j) {
        double x[4];
        double y[4];
        double phi[4];
        double dphidksi[4];
        double dphideta[4];
        double dphidx[4];
        double dphidy[4];
        int map[4];

        // Va être appelée à chaque itération
        femMeshLocal(theProblem->mesh, j, map, x, y);

        for (int i = 0; i < theProblem->rule->n; ++i) {
            double poids = theProblem->rule->weight[i];
            double ksi = theProblem->rule->xsi[i];
            double eta = theProblem->rule->eta[i];

            femDiscretePhi2(theProblem->space, ksi, eta, phi);
            femDiscreteDphi2(theProblem->space, ksi, eta, dphidksi, dphideta);

            double dxdksi = 0;
            double dydksi = 0;
            double dxdeta = 0;
            double dydeta = 0;

            for (int k = 0; k < theProblem->space->n; ++k) {
                dxdksi += x[k] * dphidksi[k];
                dxdeta += x[k] * dphideta[k];
                dydksi += y[k] * dphidksi[k];
                dydeta += y[k] * dphideta[k];
            }

            double jacob = dxdksi * dydeta - dxdeta * dydksi;
            if (jacob < 0) {
                int noeud = theProblem->mesh->elem[theProblem->mesh->nLocalNode * j];
                theProblem->mesh->elem[theProblem->mesh->nLocalNode * j] = theProblem->mesh->elem[
                        theProblem->mesh->nLocalNode * j + 2];
                theProblem->mesh->elem[theProblem->mesh->nLocalNode * j + 2] = noeud;
            }

            jacob = fabs(jacob);
            for (int l = 0; l < theProblem->space->n; ++l) {
                dphidx[l] = (dphidksi[l] * dydeta - dphideta[l] * dydksi) / jacob;
                dphidy[l] = (dphideta[l] * dxdksi - dphidksi[l] * dxdeta) / jacob;
            }

            for (int k = 0; k < theProblem->space->n; ++k) {
                theProblem->system->B[map[k]] += poids * jacob * phi[k];
                for (int l = 0; l < theProblem->space->n; ++l) {
                    theProblem->system->A[map[k]][map[l]] +=
                            jacob * poids * (dphidx[k] * dphidx[l] + dphidy[k] * dphidy[l]);
                }
            }
        }
    }

    femEdges *edges = theProblem->edges;
    double valeur = 0.0;
    for (int i = 0; i < edges->nEdge; ++i) {
        //femEdge *edge = &(edges->edges[i]);
        if (edges->edges[i].elem[1] < 0) {
            for (int j = 0; j < 2; ++j) {
                femFullSystemConstrain(theProblem->system, edges->edges[i].node[j], valeur);
            }
        }
    }
    femFullSystemEliminate(theProblem->system);
}
# endif//;

void femPoissonSolve2(femPoissonProblem *theProblem) {
    // Déclaration des variables
    femMesh *theMesh = theProblem->mesh;
    femEdges *theEdges = theProblem->edges;
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;

    if (theSpace->n > 4) {
        Error("Yeeek piyou piyou, c tro gran");
    }

    for (int j = 0; j < theMesh->nElem; ++j) { // Parcourir les éléments un par un
        double x[4];
        double y[4];
        double phi[4];
        double dphidksi[4];
        double dphideta[4];
        double dphidx[4];
        double dphidy[4];
        int map[4];

        femMeshLocal(theMesh, j, map, x, y);

        for (int j = 0; j < theRule->n; ++j) { // theProblem->rule : règle d'intégration
            double poids = theRule->weight[j];
            double ksi = theRule->xsi[j];
            double eta = theRule->eta[j];

            femDiscretePhi2(theSpace, ksi, eta, phi);
            femDiscreteDphi2(theSpace, ksi, eta, dphidksi, dphideta);

            double dxdksi = 0;
            double dydksi = 0;
            double dxdeta = 0;
            double dydeta = 0;

            for (int k = 0; k < theSpace->n; ++k) { // Voir notes manuscrites (CM4 dernière page en haut de la feuille)
                dxdksi += x[k] * dphidksi[k];
                dxdeta += x[k] * dphideta[k];
                dydksi += y[k] * dphidksi[k];
                dydeta += y[k] * dphideta[k];
            }

            double jacob = fabs(dxdksi * dydeta - dxdeta * dydksi); // Je ne sais pas si la condition de if jacob < 0 est vraiment nécéssaire
            // Cette condition est utile seulement si les nœuds sont mal orientés
            for (int l = 0; l < theSpace->n; ++l) { // Toujours voir la dernière feuille de CM4
                dphidx[l] = (dphidksi[l] * dydeta - dphideta[l] * dydksi) / jacob;
                dphidy[l] = (dphideta[l] * dxdksi - dphidksi[l] * dxdeta) / jacob;
            }

            // Assemblage
            for (int k = 0; k < theSpace->n; ++k) {
                theSystem->B[map[k]] += poids * jacob * phi[k];
                for (int l = 0; l < theSpace->n; ++l) {
                    theSystem->A[map[k]][map[l]] +=
                            jacob * poids * (dphidx[k] * dphidx[l] + dphidy[k] * dphidy[l]);
                }
            }
        }


    }
}