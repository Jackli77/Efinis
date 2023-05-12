#include "fem.h"
#include <math.h>

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)


double* femElasticitySolve(femProblem* theProblem)
{
    femFullSystem* theSystem = theProblem->system;
    femIntegration* theRule = theProblem->rule;
    femDiscrete* theSpace = theProblem->space;
    femGeo* theGeometry = theProblem->geometry;
    femNodes* theNodes = theGeometry->theNodes;
    femMesh* theMesh = theGeometry->theElements;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];

    int nLocal = theMesh->nLocalNode;

    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    double** A = theSystem->A;
    double* B = theSystem->B;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j = 0; j < nLocal; j++) {
            map[j] = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            //The factor variable is introduced to account for the axisymmetric behavior.
            //It includes the factor 2.0 * M_PI * x[i] * jac * weight, where x[i] represents the x-coordinate of the node i.
            //In the loop where the stiffness matrix A is assembled, the contributions from the x and y components are multiplied by the factor to incorporate the axisymmetric behavior.

                In the loop where the load vector B is assembled, the contribution is multiplied by the factor to account for the axisymmetric effect.
            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {
                    // Axisymmetric modifications
                    double factor = 2.0 * M_PI * x[i] * jac * weight;
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * factor;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * factor;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * factor;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * factor;
                }
            }
            for (i = 0; i < theSpace->n; i++) {
                // Axisymmetric modification
                B[mapY[i]] -= phi[i] * g * rho * jac * weight * 2.0 * M_PI * x[i];
            }
        }
    }

    int* theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }
    // Apply Neumann boundary conditions
    int* neumannEdges = theProblem->neumannEdges;
    double* neumannValues = theProblem->neumannValues;
    for (iEdge = 0; iEdge < theMesh->nEdge; iEdge++) {
        if (neumannEdges[iEdge]) {
            for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                double xsi = theRule->xsi[iInteg];
                double eta = theRule->eta[iInteg];
                double weight = theRule->weight[iInteg];
                femDiscretePhi1(theSpace, xsi, phi);
                double length = theMesh->edges[iEdge]->length;
                double neumannFactor = length * weight;
                for (i = 0; i < theSpace->n; i++) {
                    B[mapY[i]] += phi[i] * neumannValues[iEdge] * neumannFactor;
                }
            }
        }
    }


    return femFullSystemEliminate(theSystem);
}
