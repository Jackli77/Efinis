#include"fem.h"


#ifndef NORENUMBER 
femMesh* theMeshglob;
int comparateurx(const void* n1, const void* n2) {
    int* un = (int*)n1;
    int* deux = (int*)n2;
    double test = theMeshglob->X[*un] - theMeshglob->X[*deux];
    if (test < 0) {
        return 1;
    }
    if (test > 0) {
        return -1;
    }
    return 0;

}
int comparateury(const void* n1, const void* n2) {
    int* un = (int*)n1;
    int* deux = (int*)n2;
    double test = theMeshglob->Y[*un] - theMeshglob->Y[*deux];
    if (test < 0) {
        return 1;
    }
    if (test > 0) {
        return -1;
    }
    return 0;

}
void femMeshRenumber(femMesh* theMesh, femRenumType renumType)
{
    int i;
    theMeshglob = theMesh;
    int* sch = malloc(theMesh->nNode * sizeof(int));

    switch (renumType) {
    case FEM_NO:
        for (i = 0; i < theMesh->nNode; i++)
            theMesh->number[i] = i;
        break;
        // 
        // A modifier :-)
        // debut
        //
    case FEM_XNUM:

        for (i = 0; i < theMesh->nNode; i++) {
            sch[i] = i;
        }
        qsort(theMesh->number, theMesh->nNode, sizeof(int), comparateurx);
        for (i = 0; i < theMesh->nNode; i++)
        {
            theMesh->number[sch[i]] = i;
        }
        free(sch);
        break;

    case FEM_YNUM:

        for (i = 0; i < theMesh->nNode; i++) {
            sch[i] = i;
        }
        qsort(theMesh->number, theMesh->nNode, sizeof(int), comparateury);
        for (i = 0; i < theMesh->nNode; i++)
        {
            theMesh->number[sch[i]] = i;
        }
        free(sch);
        break;
        // 
        // end
        //

    default: Error("Unexpected renumbering option");
    }
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh* theMesh)
{
    int myBand = 0;
    for (int i = 0; i < theMesh->nElem; i++) {
        for (int j = 0; j < theMesh->nLocalNode; j++) {
            for (int k = 0; k < theMesh->nLocalNode; k++) {
                int diff = fabs(theMesh->number[theMesh->elem[theMesh->nLocalNode * i + j]] - theMesh->number[theMesh->elem[theMesh->nLocalNode * i + k]]);
                if (diff > myBand) {
                    myBand = diff;
                }
            }
        }
    }
    return(myBand + 1);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double* Aloc, double* Bloc, int* map, int nLoc)
{
    for (int i = 0; i < nLoc; i++) {
        int row = map[i];
        for (int j = 0; j < nLoc; j++) {
            int col = map[j];
            if (col >= row) {
                myBandSystem->A[row][col] += Aloc[i * nLoc + j];
            }
            myBandSystem->B[row] += Bloc[i];
        }


    }
}


#endif
#ifndef NOBANDELIMINATE


double* femBandSystemEliminate(femBandSystem* myBand)
{
    double** A, * B, factor;
    int     i, j, k, jend, size, band;
    A = myBand->A;
    B = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Gauss elimination */

    for (k = 0; k < size; k++) {
        if (fabs(A[k][k]) <= 1e-8) {
            printf("Pivot index %d  ", k);
            printf("Pivot value %e  ", A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }
        for (i = k + 1; i < fmin(k + band, size); i++) {
            factor = A[k][i] / A[k][k];
            for (j = k + 1; j < fmin(k + band, size); j++)
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor;
        }
    }

    /* Back-substitution */

    for (i = size - 1; i >= 0; i--) {
        factor = 0;
        for (j = i + 1; j < fmin(i + band, size); j++)
            factor += A[i][j] * B[j];
        B[i] = (B[i] - factor) / A[i][i];
    }

    return(myBand->B);
}


#endif