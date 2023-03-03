#include "fem.h"

# ifndef NOEXPAND

void edgesExpand(femEdges* theEdges)
{
    int n = theEdges->mesh->nLocalNode;
    int i;

    for (i = 0; i < theEdges->mesh->nElem; i++) {
        for (int j = n * i; j < n * i + n; j++) {
            theEdges->edges[j].elem[0] = i;
            theEdges->edges[j].elem[1] = -1;
            theEdges->edges[j].node[0] = theEdges->mesh->elem[j];
            if (j + 1 == n * i + n) {
                theEdges->edges[j].node[1] = theEdges->mesh->elem[n * i];
            }
            else {
                theEdges->edges[j].node[1] = theEdges->mesh->elem[j + 1];
            }

        }

    }
    theEdges->nEdge = n * i;
}

# endif
# ifndef NOSORT

void edgesSort(femEdges* theEdges)
{
    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOCOMPARE

int edgesCompare(const void* e0, const void* e1)
{
    int e00 = ((femEdge*)e0)->node[0];
    int e01 = ((femEdge*)e0)->node[1];
    int e10 = ((femEdge*)e1)->node[0];
    int e11 = ((femEdge*)e1)->node[1];
    int min0 = e00;
    int min1 = e10;
    if (e00 > e01) {
        min0 = e01;
    }
    if (e10 > e11) {
        min1 = e11;
    }
    if (min0 < min1) {
        return 1;
    }
    else if (min0 > min1) {
        return -1;
    }
    else {
        int max0 = e00;
        int max1 = e10;
        if (e10 < e11) {
            max1 = e11;
        }
        if (e00 < e01) {
            max0 = e01;
        }

        if (max0 > max1) {
            return -1;
        }
        if (max0 < max1) {
            return 1;
        }

        return 0;

    }
}

# endif
# ifndef NOSHRINK

void edgesShrink(femEdges* theEdges)
{
    // Ici commence votre contribution a faire 🙂

    int n = 0;
    // Nouveau nombre total de segments : A MODIFIER
    int nBoundary = 0;  // Nombre de segments frontieres : A MODIFIER
    for (int i = 0; i < theEdges->nEdge - 1; i++) {
        int n00 = theEdges->edges[i].node[0];
        int n01 = theEdges->edges[i].node[1];
        int n10 = theEdges->edges[i + 1].node[0];
        int n11 = theEdges->edges[i + 1].node[1];
        if ((n00 == n10 & n01 == n11) || (n00 == n11 & n01 == n10)) {
            theEdges->edges[n].node[0] = n00;
            theEdges->edges[n].node[1] = n01;
            theEdges->edges[n].elem[0] = theEdges->edges[i].elem[0];
            theEdges->edges[n].elem[1] = theEdges->edges[i + 1].elem[0];
            n++;
            i++;
        }
        else {
            theEdges->edges[n].node[0] = n00;
            theEdges->edges[n].node[1] = n01;
            theEdges->edges[n].elem[0] = theEdges->edges[i].elem[0];
            theEdges->edges[n].elem[1] = theEdges->edges[i].elem[1];
            nBoundary++;
            n++;
            if (i + 2 == theEdges->nEdge) {
                theEdges->edges[n].node[0] = n10;
                theEdges->edges[n].node[1] = n11;
                theEdges->edges[n].elem[0] = theEdges->edges[i + 1].elem[0];
                theEdges->edges[n].elem[1] = theEdges->edges[i + 1].elem[1];
                nBoundary++;
                n++;


            }



        }


    }


    // Ici, finit votre contribution

    // Reallocation du tableau des edges

    theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;
}

# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges* theEdges)
{
    double L = 0.0;
    for (int i = 0; i < theEdges->nEdge; i++) {

        if (theEdges->edges[i].elem[1] == -1) {
            double x2 = theEdges->mesh->X[theEdges->edges[i].node[1]];
            double x1 = theEdges->mesh->X[theEdges->edges[i].node[0]];
            double y2 = theEdges->mesh->Y[theEdges->edges[i].node[1]];
            double y1 = theEdges->mesh->Y[theEdges->edges[i].node[0]];
            double Lx = x2 - x1;
            double Ly = y2 - y1;
            L += sqrt(Lx * Lx + Ly * Ly);

        }
    }
    return L;
}

# endif