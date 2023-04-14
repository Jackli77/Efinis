#include"fem.h"

double convergenceSource(double x, double y) {

    double lap;
    lap = (20 * sqrt(2) * (1 - x) * x * (1 - y)) / (pow(16 - 10 * sqrt(2) * (x + y), 2) + 1) +
        (20 * sqrt(2) * (1 - x) * y * (1 - y)) / (pow(16 - 10 * sqrt(2) * (x + y), 2) + 1) -
        (20 * sqrt(2) * x * y * (1 - y)) / (pow(16 - 10 * sqrt(2) * (x + y), 2) + 1) +
        (800 * (1 - x) * x * y * (1 - y) * (16 - 10 * sqrt(2) * (x + y))) /
        pow(pow(16 - 10 * sqrt(2) * (x + y), 2) + 1, 2) -
        (20 * sqrt(2) * (1 - x) * x * y) / (pow(16 - 10 * sqrt(2) * (x + y), 2) + 1) -
        2 * (x - 1) * x * atan(16 - 10 * sqrt(2) * (x + y)) - 2 * (y - 1) * y * atan(16 - 10 * sqrt(2) * (x + y));

    return -lap;
}

double convergenceSoluce(double x, double y, double* u) {
    double a = sqrt(2.0);
    u[0] = x * y * (1 - x) * (1 - y) * atan(20 * (x + y) / a - 16);
    u[1] = y * (y - 1) * (10 * a * x * (x - 1)
        + (2 * x - 1) * (4 * pow(5 * a * (x + y) - 8, 2) + 1) * atan(10 * a * (x + y) - 16))
        / (4 * pow(5 * a * (x + y) - 8, 2) + 1);
    u[2] = x * (x - 1) * (10 * a * y * (y - 1)
        + (2 * y - 1) * (4 * pow(5 * a * (x + y) - 8, 2) + 1) * atan(10 * a * (x + y) - 16))
        / (4 * pow(5 * a * (x + y) - 8, 2) + 1);
    return u[0];
}


double convergenceEstimateRate(double* errors, int n, double ratio) {


    double rate = 0.0;
    int i;
    for (i = 0; i < n - 1; i++) {
        rate += log10(errors[i] / errors[i + 1]) / log10(ratio);
    }
    rate = rate / (n - 1);
    return rate;
}


void femDiffusionComputeError(femDiffusionProblem* theProblem,
    double(*soluce)(double, double, double*)) {

    femMesh* theMesh = theProblem->mesh;
    femIntegration* theRule = theProblem->rule;
    femDiscrete* theSpace = theProblem->space;
    if (theSpace->n >= 5) {
        Error("Space size is too big !");
    }
    double Xl[4], Yl[4], phi[4], dphidxsi[4], dphideta[4], Ul[4], Uapprox[4];
    double Usol[3];
    int iElem, iInteg, i, map[4], ctr[4];


    theProblem->errorSoluceL2 = 0.0;
    theProblem->errorSoluceH1 = 0.0;
    theProblem->errorInterpolationL2 = 0.0;
    theProblem->errorInterpolationH1 = 0.0;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femDiffusionMeshLocal(theProblem, iElem, map, ctr, Xl, Yl, Ul);
        for (i = 0; i < theSpace->n; i++)
            Uapprox[i] = soluce(Xl[i], Yl[i], Usol);
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            double x = 0;
            double y = 0;
            double u = 0;
            double uapproxf = 0;
            double dudx = 0;
            double dudy = 0;
            double duapproxfdx = 0;
            double duapproxfdy = 0;
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {
                x += Xl[i] * phi[i];
                y += Yl[i] * phi[i];
                u += Ul[i] * phi[i];
                uapproxf += Uapprox[i] * phi[i];
                dxdxsi += Xl[i] * dphidxsi[i];
                dxdeta += Xl[i] * dphideta[i];
                dydxsi += Yl[i] * dphidxsi[i];
                dydeta += Yl[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {
                double dphidx = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                double dphidy = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
                dudx += Ul[i] * dphidx;
                dudy += Ul[i] * dphidy;
                duapproxfdx += Uapprox[i] * dphidx;
                duapproxfdy += Uapprox[i] * dphidy;
            }

            soluce(x, y, Usol);
            double e = (u - Usol[0]);
            double dedx = (dudx - Usol[1]);
            double dedy = (dudy - Usol[2]);
            theProblem->errorSoluceL2 += jac * e * e * weight;
            theProblem->errorSoluceH1 += jac * (e * e + dedx * dedx + dedy * dedy) * weight;

            e = (uapproxf - Usol[0]);
            dedx = (duapproxfdx - Usol[1]);
            dedy = (duapproxfdy - Usol[2]);
            theProblem->errorInterpolationL2 += jac * e * e * weight;
            theProblem->errorInterpolationH1 += jac * (e * e + dedx * dedx + dedy * dedy) * weight;
        }
    }

    theProblem->errorSoluceL2 = sqrt(theProblem->errorSoluceL2);
    theProblem->errorSoluceH1 = sqrt(theProblem->errorSoluceH1);
    theProblem->errorInterpolationL2 = sqrt(theProblem->errorInterpolationL2);
    theProblem->errorInterpolationH1 = sqrt(theProblem->errorInterpolationH1);
}

femMesh* femMeshCreateBasicSquare(int n) {
    femMesh* theMesh = malloc(sizeof(femMesh));

    theMesh->nNode = (n + 1) * (n + 1);
    theMesh->X = malloc(sizeof(double) * theMesh->nNode);
    theMesh->Y = malloc(sizeof(double) * theMesh->nNode);

    double interval = 1.0 / n;
    int ieNode = 0;
    int i, j;

    for (i = 0; i < n + 1; i++) {

        double x = i * interval;

        for (j = 0; j < n + 1; j++) {

            double y = j * interval;
            theMesh->X[ieNode] = x;
            theMesh->Y[ieNode] = y;
            ieNode++;
        }
    }
    theMesh->nElem = 2 * n * n;
    theMesh->nLocalNode = 3;
    theMesh->elem = malloc(sizeof(int) * 3 * theMesh->nElem);
    int pos = 0;
    int comp;
    for (i = 0; i < n; i++) {

        for (j = 0; j < n; j++) {

            comp = j + i * (n + 1);
            theMesh->elem[pos * 3] = comp;
            theMesh->elem[pos * 3 + 1] = comp + n + 2;
            theMesh->elem[pos * 3 + 2] = comp + 1;
            theMesh->elem[pos * 3 + 3] = comp + n + 1;
            theMesh->elem[pos * 3 + 4] = comp + n + 2;
            theMesh->elem[pos * 3 + 5] = comp;
            pos += 2;
        }
    }

    theMesh->number = malloc(sizeof(int) * theMesh->nNode);
    for (i = 0; i < theMesh->nNode; i++)
        theMesh->number[i] = i;
    return theMesh;

}