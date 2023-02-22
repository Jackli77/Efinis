/*
 *  glfem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2022 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.1.2)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include <stdlib.h>
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include "fem.h"

static double GLFEM_BLACK[4] = {0.0,0.0,0.0,1.0};
static double GLFEM_BLUE[4]  = {0.0,0.0,1.0,1.0};
static double GLFEM_RED[4]   = {1.0,0.0,0.0,1.0};


void          glfemWindowCreate(const char *windowName,int w,int h); 
void          glfemWindowFree(); 
void          glfemWindowUpdate();
int           glfemWindowShouldClose();
void          glfemReshape(femMesh *theMesh);
void          glfemPlotField(femMesh *theMesh, double *u);

void          glfemSetColor(double color[4]);
void          glfemSetTextColor(double color[4]);
void          glfemDrawMessage(char *message, double pos[2]);
void          glfemDrawNodes(double *x, double *y, int n);
void          glfemDrawElement(double *x, double *y, int n);
void          glfemDrawSolution(double *x, double *y, double* u, int n);


static void   glfemKeyCallback(GLFWwindow* self,int key,int scancode,int action,int mods);


#endif