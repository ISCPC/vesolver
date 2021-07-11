/*
 * PluginAPI.h
 *
 *  Created on: Jun 7, 2021
 *      Author: uno
 */

#ifndef COMMON_PLUGINAPI_H_
#define COMMON_PLUGINAPI_H_

typedef struct SolverPlugin SolverPlugin_t;

struct SolverPlugin {
	int (*set_option)();
	Matrix_t* (*solve_pre)(const Matrix_t* A);
	int (*solve)(const Matrix_t *A, const double* b, double* x, double res);
	int (*solve_post)(Matrix_t* A);
	void (*free)(SolverPlugin_t* solver);
};

//
// Supported Solvers
//
SolverPlugin_t* bicgstab_init();

#endif /* COMMON_PLUGINAPI_H_ */
