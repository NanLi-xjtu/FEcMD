#ifndef NPT_H
#define NPT_H

#include "main.h"

#define PCR4_NPT(r, ro, v, a, a1, a2) \
	r = ro + deltaT * v + wr * (cr[0] * a + cr[1] * a1 + cr[2] * a2)
#define PCV4_NPT(r, ro, v, a, a1, a2) \
	v = (r - ro) / deltaT + wv * (cv[0] * a + cv[1] * a1 + cv[2] * a2)

extern double extPressure, g1Sum, g2Sum, massS, massV, varS, varSa,
       varSa1, varSa2, varSo, varSv, varSvo, varV, varVa, varVa1,
       varVa2, varVo, varVv, varVvo;
extern int maxEdgeCells;

void ComputeDerivsPT ();
void PredictorStepPT ();
void CorrectorStepPT ();
void InitFeedbackVars ();
void ScaleCoords ();
void UnScaleCoords ();
void ScaleVels ();
void UnScaleVels ();
void UpdateCellSize ();

#endif
