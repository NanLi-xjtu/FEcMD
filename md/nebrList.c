#include "nebrList.h"
int nebrListFlag = 0;

void BuildNebrList ()
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset, nebri, nebrj;

	rrNebr = Sqr (rCut + rNebrShell);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
	DO_MOL{
		VSAdd (rs, mol[n].r, 0.5, region);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;

		if (c >= (VProd(cells) + nMol) || c < 0) {
			printf ("error(nebrList.c): cellList is not big enough(surface.c),\n \
				cc:(%d %d %d), cells(%d %d %d)\n \
				rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
				rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
				region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10);
			PrintMovie ();
			exit (1);
		}

		cellList[n] = cellList[c];
		cellList[c] = n;
	}
	nebrTabLen = 0;
	DO_MOL mol[n].Nnebr = 0;
	for (m1z = 0; m1z < cells.z; m1z ++) {
		for (m1y = 0; m1y < cells.y; m1y ++){
			for (m1x = 0; m1x < cells.x; m1x ++){
				VSet (m1v, m1x, m1y, m1z);
				m1 = VLinear (m1v, cells) + nMol;
				for (offset = 0; offset < N_OFFSET; offset ++){
					VAdd (m2v, m1v, vOff[offset]);
					VZero (shift);

					if (strcmp(boundaryCond, "p") == 0) {
						VCellWrapAll ();
					} else if (strcmp(boundaryCond, "n") == 0) {
						//x,y,z boundaries are no longer periodic
						if (m2v.z < 0 || m2v.z >= cells.z) continue;
						if (m2v.y < 0 || m2v.y >= cells.y) continue;
						if (m2v.x < 0 || m2v.x >= cells.x) continue;
					}

					m2 = VLinear (m2v, cells) + nMol;
					DO_CELL (j1, m1){
						DO_CELL (j2,m2){
							if (m1 != m2 || j2 < j1){
								VSub (dr, mol[j1].r, mol[j2].r);

								if (strcmp(boundaryCond, "p") == 0) VVSub (dr, shift);

								if (VLenSq (dr) < rrNebr){
									if (nebrTabLen >= nebrTabMax)
										ErrExit (ERR_TOO_MANY_NEBRS);
									nebrTab[2 * nebrTabLen] = j1;
									nebrTab[2 * nebrTabLen + 1] = j2;
									++ nebrTabLen;

									nebri = mol[j1].Nnebr;
									mol[j1].id_nebr[nebri] = j2;
									mol[j1].Nnebr ++;
									nebrj = mol[j2].Nnebr;
									mol[j2].id_nebr[nebrj] = j1;
									mol[j2].Nnebr ++;
								}
							}
						}
					}
				}
			}
		}
	}
	nebrListFlag = 1;
}

void BuildNebrListCommon ()
{
	int j1, j2, nebri, nebrj, n;
	VecR dr;
	double rrNebr;

	nebrTabLen = 0;
	rrNebr = Sqr (rCut + rNebrShell);
        DO_MOL mol[n].Nnebr = 0;
	for (j1 = 0; j1 < nMol-1; j1 ++) {
		for (j2 = j1+1; j2 < nMol; j2 ++) {
			VSub (dr, mol[j1].r, mol[j2].r);
			if (strcmp(boundaryCond, "p") == 0) VWrapAll (dr);
			if (VLenSq (dr) < rrNebr) {
				if (nebrTabLen >= nebrTabMax)
					ErrExit (ERR_TOO_MANY_NEBRS);
				nebrTab[2 * nebrTabLen] = j1;
				nebrTab[2 * nebrTabLen + 1] = j2;
				++ nebrTabLen;

				nebri = mol[j1].Nnebr;
				mol[j1].id_nebr[nebri] = j2;
				mol[j1].Nnebr ++;
				nebrj = mol[j2].Nnebr;
				mol[j2].id_nebr[nebrj] = j1;
				mol[j2].Nnebr ++;
			}
		}
	}
}

void BuildCellList ()
{
	VecR dr, invWid, rs, shift;
	VecI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	double rrNebr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset, nebri, nebrj;

	rrNebr = Sqr (rCut);
	VDiv (invWid, cells, region);
	for (n = nMol; n < nMol + VProd (cells); n ++) cellList[n] = -1;
	DO_MOL{
		VSAdd (rs, mol[n].r, 0.5, region);
		VMul (cc, rs, invWid);
		c = VLinear (cc, cells) + nMol;

		if (c >= (VProd(cells) + nMol) || c < 0) {
			printf ("error(nebrList.c): cellList is not big enough(surface.c),\n \
				cc:(%d %d %d), cells(%d %d %d)\n \
				rs(%f %f %f) A, region(%f %f %f) A\n", cc.x, cc.y, cc.z, cells.x, cells.y, cells.z, \
				rs.x * lUnit * 1.e10, rs.y * lUnit * 1.e10, rs.z * lUnit * 1.e10, \
				region.x * lUnit * 1.e10, region.y * lUnit * 1.e10, region.z * lUnit * 1.e10);
			PrintMovie ();
			exit (1);
		}

		cellList[n] = cellList[c];
		cellList[c] = n;
	}
	nebrTabLen = 0;
	DO_MOL mol[n].Nnebr = 0;
	for (m1z = 0; m1z < cells.z; m1z ++){
		for (m1y = 0; m1y < cells.y; m1y ++){
			for (m1x = 0; m1x < cells.x; m1x ++){
				VSet (m1v, m1x, m1y, m1z);
				m1 = VLinear (m1v, cells) + nMol;
				for (offset = 0; offset < N_OFFSET; offset ++){
					VAdd (m2v, m1v, vOff[offset]);
					VZero (shift);

					if (strcmp(boundaryCond, "p") == 0) {
						VCellWrapAll ();
					} else if (strcmp(boundaryCond, "n") == 0) {
						//x,y,z boundaries are no longer periodic
						if (m2v.z < 0 || m2v.z >= cells.z) continue;
						if (m2v.y < 0 || m2v.y >= cells.y) continue;
						if (m2v.x < 0 || m2v.x >= cells.x) continue;
					}

					m2 = VLinear (m2v, cells) + nMol;
					DO_CELL (j1, m1){
						DO_CELL (j2,m2){
							if (m1 != m2 || j2 < j1){
								VSub (dr, mol[j1].r, mol[j2].r);

								if (strcmp(boundaryCond, "p") == 0) VVSub (dr, shift);

								if (VLenSq (dr) < rrNebr){
									if (nebrTabLen >= nebrTabMax)
										ErrExit (ERR_TOO_MANY_NEBRS);
									nebrTab[2 * nebrTabLen] = j1;
									nebrTab[2 * nebrTabLen + 1] = j2;
									++ nebrTabLen;

									nebri = mol[j1].Nnebr;
									mol[j1].id_nebr[nebri] = j2;
									mol[j1].Nnebr ++;
									nebrj = mol[j2].Nnebr;
									mol[j2].id_nebr[nebrj] = j1;
									mol[j2].Nnebr ++;
								}
							}
						}
					}
				}
			}
		}
	}
}
