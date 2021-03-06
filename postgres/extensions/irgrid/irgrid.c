#include "IrregularGrid2D.h"

#include "postgres.h"
#include "catalog/pg_type.h"
#include "fmgr.h"
#include "utils/builtins.h"
#include "utils/array.h"

#include "liblwgeom.h" /* for GSERIALIZED */

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(createIrgrid);
PG_FUNCTION_INFO_V1(getCellnumbers);

	Datum
createIrgrid(PG_FUNCTION_ARGS)
{
	int nbrRow = PG_GETARG_INT32(1);
	int nbrCell = PG_GETARG_INT32(2);

	char *conn_str = text_to_cstring(PG_GETARG_TEXT_PP(3));
	char *rel_str = text_to_cstring(PG_GETARG_TEXT_PP(4));

	GSERIALIZED *pg_geom = (GSERIALIZED*)PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
	const GBOX *gbox = NULL;
	void *temp_irgrid2d = NULL;
	bool res = false;

	if (pg_geom != NULL) {
		LWGEOM *lwgeom = lwgeom_from_gserialized(pg_geom);
		if (lwgeom != NULL && lwgeom -> type == POLYGONTYPE) {

			gbox = lwgeom_get_bbox(lwgeom);
			if (gbox == NULL) {
				elog(ERROR, "GeometryType/POLYGONTYPE --> GBOX conversion error");
				PG_RETURN_NULL();
			}

			/*elog(INFO, "value %e", gbox->xmin);
			  elog(INFO, "value %e", gbox->ymin);
			  elog(INFO, "value %e", gbox->xmax);
			  elog(INFO, "value %e", gbox->ymax);*/

			/* Clean up memory */
			/* lwfree(lwgeom);
			   PG_FREE_IF_COPY(pg_geom, 0); */
		} else {
			elog(ERROR, "first parameter must be of POLYGONTYPE");
			PG_RETURN_NULL();
		}
	} else {
		elog(ERROR, "first parameter must be of GeometryType");
		PG_RETURN_NULL();
	}

	temp_irgrid2d = initIrgrid2D(gbox->xmin, gbox->ymin, gbox->xmax,
			gbox->ymax, nbrRow, nbrCell, conn_str, rel_str);
	res = createIrgrid2D(temp_irgrid2d);

	PG_RETURN_INT32(res);
}

	Datum
getCellnumbers(PG_FUNCTION_ARGS)
{
	GSERIALIZED *pg_geom = (GSERIALIZED*)PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
	char *conn_str = text_to_cstring(PG_GETARG_TEXT_PP(1));

	const GBOX *gbox = NULL;
	void *temp_irgrid2d = NULL;
	ArrayType *array = NULL;
	int* res = NULL;
	long unsigned int res_length = 0;
	ARInfo arInfo;

	if (pg_geom != NULL) {
		LWGEOM *lwgeom = lwgeom_from_gserialized(pg_geom);
		if (lwgeom != NULL && lwgeom -> type == POLYGONTYPE) {

			gbox = lwgeom_get_bbox(lwgeom);
			if (gbox == NULL) {
				elog(ERROR, "GeometryType/POLYGONTYPE --> GBOX conversion error");
				PG_RETURN_NULL();
			}
		} else {
			elog(ERROR, "first parameter must be of POLYGONTYPE");
			PG_RETURN_NULL();
		}
	} else {
		elog(ERROR, "first parameter must be of GeometryType");
		PG_RETURN_NULL();
	}

	temp_irgrid2d = initIrgrid2D(-1, -1, -1, -1, -1, -1,
			conn_str, NULL);

	arInfo = getCellnos(temp_irgrid2d, gbox->xmin, gbox->ymin,
			gbox->xmax, gbox->ymax, conn_str);

	res = arInfo.res;
	res_length = arInfo.res_size;
	//elog(INFO, "res_length %i", (int)res_length);

	if (res != NULL) {
		Datum elements[res_length];
		long unsigned int i;
		for (i = 0; i < res_length; i++) {
			elements[i] = Int64GetDatum(res[i]);
		}
		array = construct_array(elements, res_length, INT8OID, 8, true, 'd');
	} else {
		Datum elements[0];
		array = construct_array(elements, 0, INT8OID, 8, true, 'd');
	}

	PG_RETURN_POINTER(array);
}
