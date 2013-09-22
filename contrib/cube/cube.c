/******************************************************************************
  contrib/cube/cube.c

  This file contains routines that can be bound to a Postgres backend and
  called by the backend in the process of processing queries.  The calling
  format for these routines is dictated by Postgres architecture.
******************************************************************************/

#include "postgres.h"

#include <float.h>
#include <math.h>

#include "access/gist.h"
#include "access/skey.h"
#include "utils/array.h"
#include "utils/builtins.h"

#include "cubedata.h"

PG_MODULE_MAGIC;

/*
 * Taken from the intarray contrib header
 */
#define ARRPTR(x)  ( (double *) ARR_DATA_PTR(x) )
#define ARRNELEMS(x)  ArrayGetNItems( ARR_NDIM(x), ARR_DIMS(x))

extern int	cube_yyparse(NDBOX **result);
extern void cube_yyerror(NDBOX **result, const char *message);
extern void cube_scanner_init(const char *str);
extern void cube_scanner_finish(void);

/*
** Input/Output routines
*/
PG_FUNCTION_INFO_V1(cube_in);
PG_FUNCTION_INFO_V1(cube_a_f8_f8);
PG_FUNCTION_INFO_V1(cube_a_f8);
PG_FUNCTION_INFO_V1(cube_out);
PG_FUNCTION_INFO_V1(cube_c_f8);
PG_FUNCTION_INFO_V1(cube_c_f8_f8);
PG_FUNCTION_INFO_V1(cube_dim);
PG_FUNCTION_INFO_V1(cube_ll_coord);
PG_FUNCTION_INFO_V1(cube_ur_coord);
PG_FUNCTION_INFO_V1(cube_subset);

Datum		cube_in(PG_FUNCTION_ARGS);
Datum		cube_a_f8_f8(PG_FUNCTION_ARGS);
Datum		cube_a_f8(PG_FUNCTION_ARGS);
Datum		cube_out(PG_FUNCTION_ARGS);
Datum		cube_c_f8(PG_FUNCTION_ARGS);
Datum		cube_c_f8_f8(PG_FUNCTION_ARGS);
Datum		cube_dim(PG_FUNCTION_ARGS);
Datum		cube_ll_coord(PG_FUNCTION_ARGS);
Datum		cube_ur_coord(PG_FUNCTION_ARGS);
Datum		cube_subset(PG_FUNCTION_ARGS);

/*
 * Input/Output for typed cubes
 */

static NDBOX* cube_arr_arr(Datum ur_datum, Datum ll_datum, int type);
static NDBOX* cube_arr(Datum ur_datum, int type);
static NDBOX* cube_num(Datum x_datum, int type);
static NDBOX* cube_num_num(Datum x0_datum, Datum x1_datum, int type);

CUBE_TYPE_WRAPPER1(cube_arr, FLOAT8);
CUBE_TYPE_WRAPPER1(cube_num, FLOAT8);
CUBE_TYPE_WRAPPER2(cube_arr_arr, FLOAT8);
CUBE_TYPE_WRAPPER2(cube_num_num, FLOAT8);

CUBE_TYPE_WRAPPER1(cube_arr, FLOAT4);
CUBE_TYPE_WRAPPER1(cube_num, FLOAT4);
CUBE_TYPE_WRAPPER2(cube_arr_arr, FLOAT4);
CUBE_TYPE_WRAPPER2(cube_num_num, FLOAT4);

CUBE_TYPE_WRAPPER1(cube_arr, INT4);
CUBE_TYPE_WRAPPER1(cube_num, INT4);
CUBE_TYPE_WRAPPER2(cube_arr_arr, INT4);
CUBE_TYPE_WRAPPER2(cube_num_num, INT4);

CUBE_TYPE_WRAPPER1(cube_arr, INT2);
CUBE_TYPE_WRAPPER1(cube_num, INT2);
CUBE_TYPE_WRAPPER2(cube_arr_arr, INT2);
CUBE_TYPE_WRAPPER2(cube_num_num, INT2);

CUBE_TYPE_WRAPPER1(cube_arr, INT1);
CUBE_TYPE_WRAPPER1(cube_num, INT1);
CUBE_TYPE_WRAPPER2(cube_arr_arr, INT1);
CUBE_TYPE_WRAPPER2(cube_num_num, INT1);

/*
** GiST support methods
*/
PG_FUNCTION_INFO_V1(g_cube_consistent);
PG_FUNCTION_INFO_V1(g_cube_compress);
PG_FUNCTION_INFO_V1(g_cube_decompress);
PG_FUNCTION_INFO_V1(g_cube_penalty);
PG_FUNCTION_INFO_V1(g_cube_picksplit);
PG_FUNCTION_INFO_V1(g_cube_union);
PG_FUNCTION_INFO_V1(g_cube_same);

Datum		g_cube_consistent(PG_FUNCTION_ARGS);
Datum		g_cube_compress(PG_FUNCTION_ARGS);
Datum		g_cube_decompress(PG_FUNCTION_ARGS);
Datum		g_cube_penalty(PG_FUNCTION_ARGS);
Datum		g_cube_picksplit(PG_FUNCTION_ARGS);
Datum		g_cube_union(PG_FUNCTION_ARGS);
Datum		g_cube_same(PG_FUNCTION_ARGS);

/*
** B-tree support functions
*/
PG_FUNCTION_INFO_V1(cube_eq);
PG_FUNCTION_INFO_V1(cube_ne);
PG_FUNCTION_INFO_V1(cube_lt);
PG_FUNCTION_INFO_V1(cube_gt);
PG_FUNCTION_INFO_V1(cube_le);
PG_FUNCTION_INFO_V1(cube_ge);
PG_FUNCTION_INFO_V1(cube_cmp);

Datum		cube_eq(PG_FUNCTION_ARGS);
Datum		cube_ne(PG_FUNCTION_ARGS);
Datum		cube_lt(PG_FUNCTION_ARGS);
Datum		cube_gt(PG_FUNCTION_ARGS);
Datum		cube_le(PG_FUNCTION_ARGS);
Datum		cube_ge(PG_FUNCTION_ARGS);
Datum		cube_cmp(PG_FUNCTION_ARGS);

/*
** R-tree support functions
*/
PG_FUNCTION_INFO_V1(cube_contains);
PG_FUNCTION_INFO_V1(cube_contained);
PG_FUNCTION_INFO_V1(cube_overlap);
PG_FUNCTION_INFO_V1(cube_union);
PG_FUNCTION_INFO_V1(cube_inter);
PG_FUNCTION_INFO_V1(cube_size);

Datum		cube_contains(PG_FUNCTION_ARGS);
Datum		cube_contained(PG_FUNCTION_ARGS);
Datum		cube_overlap(PG_FUNCTION_ARGS);
Datum		cube_union(PG_FUNCTION_ARGS);
Datum		cube_inter(PG_FUNCTION_ARGS);
Datum		cube_size(PG_FUNCTION_ARGS);

/*
** miscellaneous
*/
PG_FUNCTION_INFO_V1(cube_distance);
PG_FUNCTION_INFO_V1(cube_is_point);
PG_FUNCTION_INFO_V1(cube_enlarge);

Datum		cube_distance(PG_FUNCTION_ARGS);
Datum		cube_is_point(PG_FUNCTION_ARGS);
Datum		cube_enlarge(PG_FUNCTION_ARGS);

/*
** For internal use only
*/
void		check_coord(double num, int type);
int32		cube_cmp_v0(NDBOX *a, NDBOX *b);
bool		cube_contains_v0(NDBOX *a, NDBOX *b);
bool		cube_overlap_v0(NDBOX *a, NDBOX *b);
NDBOX	   *cube_union_v0(NDBOX *a, NDBOX *b);
void		rt_cube_size(NDBOX *a, double *sz);
NDBOX	   *g_cube_binary_union(NDBOX *r1, NDBOX *r2, int *sizep);
bool		g_cube_leaf_consistent(NDBOX *key, NDBOX *query, StrategyNumber strategy);
bool		g_cube_internal_consistent(NDBOX *key, NDBOX *query, StrategyNumber strategy);

/*
** Auxiliary funxtions
*/
static double distance_1D(double a1, double a2, double b1, double b2);


/*****************************************************************************
 * Typed cube abstraction layer
 *****************************************************************************/
double get_coord(NDBOX *cube, int i)
{
	switch (TYPE(cube))
	{
		case CUBE_FLOAT4:
			return cube->coord_f4[i];
		case CUBE_INT4:
			return cube->coord_i4[i];
		case CUBE_INT2:
			return cube->coord_i2[i];
		case CUBE_INT1:
			return cube->coord_i1[i];
		case CUBE_FLOAT8:
		default:
			return cube->coord_f8[i];
	}
}

void set_coord(NDBOX *cube, int i, double value)
{
	switch (TYPE(cube))
	{
		case CUBE_FLOAT4:
			cube->coord_f4[i] = value;
			break;
		case CUBE_INT4:
			cube->coord_i4[i] = value;
			break;
		case CUBE_INT2:
			cube->coord_i2[i] = value;
			break;
		case CUBE_INT1:
			cube->coord_i1[i] = value;
			break;
		case CUBE_FLOAT8:
		default:
			cube->coord_f8[i] = value;
			break;
	}
}

NDBOX* init_cube(int dim, int point, int type)
{
	NDBOX	   *cube;
	int			size;

	switch (type)
	{
		case CUBE_FLOAT4:
			size = offsetof(NDBOX, coord_f8[0]) + dim*(!point + 1)*sizeof(float);
			break;
		case CUBE_INT4:
			size = offsetof(NDBOX, coord_f8[0]) + dim*(!point + 1)*sizeof(int);
			break;
		case CUBE_INT2:
			size = offsetof(NDBOX, coord_f8[0]) + dim*(!point + 1)*sizeof(short);
			break;
		case CUBE_INT1:
			size = offsetof(NDBOX, coord_f8[0]) + dim*(!point + 1)*sizeof(char);
			break;
		case CUBE_FLOAT8:
		default:
			size = offsetof(NDBOX, coord_f8[0]) + dim*(!point + 1)*sizeof(double);
			break;
	}

	cube = (NDBOX *) palloc0(size);
	SET_VARSIZE(cube, size);
	SET_DIM(cube, dim);
	SET_TYPE(cube, type);
	if (point)
		SET_POINT_BIT(cube);

	return cube;
}

void cube_to_point(NDBOX *cube)
{
	int new_size;

	new_size = offsetof(NDBOX, coord_f8[0]) +
		(VARSIZE(cube) - offsetof(NDBOX, coord_f8[0]))/2;

	cube = repalloc(cube, new_size);
	SET_VARSIZE(cube, new_size);
	SET_POINT_BIT(cube);
}

void check_coord(double num, int type)
{
	switch (type)
	{
		case CUBE_INT4:
			if ((num < INT_MIN) || (num > INT_MAX))
				ereport(ERROR,
					(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
					 errmsg("Cube coordinate out of requested type range"),
					 errdetail("Value (%i) out of signed int4 type", (int)num)));
			break;
		case CUBE_INT2:
			if ((num < SHRT_MIN) || (num > SHRT_MAX))
				ereport(ERROR,
					(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
					 errmsg("Cube coordinate out of requested type range"),
					 errdetail("Value (%i) out of signed int2 type", (int)num)));
			break;
		case CUBE_INT1:
			if ( (num < SCHAR_MIN) || (num > SCHAR_MAX))
				ereport(ERROR,
					(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
					 errmsg("Cube coordinate out of requested type range"),
					 errdetail("Value (%i) out of signed int1 type", (int)num)));
			break;
	}
}


/*****************************************************************************
 * Input/Output functions
 *****************************************************************************/

/* NdBox = [(lowerleft),(upperright)] */
/* [(xLL(1)...xLL(N)),(xUR(1)...xUR(n))] */
Datum
cube_in(PG_FUNCTION_ARGS)
{
	char	   *str = PG_GETARG_CSTRING(0);
	NDBOX	   *result;

	// printf("cube_in for %s\n", str);
	cube_scanner_init(str);

	if (cube_yyparse(&result) != 0)
		cube_yyerror(&result, "bogus input");

	cube_scanner_finish();

	PG_RETURN_NDBOX(result);
}

/*
** Allows the construction of a cube from 2 float[]'s
*/
static NDBOX*
cube_arr_arr(Datum ur_datum, Datum ll_datum, int type)
{
	ArrayType  *ur = DatumGetArrayTypeP(ur_datum);
	ArrayType  *ll = DatumGetArrayTypeP(ll_datum);
	NDBOX	   *result;
	int			i, dim;
	bool		point = true;
	double	   *dur, *dll;

	if (array_contains_nulls(ur) || array_contains_nulls(ll))
		ereport(ERROR,
				(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
				 errmsg("cannot work with arrays containing NULLs")));

	dim = ARRNELEMS(ur);
	if (ARRNELEMS(ll) != dim)
		ereport(ERROR,
				(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
				 errmsg("UR and LL arrays must be of same length")));

	dur = ARRPTR(ur);
	dll = ARRPTR(ll);
	result = init_cube(dim, 0, type);

	for (i = 0; i < dim; i++)
	{
		check_coord(dur[i], type);
		set_coord(result, i, dur[i]);
		check_coord(dll[i], type);
		set_coord(result, i+dim, dll[i]);
		if (dur[i] != dll[i])
			point = false;
	}

	if (point)
		cube_to_point(result);

	return result;
}


/*
** Allows the construction of a zero-volume cube from a float[]
*/
static NDBOX*
cube_arr(Datum ur_datum, int type)
{
	ArrayType  *ur = DatumGetArrayTypeP(ur_datum);
	NDBOX	   *result;
	int			i, dim;
	double	   *dur;

	if (array_contains_nulls(ur))
		ereport(ERROR,
				(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
				 errmsg("cannot work with arrays containing NULLs")));

	dim = ARRNELEMS(ur);
	dur = ARRPTR(ur);

	result = init_cube(dim, 1, type);

	for (i = 0; i < dim; i++){
		check_coord(dur[i], type);
		set_coord(result, i, dur[i]);
	}

	return result;
}

/* Create a one dimensional box with identical upper and lower coordinates */
static NDBOX*
cube_num(Datum x_datum, int type)
{
	double		x = DatumGetFloat8(x_datum);
	NDBOX	   *result;

	result = init_cube(1, 1, type);
	check_coord(x, type);
	set_coord(result, 0, x);

	return result;
}

/* Create a one dimensional box */
static NDBOX*
cube_num_num(Datum x0_datum, Datum x1_datum, int type)
{
	double		x0 = DatumGetFloat8(x0_datum);
	double		x1 = DatumGetFloat8(x1_datum);
	NDBOX		*result;

	if (x0 == x1)
	{
		result = init_cube(1, 1, type);
		check_coord(x0, type);
		set_coord(result, 0, x0);
	}
	else
	{
		result = init_cube(1, 0, type);
		check_coord(x0, type);
		set_coord(result, 0, x0);
		check_coord(x1, type);
		set_coord(result, 1, x1);
	}

	return result;
}

Datum
cube_subset(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	ArrayType  *idx = PG_GETARG_ARRAYTYPE_P(1);
	NDBOX	   *result;
	int			dim, i;
	int		   *dx;

	if (array_contains_nulls(idx))
		ereport(ERROR,
				(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
				 errmsg("cannot work with arrays containing NULLs")));

	dx = (int32 *) ARR_DATA_PTR(idx);

	dim = ARRNELEMS(idx);
	result = init_cube(dim, IS_POINT(c), TYPE(c));

	for (i = 0; i < dim; i++)
	{
		if ((dx[i] <= 0) || (dx[i] > DIM(c)))
		{
			pfree(result);
			ereport(ERROR,
					(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
					 errmsg("Index out of bounds")));
		}
		set_coord(result, i, LL_COORD(c, dx[i]-1));
		if (!IS_POINT(c))
			set_coord(result, i+dim, UR_COORD(c, dx[i]-1));
	}

	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_NDBOX(result);
}


Datum
cube_out(PG_FUNCTION_ARGS)
{
	NDBOX	   *cube = PG_GETARG_NDBOX(0);
	StringInfoData buf;
	int			dim = DIM(cube);
	int			i;
	int			ndig;

	initStringInfo(&buf);

	/*
	 * Get the number of digits to display.
	 */
	ndig = DBL_DIG + extra_float_digits;
	if (ndig < 1)
		ndig = 1;

	/*
	 * while printing the first (LL) corner, check if it is equal to the
	 * second one
	 */
	appendStringInfoChar(&buf, '(');
	for (i = 0; i < dim; i++)
	{
		if (i > 0)
			appendStringInfo(&buf, ", ");
		appendStringInfo(&buf, "%.*g", ndig, LL_COORD(cube,i));
	}
	appendStringInfoChar(&buf, ')');

	if (!IS_POINT(cube))
	{
		appendStringInfo(&buf, ",(");
		for (i = 0; i < dim; i++)
		{
			if (i > 0)
				appendStringInfo(&buf, ", ");
			appendStringInfo(&buf, "%.*g", ndig, UR_COORD(cube, i));
		}
		appendStringInfoChar(&buf, ')');
	}

	/*
	 * print type unless it is float8
	 */
	switch (TYPE(cube))
	{
		case CUBE_FLOAT4:
			appendStringInfo(&buf, ":f4");
			break;
		case CUBE_INT4:
			appendStringInfo(&buf, ":i4");
			break;
		case CUBE_INT2:
			appendStringInfo(&buf, ":i2");
			break;
		case CUBE_INT1:
			appendStringInfo(&buf, ":i1");
			break;
	}

	PG_FREE_IF_COPY(cube, 0);
	PG_RETURN_CSTRING(buf.data);
}


/*****************************************************************************
 *						   GiST functions
 *****************************************************************************/

/*
** The GiST Consistent method for boxes
** Should return false if for all data items x below entry,
** the predicate x op query == FALSE, where op is the oper
** corresponding to strategy in the pg_amop table.
*/
Datum
g_cube_consistent(PG_FUNCTION_ARGS)
{
	GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
	NDBOX	   *query = PG_GETARG_NDBOX(1);
	StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);

	/* Oid		subtype = PG_GETARG_OID(3); */
	bool	   *recheck = (bool *) PG_GETARG_POINTER(4);
	bool		res;

	/* All cases served by this function are exact */
	*recheck = false;

	/*
	 * if entry is not leaf, use g_cube_internal_consistent, else use
	 * g_cube_leaf_consistent
	 */
	if (GIST_LEAF(entry))
		res = g_cube_leaf_consistent(DatumGetNDBOX(entry->key),
									 query, strategy);
	else
		res = g_cube_internal_consistent(DatumGetNDBOX(entry->key),
										 query, strategy);

	PG_FREE_IF_COPY(query, 1);
	PG_RETURN_BOOL(res);
}


/*
** The GiST Union method for boxes
** returns the minimal bounding box that encloses all the entries in entryvec
*/
Datum
g_cube_union(PG_FUNCTION_ARGS)
{
	GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
	int		   *sizep = (int *) PG_GETARG_POINTER(1);
	NDBOX	   *out = (NDBOX *) NULL;
	NDBOX	   *tmp;
	int			i;

	tmp = DatumGetNDBOX(entryvec->vector[0].key);

	/*
	 * sizep = sizeof(NDBOX); -- NDBOX has variable size
	 */
	*sizep = VARSIZE(tmp);

	for (i = 1; i < entryvec->n; i++)
	{
		out = g_cube_binary_union(tmp,
								  DatumGetNDBOX(entryvec->vector[i].key),
								  sizep);
		tmp = out;
	}

	PG_RETURN_POINTER(out);
}


/*
** GiST Compress and Decompress methods for boxes
** do not do anything.
*/
Datum
g_cube_compress(PG_FUNCTION_ARGS)
{
	GISTENTRY			*entry = (GISTENTRY *) PG_GETARG_POINTER(0);
	GISTENTRY			*retval;
	NDBOX 				*cube = DatumGetNDBOX(PG_DETOAST_DATUM(entry->key));

	retval = palloc(sizeof(GISTENTRY));
	gistentryinit(*retval, PointerGetDatum(cube),
					entry->rel, entry->page, entry->offset, FALSE);
	PG_RETURN_POINTER(retval);
}

Datum
g_cube_decompress(PG_FUNCTION_ARGS)
{
	GISTENTRY		*entry = (GISTENTRY *) PG_GETARG_POINTER(0);
	GISTENTRY		*retval;
	NDBOX			*cube = DatumGetNDBOX(PG_DETOAST_DATUM(entry->key));

	retval = palloc(sizeof(GISTENTRY));
	gistentryinit(*retval, PointerGetDatum(cube),
				entry->rel, entry->page, entry->offset, FALSE);
	PG_RETURN_POINTER(retval);
}


/*
** The GiST Penalty method for boxes
** As in the R-tree paper, we use change in area as our penalty metric
*/
Datum
g_cube_penalty(PG_FUNCTION_ARGS)
{
	GISTENTRY  *origentry = (GISTENTRY *) PG_GETARG_POINTER(0);
	GISTENTRY  *newentry = (GISTENTRY *) PG_GETARG_POINTER(1);
	float	   *result = (float *) PG_GETARG_POINTER(2);
	NDBOX	   *ud;
	double		tmp1,
				tmp2;

	ud = cube_union_v0(DatumGetNDBOX(origentry->key),
					   DatumGetNDBOX(newentry->key));
	rt_cube_size(ud, &tmp1);
	rt_cube_size(DatumGetNDBOX(origentry->key), &tmp2);
	*result = (float) (tmp1 - tmp2);

	PG_RETURN_FLOAT8(*result);
}


/*
** The GiST PickSplit method for boxes
** We use Guttman's poly time split algorithm
*/
Datum
g_cube_picksplit(PG_FUNCTION_ARGS)
{
	GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
	GIST_SPLITVEC *v = (GIST_SPLITVEC *) PG_GETARG_POINTER(1);
	OffsetNumber i,
				j;
	NDBOX	   *datum_alpha,
			   *datum_beta;
	NDBOX	   *datum_l,
			   *datum_r;
	NDBOX	   *union_d,
			   *union_dl,
			   *union_dr;
	NDBOX	   *inter_d;
	bool		firsttime;
	double		size_alpha,
				size_beta,
				size_union,
				size_inter;
	double		size_waste,
				waste;
	double		size_l,
				size_r;
	int			nbytes;
	OffsetNumber seed_1 = 1,
				seed_2 = 2;
	OffsetNumber *left,
			   *right;
	OffsetNumber maxoff;
	 
	maxoff = entryvec->n - 2;
	nbytes = (maxoff + 2) * sizeof(OffsetNumber);
	v->spl_left = (OffsetNumber *) palloc(nbytes);
	v->spl_right = (OffsetNumber *) palloc(nbytes);

	firsttime = true;
	waste = 0.0;

	for (i = FirstOffsetNumber; i < maxoff; i = OffsetNumberNext(i))
	{
		datum_alpha = DatumGetNDBOX(entryvec->vector[i].key);
		for (j = OffsetNumberNext(i); j <= maxoff; j = OffsetNumberNext(j))
		{
			datum_beta = DatumGetNDBOX(entryvec->vector[j].key);

			/* compute the wasted space by unioning these guys */
			/* size_waste = size_union - size_inter; */
			union_d = cube_union_v0(datum_alpha, datum_beta);
			rt_cube_size(union_d, &size_union);
			inter_d = DatumGetNDBOX(DirectFunctionCall2(cube_inter,
						  entryvec->vector[i].key, entryvec->vector[j].key));
			rt_cube_size(inter_d, &size_inter);
			size_waste = size_union - size_inter;

			/*
			 * are these a more promising split than what we've already seen?
			 */

			if (size_waste > waste || firsttime)
			{
				waste = size_waste;
				seed_1 = i;
				seed_2 = j;
				firsttime = false;
			}
		}
	}

	left = v->spl_left;
	v->spl_nleft = 0;
	right = v->spl_right;
	v->spl_nright = 0;

	datum_alpha = DatumGetNDBOX(entryvec->vector[seed_1].key);
	datum_l = cube_union_v0(datum_alpha, datum_alpha);
	rt_cube_size(datum_l, &size_l);
	datum_beta = DatumGetNDBOX(entryvec->vector[seed_2].key);
	datum_r = cube_union_v0(datum_beta, datum_beta);
	rt_cube_size(datum_r, &size_r);

	/*
	 * Now split up the regions between the two seeds.	An important property
	 * of this split algorithm is that the split vector v has the indices of
	 * items to be split in order in its left and right vectors.  We exploit
	 * this property by doing a merge in the code that actually splits the
	 * page.
	 *
	 * For efficiency, we also place the new index tuple in this loop. This is
	 * handled at the very end, when we have placed all the existing tuples
	 * and i == maxoff + 1.
	 */

	maxoff = OffsetNumberNext(maxoff);
	for (i = FirstOffsetNumber; i <= maxoff; i = OffsetNumberNext(i))
	{
		/*
		 * If we've already decided where to place this item, just put it on
		 * the right list.	Otherwise, we need to figure out which page needs
		 * the least enlargement in order to store the item.
		 */

		if (i == seed_1)
		{
			*left++ = i;
			v->spl_nleft++;
			continue;
		}
		else if (i == seed_2)
		{
			*right++ = i;
			v->spl_nright++;
			continue;
		}

		/* okay, which page needs least enlargement? */
		datum_alpha = DatumGetNDBOX(entryvec->vector[i].key);
		union_dl = cube_union_v0(datum_l, datum_alpha);
		union_dr = cube_union_v0(datum_r, datum_alpha);
		rt_cube_size(union_dl, &size_alpha);
		rt_cube_size(union_dr, &size_beta);

		/* pick which page to add it to */
		if (size_alpha - size_l < size_beta - size_r)
		{
			datum_l = union_dl;
			size_l = size_alpha;
			*left++ = i;
			v->spl_nleft++;
		}
		else
		{
			datum_r = union_dr;
			size_r = size_beta;
			*right++ = i;
			v->spl_nright++;
		}
	}
	*left = *right = FirstOffsetNumber; /* sentinel value, see dosplit() */

	v->spl_ldatum = PointerGetDatum(datum_l);
	v->spl_rdatum = PointerGetDatum(datum_r);

	PG_RETURN_POINTER(v);
}


/*
** Equality method
*/
Datum
g_cube_same(PG_FUNCTION_ARGS)
{
	NDBOX	   *b1 = PG_GETARG_NDBOX(0);
	NDBOX	   *b2 = PG_GETARG_NDBOX(1);
	bool	   *result = (bool *) PG_GETARG_POINTER(2);

	if (cube_cmp_v0(b1, b2) == 0)
		*result = TRUE;
	else
		*result = FALSE;

	PG_RETURN_NDBOX(result);
}

/*
** SUPPORT ROUTINES
*/
bool
g_cube_leaf_consistent(NDBOX *key,
					   NDBOX *query,
					   StrategyNumber strategy)
{
	bool		retval;

	switch (strategy)
	{
		case RTOverlapStrategyNumber:
			retval = (bool) cube_overlap_v0(key, query);
			break;
		case RTSameStrategyNumber:
			retval = (bool) (cube_cmp_v0(key, query) == 0);
			break;
		case RTContainsStrategyNumber:
		case RTOldContainsStrategyNumber:
			retval = (bool) cube_contains_v0(key, query);
			break;
		case RTContainedByStrategyNumber:
		case RTOldContainedByStrategyNumber:
			retval = (bool) cube_contains_v0(query, key);
			break;
		default:
			retval = FALSE;
	}
	return (retval);
}


bool
g_cube_internal_consistent(NDBOX *key,
						   NDBOX *query,
						   StrategyNumber strategy)
{
	bool		retval;

	switch (strategy)
	{
		case RTOverlapStrategyNumber:
			retval = (bool) cube_overlap_v0(key, query);
			break;
		case RTSameStrategyNumber:
		case RTContainsStrategyNumber:
		case RTOldContainsStrategyNumber:
			retval = (bool) cube_contains_v0(key, query);
			break;
		case RTContainedByStrategyNumber:
		case RTOldContainedByStrategyNumber:
			retval = (bool) cube_overlap_v0(key, query);
			break;
		default:
			retval = FALSE;
	}
	return (retval);
}

NDBOX *
g_cube_binary_union(NDBOX *r1, NDBOX *r2, int *sizep)
{
	NDBOX	   *retval;

	retval = cube_union_v0(r1, r2);
	*sizep = VARSIZE(retval);

	return (retval);
}


/* cube_union_v0 */
NDBOX *
cube_union_v0(NDBOX *a, NDBOX *b)
{
	int			i;
	bool		point_result = true;
	NDBOX	   *result;

	/* let's try to guess result for same pointers */
	if (a == b)
		return a;

	/* swap the box pointers if needed */
	if (DIM(a) < DIM(b))
	{
		NDBOX	   *tmp = b;
		b = a;
		a = tmp;
	}

	result = init_cube( DIM(a), 0, Min(TYPE(a), TYPE(b)) );

	 /* compute the union */
	for (i = 0; i < DIM(b); i++)
	{
		set_coord(result, i, Min(
			Min(LL_COORD(a,i), UR_COORD(a,i)),
			Min(LL_COORD(b,i), UR_COORD(b,i))
		));
		set_coord(result, i + DIM(a), Max(
			Max(LL_COORD(a,i), UR_COORD(a,i)), 
			Max(LL_COORD(b,i), UR_COORD(b,i))
		));
		if (LL_COORD(result,i) != UR_COORD(result,i))
			point_result = false;
	}
	for (i = DIM(b); i < DIM(a); i++)
	{
		set_coord(result, i, Min(0,
			Min(LL_COORD(a,i), UR_COORD(a,i))
		));
		set_coord(result, i + DIM(a), Max(0,
			Max(LL_COORD(a,i), UR_COORD(a,i))
		));
		if (LL_COORD(result,i) != UR_COORD(result,i))
			point_result = false;
	}

	if (point_result)
		cube_to_point(result);

	return (result);
}


Datum
cube_union(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	NDBOX	   *b = PG_GETARG_NDBOX(1);
	NDBOX	   *res;

	res = cube_union_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_NDBOX(res);
}


/* cube_inter */
Datum
cube_inter(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	NDBOX	   *b = PG_GETARG_NDBOX(1);
	NDBOX	   *result;
	bool		swapped = false,
				point_result = true;
	int			i;

	/* swap the box pointers if needed */
	if (DIM(a) < DIM(b))
	{
		NDBOX	   *tmp = b;

		b = a;
		a = tmp;
		swapped = true;
	}

	result = init_cube( DIM(a), 0, Min(TYPE(a), TYPE(b)) );

	for (i = 0; i < DIM(b); i++)
	{

		set_coord(result, i, Max(
			Min(LL_COORD(a,i), UR_COORD(a,i)),
			Min(LL_COORD(b,i), UR_COORD(b,i))
		));
		set_coord(result, i + DIM(a), Min(
			Max(LL_COORD(a,i), UR_COORD(a,i)),
			Max(LL_COORD(b,i), UR_COORD(b,i))
		));
		if (LL_COORD(result, i) != UR_COORD(result, i))
			point_result = false;
	}
	for (i = DIM(b); i < DIM(a); i++)
	{
		set_coord(result, i, Max(0,
			Min(LL_COORD(a,i), UR_COORD(a,i))
		));
		set_coord(result, i + DIM(a), Min(0,
			Max(LL_COORD(a,i), UR_COORD(a,i))
		));
		if (LL_COORD(result, i) != UR_COORD(result, i))
			point_result = false;
	}

	if (point_result)
		cube_to_point(result);

	if (swapped)
	{
		PG_FREE_IF_COPY(b, 0);
		PG_FREE_IF_COPY(a, 1);
	}
	else
	{
		PG_FREE_IF_COPY(a, 0);
		PG_FREE_IF_COPY(b, 1);
	}

	/*
	 * Is it OK to return a non-null intersection for non-overlapping boxes?
	 */
	PG_RETURN_NDBOX(result);
}


/* cube_size */
Datum
cube_size(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	double		result;
	int			i;

	result = 1.0;
	for (i = 0; i < DIM(a); i++)
		result = result * Abs((LL_COORD(a,i) - UR_COORD(a,i)));

	PG_FREE_IF_COPY(a, 0);
	PG_RETURN_FLOAT8(result);
}


void
rt_cube_size(NDBOX *a, double *size)
{
	int i;

	if (a == (NDBOX *) NULL)
		*size = 0.0;
	else
	{
		*size = 1.0;
		for (i = 0; i < DIM(a); i++)
			*size = (*size) * Abs(UR_COORD(a,i) - LL_COORD(a,i));
	}
	return;
}


/* make up a metric in which one box will be 'lower' than the other
   -- this can be useful for sorting and to determine uniqueness */
int32
cube_cmp_v0(NDBOX *a, NDBOX *b)
{
	int			i;
	int			dim;

	dim = Min(DIM(a), DIM(b));

	/* compare the common dimensions */
	for (i = 0; i < dim; i++)
	{
		if (Min(LL_COORD(a, i), UR_COORD(a, i)) >
			Min(LL_COORD(b, i), UR_COORD(b, i)))
			return 1;
		if (Min(LL_COORD(a, i), UR_COORD(a, i)) <
			Min(LL_COORD(b, i), UR_COORD(b, i)))
			return -1;
	}
	for (i = 0; i < dim; i++)
	{
		if (Max(LL_COORD(a, i), UR_COORD(a, i)) >
			Max(LL_COORD(b, i), UR_COORD(b, i)))
			return 1;
		if (Max(LL_COORD(a, i), UR_COORD(a, i)) <
			Max(LL_COORD(b, i), UR_COORD(b, i)))
			return -1;
	}

	/* compare extra dimensions to zero */
	if (DIM(a) > DIM(b))
	{
		for (i = dim; i < DIM(a); i++)
		{
			if (Min(LL_COORD(a, i), UR_COORD(a, i)) > 0)
				return 1;
			if (Min(LL_COORD(a, i), UR_COORD(a, i)) < 0)
				return -1;
		}
		for (i = dim; i < DIM(a); i++)
		{
			if (Max(LL_COORD(a, i), UR_COORD(a, i)) > 0)
				return 1;
			if (Max(LL_COORD(a, i), UR_COORD(a, i)) < 0)
				return -1;
		}

		/*
		 * if all common dimensions are equal, the cube with more dimensions
		 * wins
		 */
		return 1;
	}
	if (DIM(a) < DIM(b))
	{
		for (i = dim; i < DIM(b); i++)
		{
			if (Min(LL_COORD(b, i), UR_COORD(b, i)) > 0)
				return -1;
			if (Min(LL_COORD(b, i), UR_COORD(b, i)) < 0)
				return 1;
		}
		for (i = dim; i < DIM(b); i++)
		{
			if (Max(LL_COORD(b, i), UR_COORD(b, i)) > 0)
				return -1;
			if (Max(LL_COORD(b, i), UR_COORD(b, i)) < 0)
				return 1;
		}

		/*
		 * if all common dimensions are equal, the cube with more dimensions
		 * wins
		 */
		return -1;
	}

	/* They're really equal */
	return 0;
}


Datum
cube_cmp(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_INT32(res);
}


Datum
cube_eq(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res == 0);
}


Datum
cube_ne(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res != 0);
}


Datum
cube_lt(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res < 0);
}


Datum
cube_gt(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res > 0);
}


Datum
cube_le(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res <= 0);
}


Datum
cube_ge(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	int32		res;

	res = cube_cmp_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res >= 0);
}


/* Contains */
/* Box(A) CONTAINS Box(B) IFF pt(A) < pt(B) */
bool
cube_contains_v0(NDBOX *a, NDBOX *b)
{
	int			i;

	if ((a == NULL) || (b == NULL))
		return (FALSE);

	if (DIM(a) < DIM(b))
	{
		/*
		 * the further comparisons will make sense if the excess dimensions of
		 * (b) were zeroes Since both UL and UR coordinates must be zero, we
		 * can check them all without worrying about which is which.
		 */
		for (i = DIM(a); i < DIM(b); i++)
		{
			if (LL_COORD(b,i) != 0)
				return (FALSE);
			if (UR_COORD(b, i) != 0)
				return (FALSE);
		}
	}

	/* Can't care less about the excess dimensions of (a), if any */
	for (i = 0; i < Min(DIM(a), DIM(b)); i++)
	{
		if (Min(LL_COORD(a,i), UR_COORD(a, i)) >
			Min(LL_COORD(b,i), UR_COORD(b, i)))
			return (FALSE);
		if (Max(LL_COORD(a,i), UR_COORD(a, i)) <
			Max(LL_COORD(b,i), UR_COORD(b, i)))
			return (FALSE);
	}

	return (TRUE);
}


Datum
cube_contains(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	bool		res;

	res = cube_contains_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res);
}


/* Contained */
/* Box(A) Contained by Box(B) IFF Box(B) Contains Box(A) */
Datum
cube_contained(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	bool		res;

	res = cube_contains_v0(b, a);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res);
}


/* Overlap */
/* Box(A) Overlap Box(B) IFF (pt(a)LL < pt(B)UR) && (pt(b)LL < pt(a)UR) */
bool
cube_overlap_v0(NDBOX *a, NDBOX *b)
{
	int			i;

	/*
	 * This *very bad* error was found in the source: if ( (a==NULL) ||
	 * (b=NULL) ) return(FALSE);
	 */
	if ((a == NULL) || (b == NULL))
		return (FALSE);

	/* swap the box pointers if needed */
	if (DIM(a) < DIM(b))
	{
		NDBOX	   *tmp = b;

		b = a;
		a = tmp;
	}

	/* compare within the dimensions of (b) */
	for (i = 0; i < DIM(b); i++)
	{
		if (Min(LL_COORD(a,i), UR_COORD(a,i)) > Max(LL_COORD(b,i), UR_COORD(b,i)))
			return (FALSE);
		if (Max(LL_COORD(a,i), UR_COORD(a,i)) < Min(LL_COORD(b,i), UR_COORD(b,i)))
			return (FALSE);
	}

	/* compare to zero those dimensions in (a) absent in (b) */
	for (i = DIM(b); i < DIM(a); i++)
	{
		if (Min(LL_COORD(a,i), UR_COORD(a,i)) > 0)
			return (FALSE);
		if (Max(LL_COORD(a,i), UR_COORD(a,i)) < 0)
			return (FALSE);
	}

	return (TRUE);
}


Datum
cube_overlap(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	bool		res;

	res = cube_overlap_v0(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(res);
}


/* Distance */
/* The distance is computed as a per axis sum of the squared distances
   between 1D projections of the boxes onto Cartesian axes. Assuming zero
   distance between overlapping projections, this metric coincides with the
   "common sense" geometric distance */
Datum
cube_distance(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0),
			   *b = PG_GETARG_NDBOX(1);
	bool		swapped = false;
	double		d,
				distance;
	int			i;

	/* swap the box pointers if needed */
	if (DIM(a) < DIM(b))
	{
		NDBOX	   *tmp = b;

		b = a;
		a = tmp;
		swapped = true;
	}

	distance = 0.0;
	/* compute within the dimensions of (b) */
	for (i = 0; i < DIM(b); i++)
	{
		d = distance_1D(LL_COORD(a,i), UR_COORD(a,i), LL_COORD(b,i), UR_COORD(b,i));
		distance += d * d;
	}

	/* compute distance to zero for those dimensions in (a) absent in (b) */
	for (i = DIM(b); i < DIM(a); i++)
	{
		d = distance_1D(LL_COORD(a,i), UR_COORD(a,i), 0.0, 0.0);
		distance += d * d;
	}

	if (swapped)
	{
		PG_FREE_IF_COPY(b, 0);
		PG_FREE_IF_COPY(a, 1);
	}
	else
	{
		PG_FREE_IF_COPY(a, 0);
		PG_FREE_IF_COPY(b, 1);
	}

	PG_RETURN_FLOAT8(sqrt(distance));
}


static double
distance_1D(double a1, double a2, double b1, double b2)
{
	/* interval (a) is entirely on the left of (b) */
	if ((a1 <= b1) && (a2 <= b1) && (a1 <= b2) && (a2 <= b2))
		return (Min(b1, b2) - Max(a1, a2));

	/* interval (a) is entirely on the right of (b) */
	if ((a1 > b1) && (a2 > b1) && (a1 > b2) && (a2 > b2))
		return (Min(a1, a2) - Max(b1, b2));

	/* the rest are all sorts of intersections */
	return (0.0);
}


/* Test if a box is also a point */
Datum
cube_is_point(PG_FUNCTION_ARGS)
{
	NDBOX		*cube = PG_GETARG_NDBOX(0);
	if (IS_POINT(cube))
		PG_RETURN_BOOL(TRUE);
	else
		PG_RETURN_BOOL(FALSE);
}


/* Return dimensions in use in the data structure */
Datum
cube_dim(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	int			dim = DIM(c);
	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_INT32(dim);
}


/* Return a specific normalized LL coordinate */
Datum
cube_ll_coord(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	int			n = PG_GETARG_INT16(1);
	double		result;

	if (DIM(c) >= n && n > 0)
		result = Min(LL_COORD(c, n-1), UR_COORD(c, n-1));
	else
		result = 0;

	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_FLOAT8(result);
}


/* Return a specific normalized UR coordinate */
Datum
cube_ur_coord(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	int			n = PG_GETARG_INT16(1);
	double		result;

	if (DIM(c) >= n && n > 0)
		result = Max(LL_COORD(c, n-1), UR_COORD(c, n-1));
	else
		result = 0;

	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_FLOAT8(result);
}


/* Increase or decrease box size by a radius in at least n dimensions. */
Datum
cube_enlarge(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	double		r = PG_GETARG_FLOAT8(1);
	int32		n = PG_GETARG_INT32(2);
	NDBOX	   *result;
	int			i, j,
				dim = 0,
				shrunk_coordinates = 0;

	if (n > CUBE_MAX_DIM)
		n = CUBE_MAX_DIM;
	if (r > 0 && n > 0)
		dim = n;
	if (DIM(a) > dim)
		dim = DIM(a);

	result = init_cube( dim, 0, TYPE(a) );
	
	for (i = 0, j = dim; i < DIM(a); i++, j++)
	{
		if (LL_COORD(a,i) >= UR_COORD(a,i))
		{
			set_coord(result, i, UR_COORD(a,i) - r);
			set_coord(result, j, LL_COORD(a,i) + r);
		}
		else
		{
			set_coord(result, i, LL_COORD(a,i) - r);
			set_coord(result, j, UR_COORD(a,i) + r);
		}
		if (LL_COORD(result, i) > LL_COORD(result, j))
		{
			set_coord(result, i,
				(LL_COORD(result, i) + LL_COORD(result, j)) / 2);
			set_coord(result, j, LL_COORD(result, i));
			shrunk_coordinates++;
		}
		else if (LL_COORD(result, i)  == LL_COORD(result, j))
			shrunk_coordinates++;
	}
	/* dim > a->dim only if r > 0 */
	for (; i < dim; i++, j++)
	{
		set_coord(result, i, -r);
		set_coord(result, j, r);
	}

	/* Point can arise in two cases:
	   1) When argument is point and r == 0
	   2) When all coordinates was set to their averages */
	if ( (IS_POINT(a) && r == 0) || (shrunk_coordinates == dim) )
		cube_to_point(result);

	PG_FREE_IF_COPY(a, 0);
	PG_RETURN_NDBOX(result);
}


/* Add a dimension to an existing cube with the same values for the new
   coordinate */
Datum
cube_c_f8(PG_FUNCTION_ARGS)
{
	NDBOX		*cube = PG_GETARG_NDBOX(0);
	double		x = PG_GETARG_FLOAT8(1);
	NDBOX		*result;
	int			i;

	if (IS_POINT(cube))
	{
		result = init_cube( DIM(cube) + 1, 1, TYPE(cube) );
		for (i = 0; i < DIM(cube); i++)
			set_coord(result, i, LL_COORD(cube, i));
		set_coord(result, DIM(result) - 1, x);
	}
	else
	{
		result = init_cube( DIM(cube) + 1, 0, TYPE(cube) );
		for (i = 0; i < DIM(cube); i++)
		{
			set_coord(result, i, LL_COORD(cube, i));
			set_coord(result, DIM(result) + i, UR_COORD(cube, i));
		}
		set_coord(result, DIM(result) - 1, x);
		set_coord(result, 2*DIM(result) - 1, x);
	}

	PG_FREE_IF_COPY(cube, 0);
	PG_RETURN_NDBOX(result);
}


/* Add a dimension to an existing cube */
Datum
cube_c_f8_f8(PG_FUNCTION_ARGS)
{
	NDBOX	   *cube = PG_GETARG_NDBOX(0);
	double		x1 = PG_GETARG_FLOAT8(1);
	double		x2 = PG_GETARG_FLOAT8(2);
	NDBOX	   *result;
	int			i;

	if (IS_POINT(cube) && (x1 == x2)){
		result = init_cube( DIM(cube) + 1, 1, TYPE(cube) );
		for (i = 0; i < DIM(cube); i++)
			set_coord(result, i, LL_COORD(cube, i));
		set_coord(result, DIM(result) - 1, x1);
	}
	else
	{
		result = init_cube( DIM(cube) + 1, 0, TYPE(cube) );
		for (i = 0; i < DIM(cube); i++)
		{
			set_coord(result, i, LL_COORD(cube, i));
			set_coord(result, DIM(result) + i, UR_COORD(cube, i));
		}
		set_coord(result, DIM(result) - 1, x1);
		set_coord(result, 2*DIM(result) - 1, x2);
	}

	PG_FREE_IF_COPY(cube, 0);
	PG_RETURN_NDBOX(result);
}

