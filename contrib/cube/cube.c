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
PG_FUNCTION_INFO_V1(cube_f8);
PG_FUNCTION_INFO_V1(cube_f8_f8);
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
Datum		cube_f8(PG_FUNCTION_ARGS);
Datum		cube_f8_f8(PG_FUNCTION_ARGS);
Datum		cube_c_f8(PG_FUNCTION_ARGS);
Datum		cube_c_f8_f8(PG_FUNCTION_ARGS);
Datum		cube_dim(PG_FUNCTION_ARGS);
Datum		cube_ll_coord(PG_FUNCTION_ARGS);
Datum		cube_ur_coord(PG_FUNCTION_ARGS);
Datum		cube_subset(PG_FUNCTION_ARGS);

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
static int				interval_cmp_lower(const void *i1, const void *i2);
static int				interval_cmp_upper(const void *i1, const void *i2);
static int				common_entry_cmp(const void *i1, const void *i2);
static void				guttman_split(GistEntryVector *entryvec, GIST_SPLITVEC *v);
static void				korotkov_split(GistEntryVector *entryvec, GIST_SPLITVEC *v);
static void				fallback_split(GistEntryVector *entryvec, GIST_SPLITVEC *v);

static inline int32		cube_cmp_v0(NDBOX *a, NDBOX *b);
static inline bool		cube_contains_v0(NDBOX *a, NDBOX *b);
static inline bool		cube_overlap_v0(NDBOX *a, NDBOX *b);
static inline NDBOX	   *cube_union_v0(NDBOX *a, NDBOX *b);
static inline NDBOX	   *_cube_inter(NDBOX *a, NDBOX *b);
static inline void		rt_cube_size(NDBOX *a, double *sz);
static inline double	distance_1D(double a1, double a2, double b1, double b2);
static inline double	cube_penalty(NDBOX *origentry, NDBOX *newentry);

static inline float 	non_negative(float val);
static inline NDBOX	   *adjust_box(NDBOX *b, NDBOX *addon);
static inline void		cube_consider_split(ConsiderSplitContext *context, int dimNum,
						  double rightLower, int minLeftCount,
						  double leftUpper, int maxLeftCount);

NDBOX	   *g_cube_binary_union(NDBOX *r1, NDBOX *r2, int *sizep);
bool		g_cube_leaf_consistent(NDBOX *key, NDBOX *query, StrategyNumber strategy);
bool		g_cube_internal_consistent(NDBOX *key, NDBOX *query, StrategyNumber strategy);

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

	cube_scanner_init(str);

	if (cube_yyparse(&result) != 0)
		cube_yyerror(&result, "bogus input");

	cube_scanner_finish();

	PG_RETURN_NDBOX(result);
}


/*
** Allows the construction of a cube from 2 float[]'s
*/
Datum
cube_a_f8_f8(PG_FUNCTION_ARGS)
{
	ArrayType  *ur = PG_GETARG_ARRAYTYPE_P(0);
	ArrayType  *ll = PG_GETARG_ARRAYTYPE_P(1);
	NDBOX	   *result;
	int			i;
	int			dim;
	int			size;
	double	   *dur,
			   *dll;

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

	size = offsetof(NDBOX, x[0]) +sizeof(double) * 2 * dim;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = dim;

	for (i = 0; i < dim; i++)
	{
		result->x[i] = dur[i];
		result->x[i + dim] = dll[i];
	}

	PG_RETURN_NDBOX(result);
}

/*
** Allows the construction of a zero-volume cube from a float[]
*/
Datum
cube_a_f8(PG_FUNCTION_ARGS)
{
	ArrayType  *ur = PG_GETARG_ARRAYTYPE_P(0);
	NDBOX	   *result;
	int			i;
	int			dim;
	int			size;
	double	   *dur;

	if (array_contains_nulls(ur))
		ereport(ERROR,
				(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
				 errmsg("cannot work with arrays containing NULLs")));

	dim = ARRNELEMS(ur);

	dur = ARRPTR(ur);

	size = offsetof(NDBOX, x[0]) +sizeof(double) * 2 * dim;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = dim;

	for (i = 0; i < dim; i++)
	{
		result->x[i] = dur[i];
		result->x[i + dim] = dur[i];
	}

	PG_RETURN_NDBOX(result);
}

Datum
cube_subset(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	ArrayType  *idx = PG_GETARG_ARRAYTYPE_P(1);
	NDBOX	   *result;
	int			size,
				dim,
				i;
	int		   *dx;

	if (array_contains_nulls(idx))
		ereport(ERROR,
				(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
				 errmsg("cannot work with arrays containing NULLs")));

	dx = (int32 *) ARR_DATA_PTR(idx);

	dim = ARRNELEMS(idx);
	size = offsetof(NDBOX, x[0]) +sizeof(double) * 2 * dim;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = dim;

	for (i = 0; i < dim; i++)
	{
		if ((dx[i] <= 0) || (dx[i] > c->dim))
		{
			pfree(result);
			ereport(ERROR,
					(errcode(ERRCODE_ARRAY_ELEMENT_ERROR),
					 errmsg("Index out of bounds")));
		}
		result->x[i] = c->x[dx[i] - 1];
		result->x[i + dim] = c->x[dx[i] + c->dim - 1];
	}

	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_NDBOX(result);
}

Datum
cube_out(PG_FUNCTION_ARGS)
{
	NDBOX	   *cube = PG_GETARG_NDBOX(0);
	StringInfoData buf;
	int			dim = cube->dim;
	bool		equal = true;
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
		appendStringInfo(&buf, "%.*g", ndig, cube->x[i]);
		if (cube->x[i] != cube->x[i + dim])
			equal = false;
	}
	appendStringInfoChar(&buf, ')');

	if (!equal)
	{
		appendStringInfo(&buf, ",(");
		for (i = 0; i < dim; i++)
		{
			if (i > 0)
				appendStringInfo(&buf, ", ");
			appendStringInfo(&buf, "%.*g", ndig, cube->x[i + dim]);
		}
		appendStringInfoChar(&buf, ')');
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

	/*
	 * fprintf(stderr, "union\n");
	 */
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
	PG_RETURN_DATUM(PG_GETARG_DATUM(0));
}

Datum
g_cube_decompress(PG_FUNCTION_ARGS)
{
	GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
	NDBOX	   *key = DatumGetNDBOX(PG_DETOAST_DATUM(entry->key));

	if (key != DatumGetNDBOX(entry->key))
	{
		GISTENTRY  *retval = (GISTENTRY *) palloc(sizeof(GISTENTRY));

		gistentryinit(*retval, PointerGetDatum(key),
					  entry->rel, entry->page,
					  entry->offset, FALSE);
		PG_RETURN_POINTER(retval);
	}
	PG_RETURN_POINTER(entry);
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

	*result = cube_penalty(
		DatumGetNDBOX(origentry->key),
		DatumGetNDBOX(newentry->key)
	);

	/*
	 * fprintf(stderr, "penalty\n"); fprintf(stderr, "\t%g\n", *result);
	 */
	PG_RETURN_FLOAT8(*result);
}

/*
** The GiST PickSplit method for boxes
** We use Guttman's poly time split algorithm
*/
Datum
g_cube_picksplit(PG_FUNCTION_ARGS)
{
	GistEntryVector	*entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
	GIST_SPLITVEC	*v = (GIST_SPLITVEC *) PG_GETARG_POINTER(1);

	if (DatumGetNDBOX(entryvec->vector[FirstOffsetNumber].key)->dim >= SPLIT_THRESHOLD)
		guttman_split(entryvec, v);
	else
		korotkov_split(entryvec, v);

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

	/*
	 * fprintf(stderr, "same: %s\n", (*result ? "TRUE" : "FALSE" ));
	 */
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

	/*
	 * fprintf(stderr, "leaf_consistent, %d\n", strategy);
	 */
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

	/*
	 * fprintf(stderr, "internal_consistent, %d\n", strategy);
	 */
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
static inline NDBOX *
cube_union_v0(NDBOX *a, NDBOX *b)
{
	int			i;
	NDBOX	   *result;

	if (a->dim >= b->dim)
	{
		result = palloc0(VARSIZE(a));
		SET_VARSIZE(result, VARSIZE(a));
		result->dim = a->dim;
	}
	else
	{
		result = palloc0(VARSIZE(b));
		SET_VARSIZE(result, VARSIZE(b));
		result->dim = b->dim;
	}

	/* swap the box pointers if needed */
	if (a->dim < b->dim)
	{
		NDBOX	   *tmp = b;

		b = a;
		a = tmp;
	}

	/*
	 * use the potentially smaller of the two boxes (b) to fill in the result,
	 * padding absent dimensions with zeroes
	 */
	for (i = 0; i < b->dim; i++)
	{
		result->x[i] = Min(b->x[i], b->x[i + b->dim]);
		result->x[i + a->dim] = Max(b->x[i], b->x[i + b->dim]);
	}
	for (i = b->dim; i < a->dim; i++)
	{
		result->x[i] = 0;
		result->x[i + a->dim] = 0;
	}

	/* compute the union */
	for (i = 0; i < a->dim; i++)
	{
		result->x[i] =
			Min(Min(a->x[i], a->x[i + a->dim]), result->x[i]);
		result->x[i + a->dim] = Max(Max(a->x[i],
								   a->x[i + a->dim]), result->x[i + a->dim]);
	}

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
static inline NDBOX*
_cube_inter(NDBOX *a, NDBOX *b)
{
	NDBOX	   *result, *cube1, *cube2;
	int			i;
	double		edge1, edge2;

	cube1 = (a->dim > b->dim) ? a : b;
	cube2 = (a->dim <= b->dim) ? a : b;

	result = palloc0(VARSIZE(cube1));
	SET_VARSIZE(result, VARSIZE(cube1));
	result->dim = cube1->dim;

	/* compute the intersection */
	for (i = 0; i < cube2->dim; i++)
	{
		edge1 = Max(
			Min(cube1->x[i], cube1->x[i + cube1->dim]),
			Min(cube2->x[i], cube2->x[i + cube2->dim])
		);
		edge2 = Min(
			Max(cube1->x[i], cube1->x[i + cube1->dim]),
			Max(cube2->x[i], cube2->x[i + cube2->dim])
		);


		if (edge1 <= edge2)
		{
			result->x[i] = edge1;
			result->x[i + result->dim] = edge2;
		}
		else
		{
			result->x[i] = 0;
			result->x[i + result->dim] = 0;
		}
	}

	for (; i < cube2->dim; i++)
	{
		edge1 = Max(
			Min(cube1->x[i], cube1->x[i + cube1->dim]),
			0
		);
		edge2 = Min(
			Max(cube1->x[i], cube1->x[i + cube1->dim]),
			0
		);

		if (edge1 <= edge2)
		{
			result->x[i] = edge1;
			result->x[i + result->dim] = edge2;
		}
		else
		{
			result->x[i] = 0;
			result->x[i + result->dim] = 0;
		}
	}

	return result;
}

Datum
cube_inter(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	NDBOX	   *b = PG_GETARG_NDBOX(1);
	NDBOX	   *result;

	result = _cube_inter(a, b);

	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_NDBOX(result);
}

/* cube_size */
Datum
cube_size(PG_FUNCTION_ARGS)
{
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	double		result;
	int			i,
				j;

	result = 1.0;
	for (i = 0, j = a->dim; i < a->dim; i++, j++)
		result = result * Abs((a->x[j] - a->x[i]));

	PG_FREE_IF_COPY(a, 0);
	PG_RETURN_FLOAT8(result);
}

static inline void
rt_cube_size(NDBOX *a, double *size)
{
	int			i,
				j;

	if (a == (NDBOX *) NULL)
		*size = 0.0;
	else
	{
		*size = 1.0;
		for (i = 0, j = a->dim; i < a->dim; i++, j++)
			*size = (*size) * Abs((a->x[j] - a->x[i]));
	}
	return;
}

/* make up a metric in which one box will be 'lower' than the other
   -- this can be useful for sorting and to determine uniqueness */
static inline int32
cube_cmp_v0(NDBOX *a, NDBOX *b)
{
	int			i;
	int			dim;

	dim = Min(a->dim, b->dim);

	/* compare the common dimensions */
	for (i = 0; i < dim; i++)
	{
		if (Min(a->x[i], a->x[a->dim + i]) >
			Min(b->x[i], b->x[b->dim + i]))
			return 1;
		if (Min(a->x[i], a->x[a->dim + i]) <
			Min(b->x[i], b->x[b->dim + i]))
			return -1;
	}
	for (i = 0; i < dim; i++)
	{
		if (Max(a->x[i], a->x[a->dim + i]) >
			Max(b->x[i], b->x[b->dim + i]))
			return 1;
		if (Max(a->x[i], a->x[a->dim + i]) <
			Max(b->x[i], b->x[b->dim + i]))
			return -1;
	}

	/* compare extra dimensions to zero */
	if (a->dim > b->dim)
	{
		for (i = dim; i < a->dim; i++)
		{
			if (Min(a->x[i], a->x[a->dim + i]) > 0)
				return 1;
			if (Min(a->x[i], a->x[a->dim + i]) < 0)
				return -1;
		}
		for (i = dim; i < a->dim; i++)
		{
			if (Max(a->x[i], a->x[a->dim + i]) > 0)
				return 1;
			if (Max(a->x[i], a->x[a->dim + i]) < 0)
				return -1;
		}

		/*
		 * if all common dimensions are equal, the cube with more dimensions
		 * wins
		 */
		return 1;
	}
	if (a->dim < b->dim)
	{
		for (i = dim; i < b->dim; i++)
		{
			if (Min(b->x[i], b->x[b->dim + i]) > 0)
				return -1;
			if (Min(b->x[i], b->x[b->dim + i]) < 0)
				return 1;
		}
		for (i = dim; i < b->dim; i++)
		{
			if (Max(b->x[i], b->x[b->dim + i]) > 0)
				return -1;
			if (Max(b->x[i], b->x[b->dim + i]) < 0)
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

static inline NDBOX*
adjust_box(NDBOX *b, NDBOX *addon)
{
	int			i, size;
	NDBOX	   *cube;

	if (b->dim >= addon->dim)
		cube = b;
	else
	{
		size = offsetof(NDBOX, x[0]) + 2*sizeof(double)*addon->dim;
		cube = (NDBOX *) palloc0(size);
		SET_VARSIZE(cube, size);
		cube->dim = addon->dim;

		for (i=0; i<b->dim; i++)
			cube->x[cube->dim + i] = b->x[b->dim + i];

		for (; i<cube->dim; i++)
			cube->x[i] = cube->x[cube->dim + i] = 0;
	}

	for(i=0; i < addon->dim; i++)
	{
		if (cube->x[i + cube->dim] < addon->x[i + addon->dim])
			cube->x[i + cube->dim] = addon->x[i + addon->dim];
		if (cube->x[i] > addon->x[i])
			cube->x[i] = addon->x[i];
	}

	return cube;
}

static int
interval_cmp_lower(const void *i1, const void *i2)
{
	double		lower1 = ((const SplitInterval *) i1)->lower,
				lower2 = ((const SplitInterval *) i2)->lower;

	if (lower1 < lower2)
		return -1;
	else if (lower1 > lower2)
		return 1;
	else
		return 0;
}

static int
interval_cmp_upper(const void *i1, const void *i2)
{
	double		upper1 = ((const SplitInterval *) i1)->upper,
				upper2 = ((const SplitInterval *) i2)->upper;

	if (upper1 < upper2)
		return -1;
	else if (upper1 > upper2)
		return 1;
	else
		return 0;
}

static inline float
non_negative(float val)
{
	if (val >= 0.0f)
		return val;
	else
		return 0.0f;
}

static int
common_entry_cmp(const void *i1, const void *i2)
{
	double		delta1 = ((const CommonEntry *) i1)->delta,
				delta2 = ((const CommonEntry *) i2)->delta;

	if (delta1 < delta2)
		return -1;
	else if (delta1 > delta2)
		return 1;
	else
		return 0;
}

static inline double
cube_penalty(NDBOX *origentry, NDBOX *newentry)
{
	double		tmp1,
				tmp2;

	rt_cube_size(cube_union_v0(origentry, newentry), &tmp1);
	rt_cube_size(origentry, &tmp2);
	return (tmp1 - tmp2);
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
static inline bool
cube_contains_v0(NDBOX *a, NDBOX *b)
{
	int			i;

	if ((a == NULL) || (b == NULL))
		return (FALSE);

	if (a->dim < b->dim)
	{
		/*
		 * the further comparisons will make sense if the excess dimensions of
		 * (b) were zeroes Since both UL and UR coordinates must be zero, we
		 * can check them all without worrying about which is which.
		 */
		for (i = a->dim; i < b->dim; i++)
		{
			if (b->x[i] != 0)
				return (FALSE);
			if (b->x[i + b->dim] != 0)
				return (FALSE);
		}
	}

	/* Can't care less about the excess dimensions of (a), if any */
	for (i = 0; i < Min(a->dim, b->dim); i++)
	{
		if (Min(a->x[i], a->x[a->dim + i]) >
			Min(b->x[i], b->x[b->dim + i]))
			return (FALSE);
		if (Max(a->x[i], a->x[a->dim + i]) <
			Max(b->x[i], b->x[b->dim + i]))
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
static inline bool
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
	if (a->dim < b->dim)
	{
		NDBOX	   *tmp = b;

		b = a;
		a = tmp;
	}

	/* compare within the dimensions of (b) */
	for (i = 0; i < b->dim; i++)
	{
		if (Min(a->x[i], a->x[a->dim + i]) >
			Max(b->x[i], b->x[b->dim + i]))
			return (FALSE);
		if (Max(a->x[i], a->x[a->dim + i]) <
			Min(b->x[i], b->x[b->dim + i]))
			return (FALSE);
	}

	/* compare to zero those dimensions in (a) absent in (b) */
	for (i = b->dim; i < a->dim; i++)
	{
		if (Min(a->x[i], a->x[a->dim + i]) > 0)
			return (FALSE);
		if (Max(a->x[i], a->x[a->dim + i]) < 0)
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
	if (a->dim < b->dim)
	{
		NDBOX	   *tmp = b;

		b = a;
		a = tmp;
		swapped = true;
	}

	distance = 0.0;
	/* compute within the dimensions of (b) */
	for (i = 0; i < b->dim; i++)
	{
		d = distance_1D(a->x[i], a->x[i + a->dim], b->x[i], b->x[i + b->dim]);
		distance += d * d;
	}

	/* compute distance to zero for those dimensions in (a) absent in (b) */
	for (i = b->dim; i < a->dim; i++)
	{
		d = distance_1D(a->x[i], a->x[i + a->dim], 0.0, 0.0);
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

static inline double
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
	NDBOX	   *a = PG_GETARG_NDBOX(0);
	int			i,
				j;

	for (i = 0, j = a->dim; i < a->dim; i++, j++)
	{
		if (a->x[i] != a->x[j])
			PG_RETURN_BOOL(FALSE);
	}

	PG_FREE_IF_COPY(a, 0);
	PG_RETURN_BOOL(TRUE);
}

/* Return dimensions in use in the data structure */
Datum
cube_dim(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	int			dim = c->dim;

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

	if (c->dim >= n && n > 0)
		result = Min(c->x[n - 1], c->x[c->dim + n - 1]);
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

	if (c->dim >= n && n > 0)
		result = Max(c->x[n - 1], c->x[c->dim + n - 1]);
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
	int			dim = 0;
	int			size;
	int			i,
				j,
				k;

	if (n > CUBE_MAX_DIM)
		n = CUBE_MAX_DIM;
	if (r > 0 && n > 0)
		dim = n;
	if (a->dim > dim)
		dim = a->dim;
	size = offsetof(NDBOX, x[0]) +sizeof(double) * dim * 2;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = dim;
	for (i = 0, j = dim, k = a->dim; i < a->dim; i++, j++, k++)
	{
		if (a->x[i] >= a->x[k])
		{
			result->x[i] = a->x[k] - r;
			result->x[j] = a->x[i] + r;
		}
		else
		{
			result->x[i] = a->x[i] - r;
			result->x[j] = a->x[k] + r;
		}
		if (result->x[i] > result->x[j])
		{
			result->x[i] = (result->x[i] + result->x[j]) / 2;
			result->x[j] = result->x[i];
		}
	}
	/* dim > a->dim only if r > 0 */
	for (; i < dim; i++, j++)
	{
		result->x[i] = -r;
		result->x[j] = r;
	}

	PG_FREE_IF_COPY(a, 0);
	PG_RETURN_NDBOX(result);
}

/* Create a one dimensional box with identical upper and lower coordinates */
Datum
cube_f8(PG_FUNCTION_ARGS)
{
	double		x = PG_GETARG_FLOAT8(0);
	NDBOX	   *result;
	int			size;

	size = offsetof(NDBOX, x[0]) +sizeof(double) * 2;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = 1;
	result->x[0] = result->x[1] = x;

	PG_RETURN_NDBOX(result);
}

/* Create a one dimensional box */
Datum
cube_f8_f8(PG_FUNCTION_ARGS)
{
	double		x0 = PG_GETARG_FLOAT8(0);
	double		x1 = PG_GETARG_FLOAT8(1);
	NDBOX	   *result;
	int			size;

	size = offsetof(NDBOX, x[0]) +sizeof(double) * 2;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = 1;
	result->x[0] = x0;
	result->x[1] = x1;

	PG_RETURN_NDBOX(result);
}

/* Add a dimension to an existing cube with the same values for the new
   coordinate */
Datum
cube_c_f8(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	double		x = PG_GETARG_FLOAT8(1);
	NDBOX	   *result;
	int			size;
	int			i;

	size = offsetof(NDBOX, x[0]) +sizeof(double) * (c->dim + 1) *2;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = c->dim + 1;
	for (i = 0; i < c->dim; i++)
	{
		result->x[i] = c->x[i];
		result->x[result->dim + i] = c->x[c->dim + i];
	}
	result->x[result->dim - 1] = x;
	result->x[2 * result->dim - 1] = x;

	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_NDBOX(result);
}

/* Add a dimension to an existing cube */
Datum
cube_c_f8_f8(PG_FUNCTION_ARGS)
{
	NDBOX	   *c = PG_GETARG_NDBOX(0);
	double		x1 = PG_GETARG_FLOAT8(1);
	double		x2 = PG_GETARG_FLOAT8(2);
	NDBOX	   *result;
	int			size;
	int			i;

	size = offsetof(NDBOX, x[0]) +sizeof(double) * (c->dim + 1) *2;
	result = (NDBOX *) palloc0(size);
	SET_VARSIZE(result, size);
	result->dim = c->dim + 1;
	for (i = 0; i < c->dim; i++)
	{
		result->x[i] = c->x[i];
		result->x[result->dim + i] = c->x[c->dim + i];
	}
	result->x[result->dim - 1] = x1;
	result->x[2 * result->dim - 1] = x2;

	PG_FREE_IF_COPY(c, 0);
	PG_RETURN_NDBOX(result);
}


/*
 * Split routines.
 */

static void
guttman_split(GistEntryVector *entryvec, GIST_SPLITVEC *v)
{
	OffsetNumber	i,
					j,
					nelems,
					seed_1,
					seed_2,
					selected_offset;
	NDBOX		   *cube_alpha,
				   *cube_beta,
				   *union_cube,
				   *inter_cube,
				   *cube_left,
				   *cube_right,
				   *current_cube,
				   *left_union,
				   *right_union,
				   *selected_cube,
				   *cube;
	double			size_waste,
					waste,
					volume_left,
					volume_right,
					current_increase_l,
					current_increase_r,
					current_diff_increase,
					diff_increase,
					increase_l,
					increase_r,
					sz1,
					sz2;
	int				nbytes,
					min_count,
					remaining_count;
	bool			firsttime,
				   *inserted_cubes;

	nelems = entryvec->n - 1;
	nbytes = nelems*sizeof(OffsetNumber) + 1;
	v->spl_left = (OffsetNumber *) palloc(nbytes);
	v->spl_right = (OffsetNumber *) palloc(nbytes);
	inserted_cubes = palloc0(nelems + 1);
	v->spl_nleft = v->spl_nright = 0;

	firsttime = true;
	waste = 0.0;

	/* 
	 * pick two seed elements
	 */
	for (i = FirstOffsetNumber; i <= nelems; i = OffsetNumberNext(i))
	{
		cube_alpha = DatumGetNDBOX(entryvec->vector[i].key);
		for (j = OffsetNumberNext(i); j <= nelems; j = OffsetNumberNext(j))
		{
			cube_beta = DatumGetNDBOX(entryvec->vector[j].key);

			/* compute the wasted space by unioning these guys */
			/* size_waste = size_union - size_inter; */
			union_cube = cube_union_v0(cube_alpha, cube_beta);
			inter_cube = _cube_inter(cube_alpha, cube_beta);
			rt_cube_size(union_cube, &sz1);
			rt_cube_size(inter_cube, &sz2);
			size_waste = sz1 - sz2;

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

	/* selected elements became first elements in groups */
	v->spl_left[v->spl_nleft++] = seed_1;
	v->spl_right[v->spl_nright++] = seed_2;
	inserted_cubes[seed_1] = inserted_cubes[seed_2] = true;
	
	/* bounding boxes for groups */
	cube = DatumGetNDBOX(entryvec->vector[seed_1].key);
	cube_left = palloc(VARSIZE(cube));
	memcpy(cube_left, cube, VARSIZE(cube));

	cube = DatumGetNDBOX(entryvec->vector[seed_2].key);
	cube_right = palloc(VARSIZE(cube));
	memcpy(cube_right, cube, VARSIZE(cube));

	rt_cube_size(cube_left, &volume_left);
	rt_cube_size(cube_right, &volume_right);

	/*
	 * insert cubes
	 */
	remaining_count = nelems-2;
	min_count = ceil(LIMIT_RATIO*(nelems - 1));

	while (remaining_count > 0)
	{
		diff_increase = 0.0;

		/* check if there ehough cubes to satisfy mininum elements restriction */
		if ((v->spl_nleft + remaining_count) <= min_count || (v->spl_nright + remaining_count) <= min_count )
			break;

		/* find cube with biggest diff_increase to insert it */
		for (i = FirstOffsetNumber; i <= nelems; i = OffsetNumberNext(i))
		{
			current_cube = DatumGetNDBOX(entryvec->vector[i].key);
		
			/* skip already inserted cubes */
			if (inserted_cubes[i])
				continue;

			left_union = cube_union_v0(cube_left, current_cube);
			rt_cube_size(union_cube, &sz1);
			current_increase_l = sz1 - volume_left;

			right_union = cube_union_v0(cube_right, current_cube);
			rt_cube_size(union_cube, &sz2);
			current_increase_r = sz2 - volume_right;

			current_diff_increase = Abs(current_increase_l-current_increase_r);

			if (current_diff_increase >= diff_increase)
			{
				diff_increase = current_diff_increase;
				increase_l = current_increase_l;
				increase_r = current_increase_r;
				selected_offset = i;
			}
		}

		/* now insert selected cube on leaf with minimal increase */
		selected_cube = DatumGetNDBOX(entryvec->vector[selected_offset].key);
		if (increase_l < increase_r)
		{
			v->spl_left[v->spl_nleft++] = selected_offset;
			cube_left = adjust_box(cube_left, selected_cube);
			rt_cube_size(cube_left, &volume_left);
		}
		else
		{
			v->spl_right[v->spl_nright++] = selected_offset;
			cube_right = adjust_box(cube_right, selected_cube);
			rt_cube_size(cube_right, &volume_right);
		}

		inserted_cubes[selected_offset] = true;	/* mark this cube as inserted */
		remaining_count--;
	}

	/* insert all remaining entries (if any) to one branch */
	if (remaining_count > 0)
	{
		if (v->spl_nleft < v->spl_nright)
		{
			for (i = FirstOffsetNumber; i <= nelems; i = OffsetNumberNext(i))
			{
				if (!inserted_cubes[i]){
					v->spl_left[v->spl_nleft++] = i;
					cube_left = adjust_box(cube_left, DatumGetNDBOX(entryvec->vector[i].key));
				}
			}
		}
		else
		{
			for (i = FirstOffsetNumber; i <= nelems; i = OffsetNumberNext(i))
			{
				if (!inserted_cubes[i]){
					v->spl_right[v->spl_nright++] = i;
					cube_right = adjust_box(cube_right, DatumGetNDBOX(entryvec->vector[i].key));
				}
			}
		}
	}

	v->spl_ldatum = PointerGetDatum(cube_left);
	v->spl_rdatum = PointerGetDatum(cube_right);
	return;
}

static void
korotkov_split(GistEntryVector *entryvec, GIST_SPLITVEC *v)
{
	ConsiderSplitContext context;
	OffsetNumber	i,
					maxoff;
	NDBOX		   *box,
				   *leftBox,
				   *rightBox;
	int				dim,
					commonEntriesCount;
	SplitInterval  *intervalsLower,
				   *intervalsUpper;
	CommonEntry	   *commonEntries;
	int				nentries;

	memset(&context, 0, sizeof(ConsiderSplitContext));

	maxoff = entryvec->n - 1;
	nentries = context.entriesCount = maxoff - FirstOffsetNumber + 1;

	/* Allocate arrays for intervals along axes */
	intervalsLower = (SplitInterval *) palloc(nentries * sizeof(SplitInterval));
	intervalsUpper = (SplitInterval *) palloc(nentries * sizeof(SplitInterval));

	/*
	 * Calculate the overall minimum bounding box over all the entries.
	 */
	for (i = FirstOffsetNumber; i <= maxoff; i = OffsetNumberNext(i))
	{
		box = DatumGetNDBOX(entryvec->vector[i].key);
		if (i == FirstOffsetNumber)
		{
			context.boundingBox = palloc(VARSIZE(box));
			memcpy(context.boundingBox, box, VARSIZE(box));
		}
		else
			context.boundingBox = adjust_box(context.boundingBox, box);
	}

	/*
	 * Iterate over axes for optimal split searching.
	 */
	context.first = true;		/* nothing selected yet */
	for (dim = 0; dim < box->dim; dim++)
	{
		double		leftUpper,
					rightLower;
		int			i1,
					i2;

		/* Project each entry as an interval on the selected axis. */
		for (i = FirstOffsetNumber; i <= maxoff; i = OffsetNumberNext(i))
		{
			box = DatumGetNDBOX(entryvec->vector[i].key);
			intervalsLower[i - FirstOffsetNumber].lower = box->x[dim];
			intervalsLower[i - FirstOffsetNumber].upper = box->x[dim+box->dim];
		}

		/*
		 * Make two arrays of intervals: one sorted by lower bound and another
		 * sorted by upper bound.
		 */
		memcpy(intervalsUpper, intervalsLower,
			   sizeof(SplitInterval) * nentries);
		qsort(intervalsLower, nentries, sizeof(SplitInterval),
			  interval_cmp_lower);
		qsort(intervalsUpper, nentries, sizeof(SplitInterval),
			  interval_cmp_upper);

		/*----
		 * The goal is to form a left and right interval, so that every entry
		 * interval is contained by either left or right interval (or both).
		 *
		 * For example, with the intervals (0,1), (1,3), (2,3), (2,4):
		 *
		 * 0 1 2 3 4
		 * +-+
		 *	 +---+
		 *	   +-+
		 *	   +---+
		 *
		 * The left and right intervals are of the form (0,a) and (b,4).
		 * We first consider splits where b is the lower bound of an entry.
		 * We iterate through all entries, and for each b, calculate the
		 * smallest possible a. Then we consider splits where a is the
		 * uppper bound of an entry, and for each a, calculate the greatest
		 * possible b.
		 *
		 * In the above example, the first loop would consider splits:
		 * b=0: (0,1)-(0,4)
		 * b=1: (0,1)-(1,4)
		 * b=2: (0,3)-(2,4)
		 *
		 * And the second loop:
		 * a=1: (0,1)-(1,4)
		 * a=3: (0,3)-(2,4)
		 * a=4: (0,4)-(2,4)
		 */

		/*
		 * Iterate over lower bound of right group, finding smallest possible
		 * upper bound of left group.
		 */
		i1 = 0;
		i2 = 0;
		rightLower = intervalsLower[i1].lower;
		leftUpper = intervalsUpper[i2].lower;
		while (true)
		{
			/*
			 * Find next lower bound of right group.
			 */
			while (i1 < nentries && rightLower == intervalsLower[i1].lower)
			{
				leftUpper = Max(leftUpper, intervalsLower[i1].upper);
				i1++;
			}
			if (i1 >= nentries)
				break;
			rightLower = intervalsLower[i1].lower;

			/*
			 * Find count of intervals which anyway should be placed to the
			 * left group.
			 */
			while (i2 < nentries && intervalsUpper[i2].upper <= leftUpper)
				i2++;

			/*
			 * Consider found split.
			 */
			cube_consider_split(&context, dim, rightLower, i1, leftUpper, i2);
		}

		/*
		 * Iterate over upper bound of left group finding greates possible
		 * lower bound of right group.
		 */
		i1 = nentries - 1;
		i2 = nentries - 1;
		rightLower = intervalsLower[i1].upper;
		leftUpper = intervalsUpper[i2].upper;
		while (true)
		{
			/*
			 * Find next upper bound of left group.
			 */
			while (i2 >= 0 && leftUpper == intervalsUpper[i2].upper)
			{
				rightLower = Min(rightLower, intervalsUpper[i2].lower);
				i2--;
			}
			if (i2 < 0)
				break;
			leftUpper = intervalsUpper[i2].upper;

			/*
			 * Find count of intervals which anyway should be placed to the
			 * right group.
			 */
			while (i1 >= 0 && intervalsLower[i1].lower >= rightLower)
				i1--;

			/*
			 * Consider found split.
			 */
			cube_consider_split(&context, dim,
								 rightLower, i1 + 1, leftUpper, i2 + 1);
		}
	}

	/*
	 * If we failed to find any acceptable splits, use trivial split.
	 */
	if (context.first)
	{
		fallback_split(entryvec, v);
		return;
	}

	/*
	 * Ok, we have now selected the split across one axis.
	 *
	 * While considering the splits, we already determined that there will be
	 * enough entries in both groups to reach the desired ratio, but we did
	 * not memorize which entries go to which group. So determine that now.
	 */

	/* Allocate vectors for results */
	v->spl_left = (OffsetNumber *) palloc(nentries * sizeof(OffsetNumber));
	v->spl_right = (OffsetNumber *) palloc(nentries * sizeof(OffsetNumber));
	v->spl_nleft = 0;
	v->spl_nright = 0;

	/* Allocate bounding boxes of left and right groups */
	leftBox = palloc0(VARSIZE(context.boundingBox));
	rightBox = palloc0(VARSIZE(context.boundingBox));

	/*
	 * Allocate an array for "common entries" - entries which can be placed to
	 * either group without affecting overlap along selected axis.
	 */
	commonEntriesCount = 0;
	commonEntries = (CommonEntry *) palloc(nentries * sizeof(CommonEntry));

	/* Helper macros to place an entry in the left or right group */
#define PLACE_LEFT(box, off)					\
	do {										\
		if (v->spl_nleft > 0)					\
			leftBox = adjust_box(leftBox, box);	\
		else									\
			*leftBox = *(box);					\
		v->spl_left[v->spl_nleft++] = off;		\
	} while(0)

#define PLACE_RIGHT(box, off)					\
	do {										\
		if (v->spl_nright > 0)					\
			rightBox = adjust_box(rightBox, box);\
		else									\
			*rightBox = *(box);					\
		v->spl_right[v->spl_nright++] = off;	\
	} while(0)

	/*
	 * Distribute entries which can be distributed unambiguously, and collect
	 * common entries.
	 */
	for (i = FirstOffsetNumber; i <= maxoff; i = OffsetNumberNext(i))
	{
		double		lower,
					upper;

		/*
		 * Get upper and lower bounds along selected axis.
		 */
		box = DatumGetNDBOX(entryvec->vector[i].key);
		lower = box->x[context.dim];
		upper = box->x[context.dim + box->dim];

		if (upper <= context.leftUpper)
		{
			/* Fits to the left group */
			if (lower >= context.rightLower)
			{
				/* Fits also to the right group, so "common entry" */
				commonEntries[commonEntriesCount++].index = i;
			}
			else
			{
				/* Doesn't fit to the right group, so join to the left group */
				PLACE_LEFT(box, i);
			}
		}
		else
		{
			/*
			 * Each entry should fit on either left or right group. Since this
			 * entry didn't fit on the left group, it better fit in the right
			 * group.
			 */
			Assert(lower >= context.rightLower);

			/* Doesn't fit to the left group, so join to the right group */
			PLACE_RIGHT(box, i);
		}
	}

	/*
	 * Distribute "common entries", if any.
	 */
	if (commonEntriesCount > 0)
	{
		/*
		 * Calculate minimum number of entries that must be placed in both
		 * groups, to reach LIMIT_RATIO.
		 */
		int			m = ceil(LIMIT_RATIO * (double) nentries);

		/*
		 * Calculate delta between penalties of join "common entries" to
		 * different groups.
		 */
		for (i = 0; i < commonEntriesCount; i++)
		{
			box = DatumGetNDBOX(entryvec->vector[commonEntries[i].index].key);
			commonEntries[i].delta = Abs(cube_penalty(leftBox, box) - cube_penalty(rightBox, box));
		}

		/*
		 * Sort "common entries" by calculated deltas in order to distribute
		 * the most ambiguous entries first.
		 */
		qsort(commonEntries, commonEntriesCount, sizeof(CommonEntry), common_entry_cmp);

		/*
		 * Distribute "common entries" between groups.
		 */
		for (i = 0; i < commonEntriesCount; i++)
		{
			box = DatumGetNDBOX(entryvec->vector[commonEntries[i].index].key);

			/*
			 * Check if we have to place this entry in either group to achieve
			 * LIMIT_RATIO.
			 */
			if (v->spl_nleft + (commonEntriesCount - i) <= m)
				PLACE_LEFT(box, commonEntries[i].index);
			else if (v->spl_nright + (commonEntriesCount - i) <= m)
				PLACE_RIGHT(box, commonEntries[i].index);
			else
			{
				/* Otherwise select the group by minimal penalty */
				if (cube_penalty(leftBox, box) < cube_penalty(rightBox, box))
					PLACE_LEFT(box, commonEntries[i].index);
				else
					PLACE_RIGHT(box, commonEntries[i].index);
			}
		}
	}

	v->spl_ldatum = PointerGetDatum(leftBox);
	v->spl_rdatum = PointerGetDatum(rightBox);

	return;
}

static void
fallback_split(GistEntryVector *entryvec, GIST_SPLITVEC *v)
{
	OffsetNumber i,
				maxoff;
	NDBOX	   *unionL = NULL,
			   *unionR = NULL;
	int			nbytes;

	maxoff = entryvec->n - 1;

	nbytes = (maxoff + 2) * sizeof(OffsetNumber);
	v->spl_left = (OffsetNumber *) palloc(nbytes);
	v->spl_right = (OffsetNumber *) palloc(nbytes);
	v->spl_nleft = v->spl_nright = 0;

	for (i = FirstOffsetNumber; i <= maxoff; i = OffsetNumberNext(i))
	{
		NDBOX		   *cur = DatumGetNDBOX(entryvec->vector[i].key);

		if (i <= (maxoff - FirstOffsetNumber + 1) / 2)
		{
			v->spl_left[v->spl_nleft] = i;
			if (unionL == NULL)
			{
				unionL = (NDBOX *) palloc(VARSIZE(cur));
				*unionL = *cur;
			}
			else
				unionL = adjust_box(unionL, cur);

			v->spl_nleft++;
		}
		else
		{
			v->spl_right[v->spl_nright] = i;
			if (unionR == NULL)
			{
				unionR = (NDBOX *) palloc(VARSIZE(cur));
				*unionR = *cur;
			}
			else
				unionR = adjust_box(unionR, cur);

			v->spl_nright++;
		}
	}

	v->spl_ldatum = PointerGetDatum(unionL);
	v->spl_rdatum = PointerGetDatum(unionR);
}

/*
 * Consider replacement of currently selected split with the better one.
 */
static inline void
cube_consider_split(ConsiderSplitContext *context, int dimNum,
					 double rightLower, int minLeftCount,
					 double leftUpper, int maxLeftCount)
{
	int			leftCount,
				rightCount;
	float4		ratio,
				overlap;
	double		range;

	/*
	 * Calculate entries distribution ratio assuming most uniform distribution
	 * of common entries.
	 */
	if (minLeftCount >= (context->entriesCount + 1) / 2)
	{
		leftCount = minLeftCount;
	}
	else
	{
		if (maxLeftCount <= context->entriesCount / 2)
			leftCount = maxLeftCount;
		else
			leftCount = context->entriesCount / 2;
	}
	rightCount = context->entriesCount - leftCount;

	/*
	 * Ratio of split - quotient between size of lesser group and total
	 * entries count.
	 */
	ratio = ((float4) Min(leftCount, rightCount)) /
		((float4) context->entriesCount);

	if (ratio > LIMIT_RATIO)
	{
		bool		selectthis = false;

		/*
		 * The ratio is acceptable, so compare current split with previously
		 * selected one. Between splits of one dimension we search for minimal
		 * overlap (allowing negative values) and minimal ration (between same
		 * overlaps. We switch dimension if find less overlap (non-negative)
		 * or less range with same overlap.
		 */
		range = (context->boundingBox)->x[dimNum+(context->boundingBox)->dim]
		 - (context->boundingBox)->x[dimNum];

		overlap = (leftUpper - rightLower) / range;

		/* If there is no previous selection, select this */
		if (context->first)
			selectthis = true;
		else if (context->dim == dimNum)
		{
			/*
			 * Within the same dimension, choose the new split if it has a
			 * smaller overlap, or same overlap but better ratio.
			 */
			if (overlap < context->overlap ||
				(overlap == context->overlap && ratio > context->ratio))
				selectthis = true;
		}
		else
		{
			/*
			 * Across dimensions, choose the new split if it has a smaller
			 * *non-negative* overlap, or same *non-negative* overlap but
			 * bigger range. This condition differs from the one described in
			 * the article. On the datasets where leaf MBRs don't overlap
			 * themselves, non-overlapping splits (i.e. splits which have zero
			 * *non-negative* overlap) are frequently possible. In this case
			 * splits tends to be along one dimension, because most distant
			 * non-overlapping splits (i.e. having lowest negative overlap)
			 * appears to be in the same dimension as in the previous split.
			 * Therefore MBRs appear to be very prolonged along another
			 * dimension, which leads to bad search performance. Using range
			 * as the second split criteria makes MBRs more quadratic. Using
			 * *non-negative* overlap instead of overlap as the first split
			 * criteria gives to range criteria a chance to matter, because
			 * non-overlapping splits are equivalent in this criteria.
			 */
			if (non_negative(overlap) < non_negative(context->overlap) ||
				(range > context->range &&
				 non_negative(overlap) <= non_negative(context->overlap)))
				selectthis = true;
		}

		if (selectthis)
		{
			/* save information about selected split */
			context->first = false;
			context->ratio = ratio;
			context->range = range;
			context->overlap = overlap;
			context->rightLower = rightLower;
			context->leftUpper = leftUpper;
			context->dim = dimNum;
		}
	}
}
