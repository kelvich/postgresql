/* contrib/cube/cubedata.h */

#define CUBE_MAX_DIM (100)

typedef struct NDBOX
{
	/* varlena header (do not touch directly!) */
	int32 vl_len_;

	/* 
	 * Header contains info about NDBOX. For binary
	 * compatibility with old versions it is defined
	 * as uint32.
	 * 
	 * Following information is stored:
	 *
	 *  bits 0-7  : number of cube dimensions;
	 *  bits 8-30 : not used;
	 *  bit  31   : point flag. If set, then NDBOX stores
	 *             n dimensions instead of 2*n;
	 */
	unsigned int header;
	double x[1];
} NDBOX;

#define DatumGetNDBOX(x)   ((NDBOX*)DatumGetPointer(x))
#define PG_GETARG_NDBOX(x) DatumGetNDBOX( PG_DETOAST_DATUM(PG_GETARG_DATUM(x)) )
#define PG_RETURN_NDBOX(x) PG_RETURN_POINTER(x)

#define IS_POINT(cube)		( cube->header >> 31 )
#define SET_POINT_BIT(cube)	( cube->header = cube->header + 0x80000000 )
#define DIM(cube)			( cube->header & 0x7fffffff )
#define SET_DIM(cube, _dim)	( cube->header = _dim )

#define LL_COORD(cube, i) ( cube->x[i] )
#define UR_COORD(cube, i) ( IS_POINT(cube) ? cube->x[i] : cube->x[i + DIM(cube)] )

#define POINT_SIZE(_dim) (offsetof(NDBOX, x[0]) + sizeof(double)*_dim)
#define CUBE_SIZE(_dim) (offsetof(NDBOX, x[0]) + sizeof(double)*_dim*2)
