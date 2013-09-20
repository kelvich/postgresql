/* contrib/cube/cubedata.h */

#define CUBE_MAX_DIM (100)

typedef struct NDBOX
{
	int32		vl_len_;		/* varlena header (do not touch directly!) */
	unsigned int dim;
	double		x[1];
} NDBOX;

#define DatumGetNDBOX(x)	((NDBOX*)DatumGetPointer(x))
#define PG_GETARG_NDBOX(x)	DatumGetNDBOX( PG_DETOAST_DATUM(PG_GETARG_DATUM(x)) )
#define PG_RETURN_NDBOX(x)	PG_RETURN_POINTER(x)

/*
 * Represents information about an entry that can be placed to either group
 * without affecting overlap over selected axis ("common entry").
 */
typedef struct
{
	/* Index of entry in the initial array */
	int			index;
	/* Delta between penalties of entry insertion into different groups */
	double		delta;
} CommonEntry;

/*
 * Context for g_box_consider_split. Contains information about currently
 * selected split and some general information.
 */
typedef struct
{
	int			entriesCount;	/* total number of entries being split */
	NDBOX	   *boundingBox;	/* minimum bounding box across all entries */

	/* Information about currently selected split follows */

	bool		first;			/* true if no split was selected yet */

	double		leftUpper;		/* upper bound of left interval */
	double		rightLower;		/* lower bound of right interval */

	float4		ratio;
	float4		overlap;
	int			dim;			/* axis of this split */
	double		range;			/* width of general MBR projection to the
								 * selected axis */
} ConsiderSplitContext;

/*
 * Interval represents projection of box to axis.
 */
typedef struct
{
	double		lower,
				upper;
} SplitInterval;

/* Minimum accepted ratio of split */
#define LIMIT_RATIO 0.3

/* 
 * We use Guttman split when cube dim >= SPLIT_THRESHOLD 
 * and Korotkov split otherwise.
 */
#define SPLIT_THRESHOLD 6

