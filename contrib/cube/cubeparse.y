%{
/* contrib/cube/cubeparse.y */

/* NdBox = [(lowerleft),(upperright)] */
/* [(xLL(1)...xLL(N)),(xUR(1)...xUR(n))] */

#define YYSTYPE char *
#define YYDEBUG 1

#include "postgres.h"

#include "cubedata.h"

/*
 * Bison doesn't allocate anything that needs to live across parser calls,
 * so we can easily have it use palloc instead of malloc.  This prevents
 * memory leaks if we error out during parsing.  Note this only works with
 * bison >= 2.0.  However, in bison 1.875 the default is to use alloca()
 * if possible, so there's not really much problem anyhow, at least if
 * you're building with gcc.
 */
#define YYMALLOC palloc
#define YYFREE   pfree

extern int cube_yylex(void);

static char *scanbuf;
static int	scanbuflen;

extern int	cube_yyparse(NDBOX **result);
extern void cube_yyerror(NDBOX **result, const char *message);

static int delim_count(char *s, char delim);
static int check_dim(char *str1, char *str2);
static bool check_maxdim(int dim, char *str1);
static NDBOX * write_box(int dim, char *str1, char *str2, char *typestr);
static NDBOX * write_point_as_box(char *s, int dim, char *typestr);

%}

/* BISON Declarations */
%parse-param {NDBOX **result}
%expect 0
%name-prefix="cube_yy"

%token CUBEFLOAT O_PAREN C_PAREN O_BRACKET C_BRACKET COMMA TYPMOD
%start box

/* Grammar follows */
%%
box: O_BRACKET paren_list COMMA paren_list C_BRACKET
	{
		int dim;

		if ( (dim = check_dim($2, $4)) && check_maxdim(dim, $2) )
			*((void **)result) = write_box(dim, $2, $4, "");
		else
			YYABORT;
	}

	| O_BRACKET paren_list COMMA paren_list C_BRACKET TYPMOD
	{
		int dim;

		if ( (dim = check_dim($2, $4)) && check_maxdim(dim, $2) )
			*((void **)result) = write_box(dim, $2, $4, $6);
		else
			YYABORT;
	}

	| paren_list COMMA paren_list
	{
		int dim;

		if ( (dim = check_dim($1, $3)) && check_maxdim(dim, $1) )
			*((void **)result) = write_box(dim, $1, $3, "");
		else
			YYABORT;
	}

	| paren_list COMMA paren_list TYPMOD
	{
		int dim;

		if ( (dim = check_dim($1, $3)) && check_maxdim(dim, $1) )
			*((void **)result) = write_box(dim, $1, $3, $4);
		else
			YYABORT;
	}

	| paren_list
	{
		int dim;

		dim = delim_count($1, ',') + 1;
		if (!check_maxdim(dim, $1))
			YYABORT;

		*((void **)result) = write_point_as_box($1, dim, "");
	}

	| paren_list TYPMOD
	{
		int dim;

		dim = delim_count($1, ',') + 1;
		if (!check_maxdim(dim, $1))
			YYABORT;

		*((void **)result) = write_point_as_box($1, dim, $2);
	}

	| list
	{
		int dim;

		dim = delim_count($1, ',') + 1;
		if (!check_maxdim(dim, $1))
			YYABORT;

		*((void **)result) = write_point_as_box($1, dim, "");
	}

	| list TYPMOD
	{
		int dim;

		dim = delim_count($1, ',') + 1;
		if (!check_maxdim(dim, $1))
			YYABORT;

		*((void **)result) = write_point_as_box($1, dim, $2);
	}
	;

paren_list: O_PAREN list C_PAREN
	{
		$$ = $2;
	}
	;

list: CUBEFLOAT
	{
		/* alloc enough space to be sure whole list will fit */
		$$ = palloc(scanbuflen + 1);
		strcpy($$, $1);
	}
	| list COMMA CUBEFLOAT
	{
		$$ = $1;
		strcat($$, ",");
		strcat($$, $3);
	}
	;

%%

static int
delim_count(char *s, char delim)
{
	int			ndelim = 0;

	while ((s = strchr(s, delim)) != NULL)
	{
		ndelim++;
		s++;
	}
	return (ndelim);
}

static int
check_dim(char *str1, char *str2)
{
	int dim;

	dim = delim_count(str1, ',') + 1;
	if ((delim_count(str2, ',') + 1) != dim)
	{
		ereport(ERROR,
				(errcode(ERRCODE_SYNTAX_ERROR),
				 errmsg("bad cube representation"),
				 errdetail("Different point dimensions in (%s) and (%s).",
						   str1, str2)));
		return 0;
	} else 
		return dim;
}

static bool
check_maxdim(int dim, char *str1)
{
	if (dim > CUBE_MAX_DIM) {
		ereport(ERROR,
				(errcode(ERRCODE_SYNTAX_ERROR),
				 errmsg("bad cube representation"),
				 errdetail("A cube cannot have more than %d dimensions.",
						   CUBE_MAX_DIM)));
		return 0;
	} else 
		return 1;
}

static NDBOX *
write_box(int dim, char *str1, char *str2, char *typestr)
{
	NDBOX	   *bp;
	char	   *s;
	int			i, type;
	double		x;
	bool		point = true;

	if (strcmp(typestr, ":f4") == 0)
		type = CUBE_FLOAT4;
	else if (strcmp(typestr, ":i4") == 0)
		type = CUBE_INT4;
	else if (strcmp(typestr, ":i2") == 0)
		type = CUBE_INT2;
	else if (strcmp(typestr, ":i1") == 0)
		type = CUBE_INT1;
	else
		type = CUBE_FLOAT8;

	bp = init_cube(dim, 0, type);

	s = str1;
	x = strtod(s, NULL);
	check_coord(x, type);
	set_coord(bp, i=0, x);
	while ((s = strchr(s, ',')) != NULL)
	{
		s++; i++;
		x = strtod(s, NULL);
		check_coord(x, type);
		set_coord(bp, i, x);
	}

	s = str2;
	x = strtod(s, NULL);
	check_coord(x, type);
	set_coord(bp, i=dim, x);
	while ((s = strchr(s, ',')) != NULL)
	{
		s++; i++;
		x = strtod(s, NULL);
		check_coord(x, type);
		set_coord(bp, i, x);
		if (LL_COORD(bp, i-dim) != UR_COORD(bp, i-dim))
			point = false;
	}

	if (LL_COORD(bp, 0) != UR_COORD(bp, 0))
		point = false;

	if (point)
		cube_to_point(bp);

	return(bp);
}

static NDBOX *
write_point_as_box(char *str, int dim, char *typestr)
{
	NDBOX		*bp;
	int			i, type;
	double		x;
	char		*s = str;

	if (strcmp(typestr, ":f4") == 0)
		type = CUBE_FLOAT4;
	else if (strcmp(typestr, ":i4") == 0)
		type = CUBE_INT4;
	else if (strcmp(typestr, ":i2") == 0)
		type = CUBE_INT2;
	else if (strcmp(typestr, ":i1") == 0)
		type = CUBE_INT1;
	else
		type = CUBE_FLOAT8;

	bp = init_cube(dim, 1, type);

	i = 0;
	x = strtod(s, NULL);
	check_coord(x, type);
	set_coord(bp, 0, x);
	while ((s = strchr(s, ',')) != NULL)
	{
		s++; i++;
		x = strtod(s, NULL);
		check_coord(x, type);
		set_coord(bp, i, x);
	}

	return(bp);
}

#include "cubescan.c"
