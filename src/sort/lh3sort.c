/*
 * sort - sort lines of text (with all kinds of options). Copyright (C) 88,
 * 91, 92, 93, 94, 95, 1996 Free Software Foundation
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 * 
 * Written December 1988 by Mike Haertel. The author may be reached (Email) at
 * the address mike@gnu.ai.mit.edu, or (US mail) as Mike Haertel c/o Free
 * Software Foundation.
 */

/*
  This program is modified from sort.c included in GNU textutils 1.22. I
  added mix_numeric comparison which is useful for sorting biological
  data. I also replaced the old top-down recursive mergesort with an
  bottom-up iterative mergesort. On simple test file, the new code can
  sometimes one time faster, depending on the data in use.

  Edited by Heng Li on 11 April 2008.

  Note: The heap merge makes this `sort' implementation not stable.
*/

#include "config.h"

/* Get isblank from GNU libc.  */
#define _GNU_SOURCE

#include <sys/types.h>
#include <signal.h>
#include <stdio.h>
#ifndef ENABLE_ASSERTIONS
#define NDEBUG 1
#endif
#include <assert.h>
#include "system.h"
#include "error.h"

#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#ifndef UCHAR_MAX
#define UCHAR_MAX 255
#endif
#endif
#ifndef STDC_HEADERS
char           *malloc();
char           *realloc();
void            free();
#endif

/* Undefine, to avoid warning about redefinition on some systems.  */
#undef min
#define min(a, b) ((a) < (b) ? (a) : (b))

#define UCHAR_LIM (UCHAR_MAX + 1)
#define UCHAR(c) ((unsigned char) (c))

#ifndef DEFAULT_TMPDIR
#define DEFAULT_TMPDIR "/tmp"
#endif

/*
 * Use this as exit status in case of error, not EXIT_FAILURE.  This is
 * necessary because EXIT_FAILURE is usually 1 and POSIX requires that sort
 * exit with status 1 IFF invoked with -c and the input is not properly
 * sorted.  Any other irregular exit must exit with a status code greater
 * than 1.
 */
#define SORT_FAILURE 2

/* The kind of blanks for '-b' to skip in various options. */
enum blanktype {
	bl_start, bl_end, bl_both
};

/* The character marking end of line. Default to \n. */
int eolchar = '\n';

/* Lines are held in core as counted strings. */
struct line {
	char           *text;	/* Text of the line. */
	int             length;	/* Length not including final newline. */
	char           *keybeg;	/* Start of first key. */
	char           *keylim;	/* Limit of first key. */
};

/* Arrays of lines. */
struct lines {
	struct line    *lines;	/* Dynamically allocated array of lines. */
	int             used;	/* Number of slots used. */
	int             alloc;	/* Number of slots allocated. */
	int             limit;	/* Max number of slots to allocate.  */
};

/* Input buffers. */
struct buffer {
	char           *buf;	/* Dynamically allocated buffer. */
	int             used;	/* Number of bytes used. */
	int             alloc;	/* Number of bytes allocated. */
	int             left;	/* Number of bytes left after line parsing. */
};

struct keyfield {
	int             sword; /* Zero-origin 'word' to start at. */
	int             schar; /* Additional characters to skip. */
	int             eword; /* Zero-origin first word after field. */
	int             echar; /* Additional characters in field. */
	int  const     *ignore;/* Boolean array of characters to ignore. */
	char const     *translate;	/* Translation applied to characters. */
	int             skipsblanks; /* Skip leading white space at start. */
	int             skipeblanks; /* Skip trailing white space at finish. */
	int             mixed_numeric;
	int             numeric; /* Flag for numeric comparison. Handle strings of digits with optional decimal
				              * point , but no exponential notation. */
	int             general_numeric; /* Flag for general, numeric comparison. Handle numbers
						              * in exponential notation. */
	int             month;	/* Flag for comparison by month name. */
	int             reverse;/* Reverse the sense of comparison. */
	struct keyfield *next;	/* Next keyfield to try. */
};

struct month {
	char const     *name;
	int             val;
};

/* The name this program was run with. */
char           *program_name;

/* Table of white space. */
static int      blanks[UCHAR_LIM];

/* Table of non-printing characters. */
static int      nonprinting[UCHAR_LIM];

/* Table of non-dictionary characters (not letters, digits, or blanks). */
static int      nondictionary[UCHAR_LIM];

/* Translation table folding lower case to upper. */
static char     fold_toupper[UCHAR_LIM];

/*
 * Table mapping 3-letter month names to integers. Alphabetic order allows
 * binary search.
 */
static struct month const monthtab[] =
{
	{"APR", 4},
	{"AUG", 8},
	{"DEC", 12},
	{"FEB", 2},
	{"JAN", 1},
	{"JUL", 7},
	{"JUN", 6},
	{"MAR", 3},
	{"MAY", 5},
	{"NOV", 11},
	{"OCT", 10},
	{"SEP", 9}
};

/* During the merge phase, the number of files to merge at once. */
#define NMERGE 16

/*
 * Initial buffer size for in core sorting.  Will not grow unless a line
 * longer than this is seen.
 */
static int sortalloc = 512 * 1024;

/*
 * Initial buffer size for in core merge buffers.  Bear in mind that up to
 * NMERGE * mergealloc bytes may be allocated for merge buffers.
 */
static int mergealloc = 16 * 1024;

/* Guess of average line length. */
static int linelength = 30;

/* Maximum number of elements for the array(s) of struct line's, in bytes.  */
#define LINEALLOC (256 * 1024)

/* Prefix for temporary file names. */
static char *temp_file_prefix;

/* Flag to reverse the order of all comparisons. */
static int  reverse;

/*
 * Flag for stable sort.  This turns off the last ditch bytewise comparison
 * of lines, and instead leaves lines in the same order they were read if all
 * keys compare equal.
 */
static int  stable;

/*
 * Tab character separating fields.  If NUL, then fields are separated by the
 * empty string between a non-whitespace character and a whitespace
 * character.
 */
static char tab;

/*
 * Flag to remove consecutive duplicate lines from the output. Only the last
 * of a sequence of equal lines will be output.
 */
static int  unique;

/* Nonzero if any of the input files are the standard input. */
static int  have_read_stdin;

/* Lists of key field comparisons to be tried. */
static struct keyfield keyhead;

static void usage(int status)
{
	if (status != 0)
		fprintf(stderr, _("Try `%s --help' for more information.\n"),
			program_name);
	else {
		printf(_("\
Usage: %s [OPTION]... [FILE]...\n\
"),
		       program_name);
		printf(_("\
Write sorted concatenation of all FILE(s) to standard output.\n\
\n\
  +POS1 [-POS2]    start a key at POS1, end it before POS2\n\
  -b               ignore leading blanks in sort fields or keys\n\
  -c               check if given files already sorted, do not sort\n\
  -d               consider only [a-zA-Z0-9] characters in keys\n\
  -f               fold lower case to upper case characters in keys\n\
  -g               compare according to general numerical value, imply -b\n\
  -i               consider only [\\040-\\0176] characters in keys\n\
  -k POS1[,POS2]   same as +POS1 [-POS2], but all positions counted from 1\n\
  -m               merge already sorted files, do not sort\n\
  -M               compare (unknown) < `JAN' < ... < `DEC', imply -b\n\
  -n               compare according to string numerical value, imply -b\n\
  -N               mixed string and numerical value, imply -b\n\
  -o FILE          write result on FILE instead of standard output\n\
  -r               reverse the result of comparisons\n\
  -s               stabilize sort by disabling last resort comparison\n\
  -t SEP           use SEParator instead of non- to whitespace transition\n\
  -T DIRECT        use DIRECT for temporary files, not $TMPDIR or %s\n\
  -u               with -c, check for strict ordering;\n\
                   with -m, only output the first of an equal sequence\n\
  -z               end lines with 0 byte, not newline, for find -print0\n\
  -h               display this help and exit\n\
\n\
POS is F[.C][OPTS], where F is the field number and C the character\n\
position in the field, both counted from zero. OPTS is made up of one or\n\
more of NMbdfingr, this effectively disable global -NMbdfingr settings\n\
for that key. If no key given, use the entire line as key. With no\n\
FILE, or when FILE is -, read standard input.\n\
")
		       ,DEFAULT_TMPDIR);
		puts(_("\n\
This program is a modified version of `sort' from GNU textutils-1.22. The\n\
original author is Mike Haertel. Heng Li improves the sorting codes and\n\
adds mixed string-numeric comparison. However, this `sort' does not perform\n\
a stable sort any more."));
	}
	/*
	 * Don't use EXIT_FAILURE here in case it is defined to be 1. POSIX
	 * requires that sort return 1 IFF invoked with -c and the input is
	 * not properly sorted.
	 */
	assert(status == 0 || status == SORT_FAILURE);
	exit(status);
}

/* The list of temporary files. */
static struct tempnode {
	char           *name;
	struct tempnode *next;
} temphead;

/* Clean up any remaining temporary files. */

static void cleanup(void)
{
	struct tempnode *node;
	for (node = temphead.next; node; node = node->next)
		unlink(node->name);
}

/* Allocate N bytes of memory dynamically, with error checking.  */

static char *xmalloc(unsigned int n)
{
	char *p;
	p = malloc(n);
	if (p == 0) {
		error(0, 0, _("virtual memory exhausted"));
		cleanup();
		exit(SORT_FAILURE);
	}
	return p;
}

/*
 * Change the size of an allocated block of memory P to N bytes, with error
 * checking. If P is NULL, run xmalloc. If N is 0, run free and return NULL.
 */

static char *xrealloc(char *p, unsigned int n)
{
	if (p == 0) return xmalloc(n);
	if (n == 0) {
		free(p);
		return 0;
	}
	p = realloc(p, n);
	if (p == 0) {
		error(0, 0, _("virtual memory exhausted"));
		cleanup();
		exit(SORT_FAILURE);
	}
	return p;
}

static FILE *xtmpfopen(const char *file)
{
	FILE *fp;
	int fd;
	fd = open(file, O_WRONLY | O_CREAT | O_TRUNC, 0600);
	if (fd < 0 || (fp = fdopen(fd, "w")) == NULL) {
		error(0, errno, "%s", file);
		cleanup();
		exit(SORT_FAILURE);
	}
	return fp;
}

static FILE *xfopen(const char *file, const char *how)
{
	FILE *fp;
	if (strcmp(file, "-") == 0) fp = stdin;
	else {
		if ((fp = fopen(file, how)) == NULL) {
			error(0, errno, "%s", file);
			cleanup();
			exit(SORT_FAILURE);
		}
	}
	if (fp == stdin) have_read_stdin = 1;
	return fp;
}

static void xfclose(FILE * fp)
{
	if (fp == stdin) {
		/* Allow reading stdin from tty more than once. */
		if (feof(fp)) clearerr(fp);
	} else if (fp == stdout) {
		if (fflush(fp) != 0) {
			error(0, errno, _("flushing file"));
			cleanup();
			exit(SORT_FAILURE);
		}
	} else {
		if (fclose(fp) != 0) {
			error(0, errno, _("error closing file"));
			cleanup();
			exit(SORT_FAILURE);
		}
	}
}

static void write_bytes(const char *buf, size_t n_bytes, FILE * fp)
{
	if (fwrite(buf, 1, n_bytes, fp) != n_bytes) {
		error(0, errno, _("write error"));
		cleanup();
		exit(SORT_FAILURE);
	}
}

/* Return a name for a temporary file. */

static char *tempname(void)
{
	static unsigned int seq;
	int             len = strlen(temp_file_prefix);
	char           *name = xmalloc(len + 1 + sizeof("sort") - 1 + 5 + 5 + 1);
	struct tempnode *node;

	node = (struct tempnode *)xmalloc(sizeof(struct tempnode));
	sprintf(name, "%s%ssort%5.5d%5.5d",
		temp_file_prefix,
		(len && temp_file_prefix[len - 1] != '/') ? "/" : "",
		(unsigned int) getpid() & 0xffff, seq);

	/* Make sure that SEQ's value fits in 5 digits.  */
	++seq;
	if (seq >= 100000) seq = 0;

	node->name = name;
	node->next = temphead.next;
	temphead.next = node;
	return name;
}

/*
 * Search through the list of temporary files for NAME; remove it if it is
 * found on the list.
 */

static void zaptemp(char *name)
{
	struct tempnode *node, *temp;
	for (node = &temphead; node->next; node = node->next)
		if (!strcmp(name, node->next->name))
			break;
	if (node->next) {
		temp = node->next;
		unlink(temp->name);
		free(temp->name);
		node->next = temp->next;
		free((char *) temp);
	}
}

/* Initialize the character class tables. */

static void inittables(void)
{
	int i;
	for (i = 0; i < UCHAR_LIM; ++i) {
		if (ISBLANK(i)) blanks[i] = 1;
		if (!ISPRINT(i)) nonprinting[i] = 1;
		if (!ISALNUM(i) && !ISBLANK(i)) nondictionary[i] = 1;
		if (ISLOWER(i)) fold_toupper[i] = toupper(i);
		else fold_toupper[i] = i;
	}
}

/* Initialize BUF, allocating ALLOC bytes initially. */

static void initbuf(struct buffer *buf, int alloc)
{
	buf->alloc = alloc;
	buf->buf = xmalloc(buf->alloc);
	buf->used = buf->left = 0;
}

/*
 * Fill BUF reading from FP, moving buf->left bytes from the end of buf->buf
 * to the beginning first.	If EOF is reached and the file wasn't
 * terminated by a newline, supply one.  Return a count of bytes buffered.
 */

static int fillbuf(struct buffer *buf, FILE *fp)
{
	int cc;

	memmove(buf->buf, buf->buf + buf->used - buf->left, buf->left);
	buf->used = buf->left;

	while (!feof(fp) && (buf->used == 0 || !memchr(buf->buf, eolchar, buf->used))) {
		if (buf->used == buf->alloc) {
			buf->alloc *= 2;
			buf->buf = xrealloc(buf->buf, buf->alloc);
		}
		cc = fread(buf->buf + buf->used, 1, buf->alloc - buf->used, fp);
		if (ferror(fp)) {
			error(0, errno, _("read error"));
			cleanup();
			exit(SORT_FAILURE);
		}
		buf->used += cc;
	}

	if (feof(fp) && buf->used && buf->buf[buf->used - 1] != eolchar) {
		if (buf->used == buf->alloc) {
			buf->alloc *= 2;
			buf->buf = xrealloc(buf->buf, buf->alloc);
		}
		buf->buf[buf->used++] = eolchar;
	}
	return buf->used;
}

/*
 * Initialize LINES, allocating space for ALLOC lines initially. LIMIT is the
 * maximum possible number of lines to allocate space for, ever.
 */

static void initlines(struct lines *lines, int alloc, int limit)
{
	lines->alloc = alloc;
	lines->lines = (struct line *) xmalloc(lines->alloc * sizeof(struct line));
	lines->used = 0;
	lines->limit = limit;
}

/*
 * Return a pointer to the first character of the field specified by KEY in
 * LINE.
 */

static char *begfield(const struct line *line, const struct keyfield *key)
{
	register char  *ptr = line->text, *lim = ptr + line->length;
	register int sword = key->sword, schar = key->schar;

	if (tab) {
		while (ptr < lim && sword--) {
			while (ptr < lim && *ptr != tab) ++ptr;
			if (ptr < lim) ++ptr;
		}
	} else {
		while (ptr < lim && sword--) {
			while (ptr < lim && blanks[UCHAR(*ptr)]) ++ptr;
			while (ptr < lim && !blanks[UCHAR(*ptr)]) ++ptr;
		}
	}

	if (key->skipsblanks)
		while (ptr < lim && blanks[UCHAR(*ptr)]) ++ptr;

	if (ptr + schar <= lim) ptr += schar;
	else ptr = lim;

	return ptr;
}

/*
 * Return the limit of (a pointer to the first character after) the field in
 * LINE specified by KEY.
 */

static char *limfield(const struct line * line, const struct keyfield * key)
{
	register char  *ptr = line->text, *lim = ptr + line->length;
	register int eword = key->eword, echar = key->echar;

	/*
	 * Note: from the POSIX spec: The leading field separator itself is
	 * included in a field when -t is not used.  FIXME: move this comment
	 * up...
	 */

	/*
	 * Move PTR past EWORD fields or to one past the last byte on LINE,
	 * whichever comes first.  If there are more than EWORD fields, leave
	 * PTR pointing at the beginning of the field having zero-based
	 * index, EWORD.  If a delimiter character was specified (via -t),
	 * then that `beginning' is the first character following the
	 * delimiting TAB. Otherwise, leave PTR pointing at the first `blank'
	 * character after the preceding field.
	 */
	if (tab) {
		while (ptr < lim && eword--) {
			while (ptr < lim && *ptr != tab) ++ptr;
			if (ptr < lim && (eword || echar > 0)) ++ptr;
		}
	} else {
		while (ptr < lim && eword--) {
			while (ptr < lim && blanks[UCHAR(*ptr)]) ++ptr;
			while (ptr < lim && !blanks[UCHAR(*ptr)]) ++ptr;
		}
	}

#ifdef POSIX_UNSPECIFIED
	/*
	 * The following block of code makes GNU sort incompatible with
	 * standard Unix sort, so it's ifdef'd out for now. The POSIX spec
	 * isn't clear on how to interpret this. FIXME: request
	 * clarification.
	 * 
	 * From: kwzh@gnu.ai.mit.edu (Karl Heuer) Date: Thu, 30 May 96 12:20:41
	 * -0400
	 * 
	 * [...]I believe I've found another bug in `sort'.
	 * 
	 * $ cat /tmp/sort.in a b c 2 d pq rs 1 t $ textutils-1.15/src/sort +0.6
	 * -0.7 </tmp/sort.in a b c 2 d pq rs 1 t $ /bin/sort +0.6 -0.7
	 * </tmp/sort.in pq rs 1 t a b c 2 d
	 * 
	 * Unix sort produced the answer I expected: sort on the single
	 * character in column 6.  GNU sort produced different results,
	 * because it disagrees on the interpretation of the key-end spec
	 * "-M.N".  Unix sort reads this as "skip M fields, then N
	 * characters"; but GNU sort wants it to mean "skip M fields, then
	 * either N characters or the rest of the current field, whichever
	 * comes first".  This extra clause applies only to key-ends, not
	 * key-starts.
	 */

	/* Make LIM point to the end of (one byte past) the current field.  */
	if (tab) {
		char           *newlim;
		newlim = memchr(ptr, tab, lim - ptr);
		if (newlim) lim = newlim;
	} else {
		char           *newlim;
		newlim = ptr;
		while (newlim < lim && blanks[UCHAR(*newlim)]) ++newlim;
		while (newlim < lim && !blanks[UCHAR(*newlim)]) ++newlim;
		lim = newlim;
	}
#endif

	/*
	 * If we're skipping leading blanks, don't start counting characters
	 * until after skipping past any leading blanks.
	 */
	if (key->skipsblanks)
		while (ptr < lim && blanks[UCHAR(*ptr)]) ++ptr;

	/* Advance PTR by ECHAR (if possible), but no further than LIM.  */
	if (ptr + echar <= lim) ptr += echar;
	else ptr = lim;

	return ptr;
}

/* FIXME */

void trim_trailing_blanks(const char *a_start, char **a_end)
{
	while (*a_end > a_start && blanks[UCHAR(*(*a_end - 1))]) --(*a_end);
}

/*
 * Find the lines in BUF, storing pointers and lengths in LINES. Also replace
 * newlines in BUF with NULs.
 */

static void findlines(struct buffer * buf, struct lines * lines)
{
	register char  *beg = buf->buf, *lim = buf->buf + buf->used, *ptr;
	struct keyfield *key = keyhead.next;

	lines->used = 0;

	while (beg < lim && (ptr = memchr(beg, eolchar, lim - beg))
	       && lines->used < lines->limit) {
		/*
		 * There are various places in the code that rely on a NUL
		 * being at the end of in-core lines; NULs inside the lines
		 * will not cause trouble, though.
		 */
		*ptr = '\0';

		if (lines->used == lines->alloc) {
			lines->alloc *= 2;
			lines->lines = (struct line *)xrealloc((char *) lines->lines, lines->alloc * sizeof(struct line));
		}
		lines->lines[lines->used].text = beg;
		lines->lines[lines->used].length = ptr - beg;

		/* Precompute the position of the first key for efficiency. */
		if (key) {
			if (key->eword >= 0) lines->lines[lines->used].keylim = limfield(&lines->lines[lines->used], key);
			else lines->lines[lines->used].keylim = ptr;

			if (key->sword >= 0) lines->lines[lines->used].keybeg = begfield(&lines->lines[lines->used], key);
			else {
				if (key->skipsblanks)
					while (blanks[UCHAR(*beg)]) ++beg;
				lines->lines[lines->used].keybeg = beg;
			}
			if (key->skipeblanks)
				trim_trailing_blanks(lines->lines[lines->used].keybeg, &lines->lines[lines->used].keylim);
		} else {
			lines->lines[lines->used].keybeg = 0;
			lines->lines[lines->used].keylim = 0;
		}

		++lines->used;
		beg = ptr + 1;
	}

	buf->left = lim - beg;
}

/*
 * Compare strings A and B containing decimal fractions < 1.  Each string
 * should begin with a decimal point followed immediately by the digits of
 * the fraction.  Strings not of this form are considered to be zero.
 */

static int fraccompare(register const char *a, register const char *b)
{
	register int    tmpa = *a;
	register int    tmpb = *b;

	if (tmpa == '.' && tmpb == '.') {
		do tmpa = *++a, tmpb = *++b; while (tmpa == tmpb && ISDIGIT(tmpa));
		if (ISDIGIT(tmpa) && ISDIGIT(tmpb)) return tmpa - tmpb;
		if (ISDIGIT(tmpa)) {
			while (tmpa == '0') tmpa = *++a;
			if (ISDIGIT(tmpa)) return 1;
			return 0;
		}
		if (ISDIGIT(tmpb)) {
			while (tmpb == '0') tmpb = *++b;
			if (ISDIGIT(tmpb)) return -1;
			return 0;
		}
		return 0;
	} else if (tmpa == '.') {
		do tmpa = *++a; while (tmpa == '0');
		if (ISDIGIT(tmpa)) return 1;
		return 0;
	} else if (tmpb == '.') {
		do tmpb = *++b; while (tmpb == '0');
		if (ISDIGIT(tmpb)) return -1;
		return 0;
	}
	return 0;
}

/*
 * Compare strings A and B as numbers without explicitly converting them to
 * machine numbers.  Comparatively slow for short strings, but asymptotically
 * hideously fast.
 */

static int numcompare(register const char *a, register const char *b)
{
	register int    tmpa, tmpb, loga, logb, tmp;

	tmpa = UCHAR(*a);
	tmpb = UCHAR(*b);

	while (blanks[tmpa]) tmpa = UCHAR(*++a);
	while (blanks[tmpb]) tmpb = UCHAR(*++b);

	if (tmpa == '-') {
		do tmpa = *++a; while (tmpa == '0');
		if (tmpb != '-') {
			if (tmpa == '.')
				do tmpa = *++a; while (tmpa == '0');
			if (ISDIGIT(tmpa)) return -1;
			while (tmpb == '0') tmpb = *++b;
			if (tmpb == '.')
				do tmpb = *++b; while (tmpb == '0');
			if (ISDIGIT(tmpb)) return -1;
			return 0;
		}
		do tmpb = *++b; while (tmpb == '0');

		while (tmpa == tmpb && ISDIGIT(tmpa))
			tmpa = *++a, tmpb = *++b;

		if ((tmpa == '.' && !ISDIGIT(tmpb))
		    || (tmpb == '.' && !ISDIGIT(tmpa)))
			return -fraccompare(a, b);

		if (ISDIGIT(tmpa)) for (loga = 1; ISDIGIT(*++a); ++loga);
		else loga = 0;
		if (ISDIGIT(tmpb)) for (logb = 1; ISDIGIT(*++b); ++logb);
		else logb = 0;
		if ((tmp = logb - loga) != 0) return tmp;
		if (!loga) return 0;
		return tmpb - tmpa;
	} else if (tmpb == '-') {
		do tmpb = *++b; while (tmpb == '0');
		if (tmpb == '.')
			do tmpb = *++b; while (tmpb == '0');
		if (ISDIGIT(tmpb)) return 1;
		while (tmpa == '0') tmpa = *++a;
		if (tmpa == '.')
			do tmpa = *++a; while (tmpa == '0');
		if (ISDIGIT(tmpa)) return 1;
		return 0;
	} else {
		while (tmpa == '0') tmpa = *++a;
		while (tmpb == '0') tmpb = *++b;

		while (tmpa == tmpb && ISDIGIT(tmpa))
			tmpa = *++a, tmpb = *++b;

		if ((tmpa == '.' && !ISDIGIT(tmpb))
		    || (tmpb == '.' && !ISDIGIT(tmpa)))
			return fraccompare(a, b);

		if (ISDIGIT(tmpa)) for (loga = 1; ISDIGIT(*++a); ++loga);
		else loga = 0;
		if (ISDIGIT(tmpb)) for (logb = 1; ISDIGIT(*++b); ++logb);
		else logb = 0;
		if ((tmp = loga - logb) != 0) return tmp;
		if (!loga) return 0;
		return tmpa - tmpb;
	}
}

/* lh3: This function is copied from sort in GNU coreutils 6.9. It is more general than the previous version. */
static int general_numcompare(const char *sa, const char *sb)
{
	/* FIXME: add option to warn about failed conversions.  */
	/*
	 * FIXME: maybe add option to try expensive FP conversion only if A
	 * and B can't be compared more cheaply/accurately.
	 */

	char *ea;
	char *eb;
	double a = strtod(sa, &ea);
	double b = strtod(sb, &eb);

	/* Put conversion errors at the start of the collating sequence.  */
	if (sa == ea) return sb == eb ? 0 : -1;
	if (sb == eb) return 1;

	/*
	 * Sort numbers in the usual way, where -0 == +0.  Put NaNs after
	 * conversion errors but before numbers; sort them by internal
	 * bit-pattern, for lack of a more portable alternative.
	 */
	return (a < b ? -1
		: a > b ? 1
		: a == b ? 0
		: b == b ? -1
		: a == a ? 1
		: memcmp((char *) &a, (char *) &b, sizeof a));
}
/* lh3: if a<b and b<c, then we must have a<c (transitivity). */
static int mixed_numcompare(const char *_a, const char *_b)
{
	const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
	const unsigned char *pa = a, *pb = b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			while (*pa == '0') ++pa;
			while (*pb == '0') ++pb;
			while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
			if (isdigit(*pa) && isdigit(*pb)) {
				int i = 0;
				while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
				return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
			} else if (isdigit(*pa)) return 1;
			else if (isdigit(*pb)) return -1;
			else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
		} else {
			if (*pa != *pb) return (int)*pa - (int)*pb;
			++pa; ++pb;
		}
	}
	return *pa? 1 : *pb? -1 : 0;
}
/*
static int mixed_numcompare(const char *a, const char *b)
{
	char *pa, *pb;
	pa = (char*)a; pb = (char*)b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai < bi? -1 : 1;
			if (pa - a != pb - b) return pa - a > pb - b? -1 : 1;
		} else {
			if (*pa != *pb) return *pa < *pb? -1 : 1;
			++pa; ++pb;
		}
	}
	return *pa? 1 : *pb? -1 : 0;
}
*/
/*
 * Return an integer <= 12 associated with month name S with length LEN, 0 if
 * the name in S is not recognized.
 */

static int getmonth(const char *s, int len)
{
	char            month[4];
	register int i, lo = 0, hi = 12;

	while (len > 0 && blanks[UCHAR(*s)])
		++s, --len;

	if (len < 3) return 0;

	for (i = 0; i < 3; ++i) month[i] = fold_toupper[UCHAR(s[i])];
	month[3] = '\0';

	while (hi - lo > 1)
		if (strcmp(month, monthtab[(lo + hi) / 2].name) < 0) hi = (lo + hi) / 2;
		else lo = (lo + hi) / 2;
	if (!strcmp(month, monthtab[lo].name))
		return monthtab[lo].val;
	return 0;
}

/*
 * Compare two lines A and B trying every key in sequence until there are no
 * more keys or a difference is found.
 */

static int keycompare(const struct line *a, const struct line *b)
{
	register char  *texta, *textb, *lima, *limb;
	register unsigned char *translate;
	register int  const *ignore;
	struct keyfield *key;
	int diff = 0, iter = 0, lena, lenb;

	for (key = keyhead.next; key; key = key->next, ++iter) {
		ignore = key->ignore;
		translate = (unsigned char *) key->translate;

		/* Find the beginning and limit of each field. */
		if (iter || a->keybeg == NULL || b->keybeg == NULL) {
			if (key->eword >= 0) lima = limfield(a, key), limb = limfield(b, key);
			else lima = a->text + a->length, limb = b->text + b->length;

			if (key->sword >= 0) texta = begfield(a, key), textb = begfield(b, key);
			else {
				texta = a->text, textb = b->text;
				if (key->skipsblanks) {
					while (texta < lima && blanks[UCHAR(*texta)]) ++texta;
					while (textb < limb && blanks[UCHAR(*textb)]) ++textb;
				}
			}
		} else {
			/*
			 * For the first iteration only, the key positions
			 * have been precomputed for us.
			 */
			texta = a->keybeg, lima = a->keylim;
			textb = b->keybeg, limb = b->keylim;
		}

		/* Find the lengths. */
		lena = lima - texta, lenb = limb - textb;
		if (lena < 0) lena = 0;
		if (lenb < 0) lenb = 0;

		if (key->skipeblanks) {
			char *a_end = texta + lena;
			char *b_end = textb + lenb;
			trim_trailing_blanks(texta, &a_end);
			trim_trailing_blanks(textb, &b_end);
			lena = a_end - texta;
			lenb = b_end - textb;
		}
		/* Actually compare the fields. */
		if (key->mixed_numeric) {
			if (*lima || *limb) {
				char savea = *lima, saveb = *limb;
				*lima = *limb = '\0';
				diff = mixed_numcompare(texta, textb);
				*lima = savea, *limb = saveb;
			} else diff = mixed_numcompare(texta, textb);
			if (diff) return key->reverse ? -diff : diff;
			continue;
		} else if (key->numeric) {
			if (*lima || *limb) {
				char savea = *lima, saveb = *limb;
				*lima = *limb = '\0';
				diff = numcompare(texta, textb);
				*lima = savea, *limb = saveb;
			} else diff = numcompare(texta, textb);
			if (diff) return key->reverse ? -diff : diff;
			continue;
		} else if (key->general_numeric) {
			if (*lima || *limb) {
				char savea = *lima, saveb = *limb;
				*lima = *limb = '\0';
				diff = general_numcompare(texta, textb);
				*lima = savea, *limb = saveb;
			} else diff = general_numcompare(texta, textb);
			if (diff) return key->reverse ? -diff : diff;
			continue;
		} else if (key->month) {
			diff = getmonth(texta, lena) - getmonth(textb, lenb);
			if (diff) return key->reverse ? -diff : diff;
			continue;
		} else if (ignore && translate)
#define CMP_WITH_IGNORE(A, B) \
	do { \
	  while (texta < lima && textb < limb) { \
	      while (texta < lima && ignore[UCHAR (*texta)]) ++texta; \
	      while (textb < limb && ignore[UCHAR (*textb)]) ++textb; \
	      if (texta < lima && textb < limb)	{ \
		  if ((A) != (B)) { \
		      diff = (A) - (B); \
		      break; \
		  } \
		  ++texta; ++textb; \
		} \
	    if (texta == lima && textb < limb && !ignore[UCHAR (*textb)]) diff = -1; \
	    else if (texta < lima && textb == limb && !ignore[UCHAR (*texta)]) diff = 1; \
	  } \
	  \
	  if (diff == 0) { \
	      while (texta < lima && ignore[UCHAR (*texta)]) ++texta; \
	      while (textb < limb && ignore[UCHAR (*textb)]) ++textb; \
		  \
	      if (texta == lima && textb < limb) diff = -1; \
	      else if (texta < lima && textb == limb) diff = 1; \
	  } \
	  /* Relative lengths are meaningless if characters were ignored. \
	     Handling this case here avoids what might be an invalid length \
	     comparison below.  */ \
	  if (diff == 0 && texta == lima && textb == limb) return 0; \
    } while (0)

			CMP_WITH_IGNORE(translate[UCHAR(*texta)], translate[UCHAR(*textb)]);
		else if (ignore)
			CMP_WITH_IGNORE(UCHAR(*texta), UCHAR(*textb));
		else if (translate)
			while (texta < lima && textb < limb) {
				if (translate[UCHAR(*texta++)] != translate[UCHAR(*textb++)]) {
					diff = (translate[UCHAR(*--texta)] - translate[UCHAR(*--textb)]);
					break;
				}
			}
		else diff = memcmp(texta, textb, min(lena, lenb));
		if (diff) return key->reverse ? -diff : diff;
		if ((diff = lena - lenb) != 0) return key->reverse ? -diff : diff;
	}
	return 0;
}

/*
 * Compare two lines A and B, returning negative, zero, or positive depending
 * on whether A compares less than, equal to, or greater than B.
 */

static inline int compare(register const struct line * a, register const struct line * b)
{
	int             diff, tmpa, tmpb, mini;

	/*
	 * First try to compare on the specified keys (if any). The only two
	 * cases with no key at all are unadorned sort, and unadorned sort
	 * -r.
	 */
	if (keyhead.next) {
		diff = keycompare(a, b);
		if (diff != 0) return diff;
		if (unique || stable) return 0;
	}
	/*
	 * If the keys all compare equal (or no keys were specified) fall
	 * through to the default byte-by-byte comparison.
	 */
	tmpa = a->length, tmpb = b->length;
	mini = min(tmpa, tmpb);
	if (mini == 0) diff = tmpa - tmpb;
	else {
		char *ap = a->text, *bp = b->text;
		diff = UCHAR(*ap) - UCHAR(*bp);
		if (diff == 0) {
			diff = memcmp(ap, bp, mini);
			if (diff == 0) diff = tmpa - tmpb;
		}
	}
	return reverse ? -diff : diff;
}

/*
 * Check that the lines read from the given FP come in order.  Return 1 if
 * they do and 0 if there is a disorder. FIXME: return number of first
 * out-of-order line if not sorted.
 */

static int checkfp(FILE * fp)
{
	struct buffer   buf;	/* Input buffer. */
	struct lines    lines;	/* Lines scanned from the buffer. */
	struct line     temp;	/* Copy of previous line. */
	int             cc;	/* Character count. */
	int             alloc, sorted = 1;

	initbuf(&buf, mergealloc);
	initlines(&lines, mergealloc / linelength + 1, LINEALLOC / ((NMERGE + NMERGE) * sizeof(struct line)));
	alloc = linelength;
	temp.text = xmalloc(alloc);

	cc = fillbuf(&buf, fp);
	if (cc == 0) goto finish;

	findlines(&buf, &lines);

	while (1) {
		struct line    *prev_line;	/* Pointer to previous line. */
		int             cmp;	/* Result of calling compare. */
		int             i;

		/* Compare each line in the buffer with its successor. */
		for (i = 0; i < lines.used - 1; ++i) {
			cmp = compare(&lines.lines[i], &lines.lines[i + 1]);
			if ((unique && cmp >= 0) || (cmp > 0)) {
				sorted = 0;
				goto finish;
			}
		}

		/* Save the last line of the buffer and refill the buffer. */
		prev_line = lines.lines + (lines.used - 1);
		if (prev_line->length + 1 > alloc) {
			do {
				alloc *= 2;
			}
			while (alloc < prev_line->length + 1);
			temp.text = xrealloc(temp.text, alloc);
		}
		assert(prev_line->length + 1 <= alloc);
		memcpy(temp.text, prev_line->text, prev_line->length + 1);
		temp.length = prev_line->length;
		temp.keybeg = temp.text + (prev_line->keybeg - prev_line->text);
		temp.keylim = temp.text + (prev_line->keylim - prev_line->text);

		cc = fillbuf(&buf, fp);
		if (cc == 0) break;

		findlines(&buf, &lines);
		/*
		 * Make sure the line saved from the old buffer contents is
		 * less than or equal to the first line of the new buffer.
		 */
		cmp = compare(&temp, &lines.lines[0]);
		if ((unique && cmp >= 0) || (cmp > 0)) {
			sorted = 0;
			break;
		}
	}

finish:
	xfclose(fp);
	free(buf.buf);
	free((char *) lines.lines);
	free(temp.text);
	return sorted;
}

/* lh3: for adjusting a heap structure */
static void algo_heap_adjust(struct line *ls[], int l[], int i, int n)
{
	int k, tmp = l[i];
	for (;;) {
		k = (i << 1) + 1; // the left child
		if (k >= n) { l[i] = tmp; return; }
		if (k < n - 1 && compare(ls[l[k]], ls[l[k+1]]) > 0) ++k;
		if (compare(ls[tmp], ls[l[k]]) > 0) { l[i] = l[k]; i = k; }
		else { l[i] = tmp; return; }
	}
}

/*
 * Merge lines from FPS onto OFP.  NFPS cannot be greater than NMERGE. Close
 * FPS before returning.
 */
/* lh3: Marginally improve the speed by using a heap structure */
static void mergefps(FILE ** fps, register int nfps, FILE * ofp)
{
	struct buffer   buffer[NMERGE];	/* Input buffers for each file. */
	struct lines    lines[NMERGE];	/* Line tables for each buffer. */
	struct line     saved;	/* Saved line for unique check. */
	struct line    *ls[NMERGE];
	int             savedflag = 0;	/* True if there is a saved line. */
	int             savealloc = 0;	/* Size allocated for the saved line. */
	int             cur[NMERGE];	/* Current line in each line table. */
	int             ord[NMERGE];	/* Table representing a permutation
					 * of fps, such that
					 * lines[ord[0]].lines[cur[ord[0]]] is the smallest
					 * line and will be next output. */
	register int    i, j;

	/* Allocate space for a saved line if necessary. */
	if (unique) {
		savealloc = linelength;
		saved.text = xmalloc(savealloc);
	}
	/* Read initial lines from each input file. */
	for (i = 0; i < nfps; ++i) {
		initbuf(&buffer[i], mergealloc);
		/* If a file is empty, eliminate it from future consideration. */
		while (i < nfps && !fillbuf(&buffer[i], fps[i])) {
			xfclose(fps[i]);
			--nfps;
			for (j = i; j < nfps; ++j) fps[j] = fps[j + 1];
		}
		if (i == nfps) free(buffer[i].buf);
		else {
			initlines(&lines[i], mergealloc / linelength + 1,
				  LINEALLOC / ((NMERGE + NMERGE) * sizeof(struct line)));
			findlines(&buffer[i], &lines[i]);
			cur[i] = 0;
		}
	}

	/* make a heap structure */
	for (i = 0; i < nfps; ++i) {
		ord[i] = i;
		ls[i] = &lines[i].lines[cur[i]];
	}
	for (i = (nfps >> 1) - 1; i >= 0; --i) 
		algo_heap_adjust(ls, ord, i, nfps);

	/* Repeatedly output the smallest line until no input remains. */
	while (nfps) {
		/*
		 * If uniqified output is turned on, output only the first of
		 * an identical series of lines.
		 */
		if (unique) {
			if (savedflag && compare(&saved, &lines[ord[0]].lines[cur[ord[0]]])) {
				write_bytes(saved.text, saved.length, ofp);
				putc(eolchar, ofp);
				savedflag = 0;
			}
			if (!savedflag) {
				if (savealloc < lines[ord[0]].lines[cur[ord[0]]].length + 1) {
					while (savealloc < lines[ord[0]].lines[cur[ord[0]]].length + 1)
						savealloc *= 2;
					saved.text = xrealloc(saved.text, savealloc);
				}
				saved.length = lines[ord[0]].lines[cur[ord[0]]].length;
				memcpy(saved.text, lines[ord[0]].lines[cur[ord[0]]].text,
				       saved.length + 1);
				if (lines[ord[0]].lines[cur[ord[0]]].keybeg != NULL) {
					saved.keybeg = saved.text +
						(lines[ord[0]].lines[cur[ord[0]]].keybeg
						 - lines[ord[0]].lines[cur[ord[0]]].text);
				}
				if (lines[ord[0]].lines[cur[ord[0]]].keylim != NULL) {
					saved.keylim = saved.text +
						(lines[ord[0]].lines[cur[ord[0]]].keylim
						 - lines[ord[0]].lines[cur[ord[0]]].text);
				}
				savedflag = 1;
			}
		} else {
			write_bytes(lines[ord[0]].lines[cur[ord[0]]].text,
			      lines[ord[0]].lines[cur[ord[0]]].length, ofp);
			putc(eolchar, ofp);
		}

		/* Check if we need to read more lines into core. */
		if (++cur[ord[0]] == lines[ord[0]].used) {
			if (fillbuf(&buffer[ord[0]], fps[ord[0]])) {
				findlines(&buffer[ord[0]], &lines[ord[0]]);
				cur[ord[0]] = 0;
			} else {
				/* We reached EOF on fps[ord[0]]. */
				for (i = 1; i < nfps; ++i)
					if (ord[i] > ord[0]) --ord[i];
				--nfps;
				xfclose(fps[ord[0]]);
				free(buffer[ord[0]].buf);
				free((char *) lines[ord[0]].lines);
				for (i = ord[0]; i < nfps; ++i) {
					fps[i] = fps[i + 1];
					buffer[i] = buffer[i + 1];
					lines[i] = lines[i + 1];
					cur[i] = cur[i + 1];
					ls[i] = ls[i + 1];
				}
				ord[0] = ord[nfps];
			}
		}
		ls[ord[0]] = &lines[ord[0]].lines[cur[ord[0]]];
		algo_heap_adjust(ls, ord, 0, nfps);
	}

	if (unique && savedflag) {
		write_bytes(saved.text, saved.length, ofp);
		putc(eolchar, ofp);
		free(saved.text);
	}
}

/* Sort the array LINES with NLINES members, using TEMP for temporary space. */

/* lh3: The following codes implement a bottom-up iterative
 * MergeSort. The original codes implemented a top-down recursive
 * MergeSort which is quite inefficient given large dataset. Note that
 * for sort in textutils-1.22 or before, sorting was not the speed
 * bottleneck as sort at the old age would not put too many lines in the
 * main memory. However, as the current sort is trying to use all
 * available memory, recursive codes will greatly harm the
 * performance. Actually, I would recommend the latest sort to use my
 * codes instead. */

static void sortlines(struct line *array, int nlines, struct line *temp)
{
	struct line *a2[2], *a, *b;
	int curr, shift;
	size_t i, j, k, n;

	n = nlines;
	a2[0] = array; a2[1] = temp;
	for (curr = 0, shift = 0; (1ul<<shift) < n; ++shift) {
		a = a2[curr]; b = a2[1-curr];
		if (shift == 0) {
			for (i = 0; i < n; i += 2) {
				if (i == n - 1) {
					b[i] = a[i];
					break;
				}
				if (compare(&a[i+1], &a[i]) < 0) {
					b[i] = a[i+1]; b[i+1] = a[i];
				} else {
					b[i] = a[i]; b[i+1] = a[i+1];
				}
			}
		} else {
			size_t step = 1ul<<shift;
			for (i = 0; i < n; i += step<<1) {
				size_t ea, eb;
				struct line *p;
				if (n < i + step) {
					ea = n; eb = 0;
				} else {
					ea = i + step;
					eb = n < i + (step<<1)? n : i + (step<<1);
				}
				j = i; k = i + step; p = b + i;
				while (j < ea && k < eb) {
					if (compare(&a[j], &a[k]) < 0) *p++ = a[j++];
					else *p++ = a[k++];
				}
				while (j < ea) *p++ = a[j++];
				while (k < eb) *p++ = a[k++];
			}
		}
		curr = 1 - curr;
	}
	if (curr == 1)
		for (i = 0; i < n; ++i) array[i] = temp[i];
}

/*
 * Check that each of the NFILES FILES is ordered. Return a count of
 * disordered files.
 */

static int check(char **files, int nfiles)
{
	int             i, disorders = 0;
	FILE           *fp;

	for (i = 0; i < nfiles; ++i) {
		fp = xfopen(files[i], "r");
		if (!checkfp(fp)) {
			fprintf(stderr, _("%s: disorder on %s\n"), program_name, files[i]);
			++disorders;
		}
	}
	return disorders;
}

/* Merge NFILES FILES onto OFP. */

static void merge(char **files, int nfiles, FILE * ofp)
{
	int             i, j, t;
	char           *temp;
	FILE           *fps[NMERGE], *tfp;

	while (nfiles > NMERGE) {
		t = 0;
		for (i = 0; i < nfiles / NMERGE; ++i) {
			for (j = 0; j < NMERGE; ++j)
				fps[j] = xfopen(files[i * NMERGE + j], "r");
			tfp = xtmpfopen(temp = tempname());
			mergefps(fps, NMERGE, tfp);
			xfclose(tfp);
			for (j = 0; j < NMERGE; ++j)
				zaptemp(files[i * NMERGE + j]);
			files[t++] = temp;
		}
		for (j = 0; j < nfiles % NMERGE; ++j)
			fps[j] = xfopen(files[i * NMERGE + j], "r");
		tfp = xtmpfopen(temp = tempname());
		mergefps(fps, nfiles % NMERGE, tfp);
		xfclose(tfp);
		for (j = 0; j < nfiles % NMERGE; ++j)
			zaptemp(files[i * NMERGE + j]);
		files[t++] = temp;
		nfiles = t;
	}
	for (i = 0; i < nfiles; ++i) fps[i] = xfopen(files[i], "r");
	mergefps(fps, i, ofp);
	for (i = 0; i < nfiles; ++i)
		zaptemp(files[i]);
}

/* Sort NFILES FILES onto OFP. */

static void sort(char **files, int nfiles, FILE *ofp)
{
	struct buffer   buf;
	struct lines    lines;
	struct line    *tmp;
	int             i, ntmp;
	FILE           *fp, *tfp;
	struct tempnode *node;
	int             n_temp_files = 0;
	char          **tempfiles;

	initbuf(&buf, sortalloc);
	initlines(&lines, sortalloc / linelength + 1, LINEALLOC / sizeof(struct line));
	ntmp = lines.alloc;
	tmp = (struct line *) xmalloc(ntmp * sizeof(struct line));

	while (nfiles--) {
		fp = xfopen(*files++, "r");
		while (fillbuf(&buf, fp)) {
			findlines(&buf, &lines);
			if (lines.used > ntmp) {
				while (lines.used > ntmp) ntmp *= 2;
				tmp = (struct line *)xrealloc((char *) tmp, ntmp * sizeof(struct line));
			}
			sortlines(lines.lines, lines.used, tmp);
			if (feof(fp) && !nfiles && !n_temp_files && !buf.left) tfp = ofp;
			else {
				++n_temp_files;
				tfp = xtmpfopen(tempname());
			}
			for (i = 0; i < lines.used; ++i)
				if (!unique || i == 0 || compare(&lines.lines[i], &lines.lines[i - 1])) {
					write_bytes(lines.lines[i].text, lines.lines[i].length, tfp);
					putc(eolchar, tfp);
				}
			if (tfp != ofp) xfclose(tfp);
		}
		xfclose(fp);
	}

	free(buf.buf);
	free((char *) lines.lines);
	free((char *) tmp);

	if (n_temp_files) {
		tempfiles = (char **) xmalloc(n_temp_files * sizeof(char *));
		i = n_temp_files;
		for (node = temphead.next; i > 0; node = node->next)
			tempfiles[--i] = node->name;
		merge(tempfiles, n_temp_files, ofp);
		free((char *) tempfiles);
	}
}

/* Insert key KEY at the end of the list (`keyhead'). */

static void insertkey(struct keyfield * key)
{
	struct keyfield *k = &keyhead;
	while (k->next) k = k->next;
	k->next = key;
	key->next = NULL;
}

static void badfieldspec(const char *s)
{
	error(SORT_FAILURE, 0, _("invalid field specification `%s'"), s);
}

/* Handle interrupts and hangups. */

static void sighandler(int sig)
{
#ifdef SA_INTERRUPT
	struct sigaction sigact;

	sigact.sa_handler = SIG_DFL;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(sig, &sigact, NULL);
#else /* !SA_INTERRUPT */
	signal(sig, SIG_DFL);
#endif /* SA_INTERRUPT */
	cleanup();
	kill(getpid(), sig);
}

/*
 * Set the ordering options for KEY specified in S. Return the address of the
 * first character in S that is not a valid ordering option. BLANKTYPE is the
 * kind of blanks that 'b' should skip.
 */

static char *set_ordering(register const char *s, struct keyfield * key,
						  enum blanktype blanktype)
{
	while (*s) {
		switch (*s) {
		case 'b':
			if (blanktype == bl_start || blanktype == bl_both)
				key->skipsblanks = 1;
			if (blanktype == bl_end || blanktype == bl_both)
				key->skipeblanks = 1;
			break;
		case 'd': key->ignore = nondictionary; break;
		case 'f': key->translate = fold_toupper; break;
		case 'g': key->general_numeric = 1; break;
		case 'N': key->mixed_numeric = 1; break;
		case 'i': key->ignore = nonprinting; break;
		case 'M': key->month = 1; break;
		case 'n': key->numeric = 1; break;
		case 'r': key->reverse = 1; break;
		default: return (char *) s;
		}
		++s;
	}
	return (char *) s;
}

static void key_init(struct keyfield * key)
{
	memset(key, 0, sizeof(*key));
	key->eword = -1;
}

int main(int argc, char **argv)
{
	struct keyfield *key = NULL, gkey;
	char           *s;
	int             i, t, t2;
	int             checkonly = 0, mergeonly = 0, nfiles = 0;
	char           *minus = "-", *outfile = minus, **files, *tmp;
	FILE           *ofp;
#ifdef SA_INTERRUPT
	struct sigaction oldact, newact;
#endif /* SA_INTERRUPT */

	program_name = argv[0];
	setlocale(LC_ALL, "");
	bindtextdomain(PACKAGE, LOCALEDIR);
	textdomain(PACKAGE);

	have_read_stdin = 0;
	inittables();

	temp_file_prefix = getenv("TMPDIR");
	if (temp_file_prefix == NULL)
		temp_file_prefix = DEFAULT_TMPDIR;

#ifdef SA_INTERRUPT
	newact.sa_handler = sighandler;
	sigemptyset(&newact.sa_mask);
	newact.sa_flags = 0;

	sigaction(SIGINT, NULL, &oldact);
	if (oldact.sa_handler != SIG_IGN)
		sigaction(SIGINT, &newact, NULL);
	sigaction(SIGHUP, NULL, &oldact);
	if (oldact.sa_handler != SIG_IGN)
		sigaction(SIGHUP, &newact, NULL);
	sigaction(SIGPIPE, NULL, &oldact);
	if (oldact.sa_handler != SIG_IGN)
		sigaction(SIGPIPE, &newact, NULL);
	sigaction(SIGTERM, NULL, &oldact);
	if (oldact.sa_handler != SIG_IGN)
		sigaction(SIGTERM, &newact, NULL);
#else				/* !SA_INTERRUPT */
	if (signal(SIGINT, SIG_IGN) != SIG_IGN)
		signal(SIGINT, sighandler);
	if (signal(SIGHUP, SIG_IGN) != SIG_IGN)
		signal(SIGHUP, sighandler);
	if (signal(SIGPIPE, SIG_IGN) != SIG_IGN)
		signal(SIGPIPE, sighandler);
	if (signal(SIGTERM, SIG_IGN) != SIG_IGN)
		signal(SIGTERM, sighandler);
#endif				/* !SA_INTERRUPT */

	gkey.sword = gkey.eword = -1;
	gkey.ignore = NULL;
	gkey.translate = NULL;
	gkey.numeric = gkey.mixed_numeric = gkey.general_numeric = gkey.month = gkey.reverse = 0;
	gkey.skipsblanks = gkey.skipeblanks = 0;

	files = (char **) xmalloc(sizeof(char *) * argc);

	for (i = 1; i < argc; ++i) {
		if (argv[i][0] == '+') {
			if (key) insertkey(key);
			key = (struct keyfield *) xmalloc(sizeof(struct keyfield));
			key_init(key);
			s = argv[i] + 1;
			if (!(ISDIGIT(*s) || (*s == '.' && ISDIGIT(s[1]))))
				badfieldspec(argv[i]);
			for (t = 0; ISDIGIT(*s); ++s) t = 10 * t + *s - '0';
			t2 = 0;
			if (*s == '.')
				for (++s; ISDIGIT(*s); ++s)
					t2 = 10 * t2 + *s - '0';
			if (t2 || t) {
				key->sword = t;
				key->schar = t2;
			} else
				key->sword = -1;
			s = set_ordering(s, key, bl_start);
			if (*s) badfieldspec(argv[i]);
		} else if (argv[i][0] == '-' && argv[i][1]) { /* parse command line options */
			s = argv[i] + 1;
			if (ISDIGIT(*s) || (*s == '.' && ISDIGIT(s[1]))) {
				if (!key) {
					/* Provoke with `sort -9'.  */
					error(0, 0, _("when using the old-style +POS and -POS key specifiers,\nthe +POS specifier must come first"));
					usage(SORT_FAILURE);
				}
				for (t = 0; ISDIGIT(*s); ++s) t = t * 10 + *s - '0';
				t2 = 0;
				if (*s == '.')
					for (++s; ISDIGIT(*s); ++s)
						t2 = t2 * 10 + *s - '0';
				key->eword = t;
				key->echar = t2;
				s = set_ordering(s, key, bl_end);
				if (*s) badfieldspec(argv[i]);
				insertkey(key);
				key = NULL;
			} else {
				while (*s) {
					s = set_ordering(s, &gkey, bl_both);
					switch (*s) {
					case '\0': break;
					case 'c': checkonly = 1; break;
					case 'h': usage(0); break;
					case 'k':
						if (s[1]) ++s;
						else {
							if (i == argc - 1) error(SORT_FAILURE, 0, _("option `-k' requires an argument"));
							else s = argv[++i];
						}
						if (key) insertkey(key);
						key = (struct keyfield *)xmalloc(sizeof(struct keyfield));
						key_init(key);
						/* Get POS1. */
						if (!ISDIGIT(*s)) badfieldspec(argv[i]);
						for (t = 0; ISDIGIT(*s); ++s) t = 10 * t + *s - '0';
						if (t == 0) { /* Provoke with `sort -k0' */
							error(0, 0, _("the starting field number argument to the `-k' option must be positive"));
							badfieldspec(argv[i]);
						}
						--t;
						t2 = 0;
						if (*s == '.') {
							if (!ISDIGIT(s[1])) { /* Provoke with `sort -k1.' */
								error(0, 0, _("starting field spec has `.' but lacks following character offset"));
								badfieldspec(argv[i]);
							}
							for (++s; ISDIGIT(*s); ++s) t2 = 10 * t2 + *s - '0';
							if (t2 == 0) { /* Provoke with `sort -k1.0' */
								error(0, 0, _("starting field character offset argument to the `-k' option\nmust be positive"));
								badfieldspec(argv[i]);
							}
							--t2;
						}
						if (t2 || t) {
							key->sword = t;
							key->schar = t2;
						} else key->sword = -1;
						s = set_ordering(s, key, bl_start);
						if (*s == 0) {
							key->eword = -1;
							key->echar = 0;
						} else if (*s != ',')
							badfieldspec(argv[i]);
						else if (*s == ',') {
							/* Skip over comma.  */
							++s;
							if (*s == 0) { /* Provoke with `sort -k1,' */
								error(0, 0, _("field specification has `,' but lacks following field spec"));
								badfieldspec(argv[i]);
							}
							/* Get POS2. */
							for (t = 0; ISDIGIT(*s); ++s)
								t = t * 10 + *s - '0';
							if (t == 0) { /* Provoke with `sort -k1,0' */
								error(0, 0, _("ending field number argument to the `-k' option must be positive"));
								badfieldspec(argv[i]);
							}
							--t;
							t2 = 0;
							if (*s == '.') {
								if (!ISDIGIT(s[1])) { /* Pro vok e wit h `so rt -k1,1.' */
									error(0, 0, _("ending field spec has `.' but lacks following character offset"));
									badfieldspec(argv[i]);
								}
								for (++s; ISDIGIT(*s); ++s)
									t2 = t2 * 10 + *s - '0';
							} else ++t; /* `-k 2,3' is equivalent to `+1 -3'. */
							key->eword = t;
							key->echar = t2;
							s = set_ordering(s, key, bl_end);
							if (*s)
								badfieldspec(argv[i]);
						}
						insertkey(key);
						key = NULL;
						goto outer;
					case 'm': mergeonly = 1; break;
					case 'o':
						if (s[1]) outfile = s + 1;
						else {
							if (i == argc - 1) error(SORT_FAILURE, 0, _("option `-o' requires an argument"));
							else outfile = argv[++i];
						}
						goto outer;
					case 's': stable = 1; break;
					case 't':
						if (s[1]) tab = *++s;
						else if (i < argc - 1) {
							tab = *argv[++i];
							goto outer;
						} else error(SORT_FAILURE, 0, _("option `-t' requires an argument"));
						break;
					case 'T':
						if (s[1]) temp_file_prefix = ++s;
						else {
							if (i < argc - 1) temp_file_prefix = argv[++i];
							else error(SORT_FAILURE, 0, _("option `-T' requires an argument"));
						}
						goto outer;
					case 'u': unique = 1; break;
					case 'z': eolchar = 0; break;
					case 'y': goto outer; /* Accept and ignore e.g. -y0 for compatibility with Solaris 2. */
					default:
						fprintf(stderr, _("%s: unrecognized option `-%c'\n"), argv[0], *s);
						usage(SORT_FAILURE);
					}
					if (*s) ++s;
				}
			}
		} else {	/* Not an option. */
			files[nfiles++] = argv[i];
		}
outer:	;
	}

	if (key) insertkey(key);

	/* Inheritance of global options to individual keys. */
	for (key = keyhead.next; key; key = key->next) {
		if (!key->ignore && !key->translate && !key->skipsblanks && !key->reverse
		    && !key->skipeblanks && !key->month && !key->numeric
		    && !key->general_numeric && !key->mixed_numeric) {
			key->ignore = gkey.ignore;
			key->translate = gkey.translate;
			key->skipsblanks = gkey.skipsblanks;
			key->skipeblanks = gkey.skipeblanks;
			key->month = gkey.month;
			key->numeric = gkey.numeric;
			key->general_numeric = gkey.general_numeric;
			key->mixed_numeric = gkey.mixed_numeric;
			key->reverse = gkey.reverse;
		}
	}
	if (!keyhead.next && (gkey.ignore || gkey.translate || gkey.skipsblanks
						  || gkey.skipeblanks || gkey.month || gkey.numeric
						  || gkey.general_numeric || gkey.mixed_numeric))
	{
		insertkey(&gkey);
	}
	reverse = gkey.reverse;

	if (nfiles == 0) {
		nfiles = 1;
		files = &minus;
	}
	if (checkonly) {
		/*
		 * POSIX requires that sort return 1 IFF invoked with -c and the
		 * input is not properly sorted.
		 */
		exit(check(files, nfiles) == 0 ? 0 : 1);
	}
	if (strcmp(outfile, "-")) {
		struct stat     outstat;
		if (stat(outfile, &outstat) == 0) {
			/*
			 * The following code prevents a race condition when people
			 * use the brain dead shell programming idiom: cat file |
			 * sort -o file This feature is provided for historical
			 * compatibility, but we strongly discourage ever relying on
			 * this in new shell programs.
			 */

			/*
			 * Temporarily copy each input file that might be another
			 * name for the output file.  When in doubt (e.g. a pipe),
			 * copy.
			 */
			for (i = 0; i < nfiles; ++i) {
				char            buf[8192];
				FILE           *fp;
				int             cc;

				if (S_ISREG(outstat.st_mode) && strcmp(outfile, files[i])) {
					struct stat     instat;
					if ((strcmp(files[i], "-")
					     ? stat(files[i], &instat)
					     : fstat(STDIN_FILENO, &instat)) != 0)
					{
						error(0, errno, "%s", files[i]);
						cleanup();
						exit(SORT_FAILURE);
					}
					if (S_ISREG(instat.st_mode)
						&& (instat.st_ino != outstat.st_ino
							|| instat.st_dev != outstat.st_dev))
					{
						/* We know the files are distinct. */
						continue;
					}
				}
				fp = xfopen(files[i], "r");
				tmp = tempname();
				ofp = xtmpfopen(tmp);
				while ((cc = fread(buf, 1, sizeof buf, fp)) > 0)
					write_bytes(buf, cc, ofp);
				if (ferror(fp)) {
					error(0, errno, "%s", files[i]);
					cleanup();
					exit(SORT_FAILURE);
				}
				xfclose(ofp);
				xfclose(fp);
				files[i] = tmp;
			}
		}
		ofp = xfopen(outfile, "w");
	} else ofp = stdout;

	if (mergeonly) merge(files, nfiles, ofp);
	else sort(files, nfiles, ofp);
	cleanup();

	/*
	 * If we wait for the implicit flush on exit, and the parent process
	 * has closed stdout (e.g., exec >&- in a shell), then the output
	 * file winds up empty.  I don't understand why.  This is under
	 * SunOS, Solaris, Ultrix, and Irix.  This premature fflush makes
	 * the output reappear. --karl@cs.umb.edu
	 */
	if (fflush(ofp) < 0) error(SORT_FAILURE, errno, _("%s: write error"), outfile);

	if (have_read_stdin && fclose(stdin) == EOF) error(SORT_FAILURE, errno, outfile);
	if (ferror(stdout) || fclose(stdout) == EOF)
		error(SORT_FAILURE, errno, _("%s: write error"), outfile);

	exit(EXIT_SUCCESS);
}
