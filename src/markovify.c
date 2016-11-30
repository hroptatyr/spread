#if defined HAVE_CONFIG_H
# include "config.h"
#endif	/* HAVE_CONFIG_H */
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#if defined HAVE_DFP754_H
# include <dfp754.h>
#endif	/* HAVE_DFP754_H */
#include "nifty.h"
#include "dfp754_d32.h"
#include "dfp754_d64.h"

typedef _Decimal32 px_t;
typedef _Decimal64 qx_t;
#define strtopx		strtod32
#define pxtostr		d32tostr
#define strtoqx		strtod64
#define qxtostr		d64tostr

/* hashed prices, most significant digit plus exponent, [0,31) */
typedef unsigned int hx_t;

static size_t arity = 1U;


static __attribute__((format(printf, 1, 2))) void
serror(const char *fmt, ...)
{
	va_list vap;
	va_start(vap, fmt);
	vfprintf(stderr, fmt, vap);
	va_end(vap);
	if (errno) {
		fputc(':', stderr);
		fputc(' ', stderr);
		fputs(strerror(errno), stderr);
	}
	fputc('\n', stderr);
	return;
}

static __attribute__((const, pure)) hx_t
hashpx(px_t p)
{
/* use most significant (base10) digit and exponent, alongside quantum */
	int e = 0;
	int q = quantexpd32(p);
	int v = (int)scalbnd32(frexpd32(p, &e), 1);
	const int f = (e - q - 1);
	int x = 10 * f + v - f;
	return ~(hx_t)((v > 0 && x < 32) - 1U) & (hx_t)x;
}


/* limit to arity of 8 (which would consume 8TB of RAM) */
static hx_t hist[8U];
static size_t nhist;
/* big counting array */
static size_t *chain;

static int
push_beef(char *ln, size_t UNUSED(lz))
{
	char *on;
	px_t s = strtopx(ln, &on);

	if (UNLIKELY(*on != '\n')) {
		/* ignore */
		return -1;
	}
	/* hash him */
	hist[nhist] = hashpx(s);
	nhist = (nhist + 1U) % arity;
	return 0;
}

static int
addchain(void)
{
	size_t k = 0U;

	for (size_t i = nhist; i > 0U;) {
		k *= 32U;
		k += hist[--i];
	}
	for (size_t i = arity; i > nhist;) {
		k *= 32U;
		k += hist[--i];
	}
	/* actually increment the count now */
	chain[k]++;
	return 0;
}

static int
prchain(void)
{
	for (size_t i = 0U, n = 1ULL << (5U * arity); i < n; i++) {
		if (LIKELY(!chain[i])) {
			continue;
		}
		for (size_t j = 0, k = 1ULL; j < arity; j++, k <<= 5U) {
			printf("%zu\t", (i / k) % 32U);
		}
		printf("%zu\n", chain[i]);
	}
	return 0;
}


#include "markovify.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	}

	if (argi->arity_arg) {
		if (UNLIKELY(!(arity = strtoul(argi->arity_arg, NULL, 10)))) {
			errno = 0, serror("\
Error: arity must be strictly positive");
			rc = 1;
			goto out;
		} else if (arity > countof(hist)) {
			errno = 0, serror("\
Error: arity value too big\n\
Info: arity=8 needs 8TB RAM, arity=9 would need 256TB RAM");
			rc = 1;
			goto out;
		}
	}

	/* initialise counting array */
	with (size_t z = 1ULL << (5U * arity)) {
		chain = calloc(z, sizeof(*chain));

		if (UNLIKELY(chain == NULL)) {
			rc = 1;
			goto out;
		}
	}

	{
		char *line = NULL;
		size_t llen = 0UL;
		ssize_t nrd;

		while ((nrd = getline(&line, &llen, stdin)) > 0) {
			if (UNLIKELY(push_beef(line, nrd) < 0)) {
				;
			} else {
				addchain();
			}
		}

		free(line);

		/* print results */
		prchain();
	}

	free(chain);

out:
	yuck_free(argi);
	return rc;
}

/* markovify.c ends here */
