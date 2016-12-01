#if defined HAVE_CONFIG_H
# include "config.h"
#endif	/* HAVE_CONFIG_H */
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
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

/* chain count */
typedef size_t cc_t;

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
pxtohx(px_t p, int q)
{
/* use most significant (base10) digit and exponent, alongside quantum */
	int e = 0;
	int v = (int)scalbnd32(frexpd32(p, &e), 1);
	const int f = (e - q - 1);
	int x = 10 * f + v - f;
	return ~(hx_t)((v > 0 && x < 32) - 1U) & (hx_t)x;
}

static __attribute__((const, pure)) px_t
hxtopx(hx_t h, int q)
{
/* try and reconstruct the upper bound representant that would hash to H */
	const int e = (h - 1) / 9U;
	const int v = (h + e) % 10U;
	const px_t x = scalbnd32(0.df, q);
	return h ? quantized32(ldexpd32((px_t)v, q + e), x) : nand32("");
}

static ssize_t
cctostr(char *restrict buf, size_t bsz, cc_t c)
{
	return snprintf(buf, bsz, "%zu", c);
}


/* limit to arity of 8 (which would consume 8TB of RAM) */
static hx_t hist[8U];
static size_t nhist;
/* big counting array */
static cc_t *chain;
/* dataset quantum */
static int qunt;

static int
push_beef(char *ln, size_t UNUSED(lz))
{
	char *on;
	px_t s = strtopx(ln, &on);
	int q;

	if (UNLIKELY(*on != '\n')) {
		/* ignore */
		return -1;
	} else if (UNLIKELY(qunt != (q = quantexpd32(s)))) {
		if (LIKELY(qunt)) {
			errno = 0, serror("\
Warning: quantum has changed");
			return -1;
		}
		qunt = q;
	}
	/* hash him */
	hist[nhist] = pxtohx(s, qunt);
	nhist = (nhist + 1U) % arity;
	return 0;
}

static int
add_chain(void)
{
	size_t k = 0U;

	for (size_t i = nhist; i < arity; i++) {
		k *= 32U;
		k += hist[i];
	}
	for (size_t i = 0U; i < nhist; i++) {
		k *= 32U;
		k += hist[i];
	}
	/* actually increment the count now */
	chain[k]++;
	return 0;
}

static int
prnt_chain(bool rawp)
{
	ssize_t(*ixtostr)(char*restrict, size_t, size_t idx);

	static ssize_t ixtopstr(char *restrict buf, size_t bsz, size_t i)
	{
		size_t len = 0U;

		for (size_t j = arity; --j;) {
			const hx_t h = (i >> (5U * j)) % 32U;
			const px_t p = hxtopx(h, qunt);

			len += pxtostr(buf + len, bsz - len, p);
			buf[len++] = '\t';
		}
		/* last one is special */
		if (i % 32U) {
			const hx_t h = i % 32U;
			const px_t p = hxtopx(h, qunt);

			len += pxtostr(buf + len, bsz - len, p);
			buf[len++] = '\t';
		} else {
			memcpy(buf + len, "tot\t", 4U);
			len += 4U;
		}
		return len;
	}

	static ssize_t ixtoistr(char *restrict buf, size_t bsz, size_t i)
	{
		size_t len = 0U;

		for (size_t j = arity; j--;) {
			const hx_t h = (i >> (5U * j)) % 32U;

			len += snprintf(buf + len, bsz - len, "%u", h);
			buf[len++] = '\t';
		}
		return len;
	}

	/* which index printer to use */
	ixtostr = !rawp ? ixtopstr : ixtoistr;

	/* calculate totals first, use 0 slot to store them */
	for (size_t i = 0U, n = 1ULL << (5U * arity); i < n; i += 32U) {
		chain[i] = 0U;
		for (size_t j = 1U; j < 32U; j++) {
			chain[i] += chain[i + j];
		}
	}
	for (size_t i = 0U, n = 1ULL << (5U * arity); i < n; i++) {
		char buf[256U];
		size_t len = 0U;

		if (LIKELY(!chain[i])) {
			continue;
		}
		len += ixtostr(buf + len, sizeof(buf) - len, i);
		len += cctostr(buf + len, sizeof(buf) - len, chain[i]);
		buf[len++] = '\n';
		fwrite(buf, 1, len, stdout);
	}
	return 0;
}


#include "mcfit.yucc"

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
				add_chain();
			}
		}

		free(line);

		/* print results */
		rc = prnt_chain(!!argi->raw_flag) < 0;
	}

	/* no need for that chain no more */
	free(chain);

out:
	yuck_free(argi);
	return rc;
}

/* mcfit.c ends here */
