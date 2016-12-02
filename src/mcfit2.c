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

typedef long unsigned int tv_t;
#define NOT_A_TIME	((tv_t)-1ULL)
#define MSECS		(1000)

/* hashed stamps (ilog2) and spreads (most significant digit plus exponent) */
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

static tv_t
strtotv(const char *ln, char **endptr)
{
	char *on;
	tv_t r;

	/* time value up first */
	with (long unsigned int s, x) {
		if (UNLIKELY((s = strtoul(ln, &on, 10), on == ln))) {
			return NOT_A_TIME;
		} else if (*on == '.') {
			char *moron;

			x = strtoul(++on, &moron, 10);
			if (UNLIKELY(moron - on > 9U)) {
				return NOT_A_TIME;
			} else if ((moron - on) % 3U) {
				/* huh? */
				return NOT_A_TIME;
			}
			switch (moron - on) {
			case 9U:
				x /= MSECS;
			case 6U:
				x /= MSECS;
			case 3U:
				break;
			case 0U:
			default:
				break;
			}
			on = moron;
		} else {
			x = 0U;
		}
		r = s * MSECS + x;
	}
	if (LIKELY(endptr != NULL)) {
		*endptr = on;
	}
	return r;
}

static ssize_t
tvtostr(char *restrict buf, size_t bsz, tv_t t)
{
	return snprintf(buf, bsz, "%lu.%03lu000000", t / MSECS, t % MSECS);
}

static __attribute__((const, pure)) hx_t
tvtohx(tv_t t)
{
/* just use the base2 log of t (resolution is milliseconds) */
	unsigned int x = 64U - __builtin_clzll(t);
	return -(x < 32U) & (hx_t)x;
}

static __attribute__((const, pure)) tv_t
hxtotv(hx_t h)
{
#if 0
/* round them off to good numbers */
	static tv_t r[] = {
		0U, 1U, 5U, 10U, 20U, 30U, 60U, 100U,
		200U, 500U, 1000U, 2000U, 5000U, 10000U, 15000U, 30000U,
		60000U, 120000U, 300000U, 600000U,
		1200000U, 1800000U, 3600000U, 7200000U,
		14400000U, 28800000U, 43200000U, 86400000U,
		259200000U, 604800000U, 1209600000U, 1814400000U,
	};
	return r[h];
#else
/* turns out we need to read the model back in */
	return (1ULL << h) - 1ULL;
#endif
}

static __attribute__((const, pure)) hx_t
pxtohx(px_t p, int q)
{
/* use most significant (base10) digit and exponent, alongside quantum */
	int e = 0;
	int v = (int)scalbnd32(frexpd32(p, &e), 1);
	const int f = (e - q - 1);
	int x = 10 * f + v - f;
	return -(v > 0 && x < 32) & (hx_t)x;
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


/* limit to arity of 4 (which would consume 8TB of RAM) */
static hx_t hist[4U];
static size_t nhist;
/* big counting array */
static cc_t *chain;
/* dataset quantum */
static int qunt;

static inline void
add_hist(hx_t ht, hx_t hs)
{
	hist[nhist++] = ht << 5U ^ hs;
	nhist %= arity;
	return;
}

static int
push_beef(char *ln, size_t UNUSED(lz))
{
	static tv_t last;
	char *on;
	tv_t t;
	px_t s;
	int q;

	if (UNLIKELY((t = strtotv(ln, &on)) == NOT_A_TIME)) {
		/* yeah right */
		return -1;
	} else if (UNLIKELY(*on++ != '\t')) {
		/* still not ok */
		return -1;
	}

	if (UNLIKELY((s = strtopx(on, &on), *on != '\n'))) {
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
	with (hx_t ht = tvtohx(t - last), hs = pxtohx(s, qunt)) {
		add_hist(ht, hs);
	}
	/* store stamp */
	last = t;
	return 0;
}

static int
add_chain(void)
{
	size_t k = 0U;

	for (size_t i = nhist; i < arity; i++) {
		k *= 32U * 32U;
		k += hist[i];
	}
	for (size_t i = 0U; i < nhist; i++) {
		k *= 32U * 32U;
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
			const hx_t h = (i >> (5U * 2U * j)) % (32U * 32U);
			const tv_t t = hxtotv(h / 32U);
			const px_t p = hxtopx(h % 32U, qunt);

			len += tvtostr(buf + len, bsz - len, t);
			buf[len++] = '\t';
			len += pxtostr(buf + len, bsz - len, p);
			buf[len++] = '\t';
		}
		/* last one is special */
		if (i % (32U * 32U)) {
			const hx_t h = i % (32U * 32U);
			const tv_t t = hxtotv(h / 32U);
			const px_t p = hxtopx(h % 32U, qunt);

			len += tvtostr(buf + len, bsz - len, t);
			buf[len++] = '\t';
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
			const hx_t h = (i >> (5U * 2U * j)) % (32U * 32U);

			len += snprintf(buf + len, bsz - len, "%u", h);
			buf[len++] = '\t';
		}
		return len;
	}

	/* which index printer to use */
	ixtostr = !rawp ? ixtopstr : ixtoistr;

	/* calculate totals first, use 0 slot to store them */
	for (size_t i = 0U, n = 1ULL << (5U * 2U * arity); i < n;
	     i += 32U * 32U) {
		chain[i] = 0U;
		for (size_t j = 1U; j < 32U * 32U; j++) {
			chain[i] += chain[i + j];
		}
	}
	for (size_t i = 0U, n = 1ULL << (5U * 2U * arity); i < n; i++) {
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


#include "mcfit2.yucc"

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
Info: arity=4 needs 8TB RAM, arity=5 would need 8PB RAM");
			rc = 1;
			goto out;
		}
	}

	/* initialise counting array */
	with (size_t z = 1ULL << (5U * 2U * arity)) {
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

/* mcfit2.c ends here */
