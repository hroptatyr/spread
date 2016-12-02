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
#include "pcg_basic.h"

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


/* limit to arity of 8 (which would consume 8TB of RAM) */
static hx_t hist[8U];
static size_t nhist;
/* big counting array */
static cc_t *chain;
/* dataset quantum */
static int qunt;

static int
make_chain(void)
{
/* uses ARITY global */
	size_t z = 1ULL << (5U * arity);
	chain = calloc(z, sizeof(*chain));
	return (chain != NULL) - 1U;
}

static void
free_chain(void)
{
	if (LIKELY(chain != NULL)) {
		free(chain);
	}
	return;
}

static inline void
add_hist(hx_t h)
{
	hist[nhist] = h;
	nhist = (nhist + 1U) % arity;
	return;
}

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
	add_hist(pxtohx(s, qunt));
	return 0;
}

static int
add_chain(size_t cnt)
{
	size_t k = 0U;

	for (size_t i = 0U; i < arity; i++) {
		k *= 32U;
		k += hist[i];
	}
	/* actually set the count now */
	chain[k] = cnt;
	return 0;
}

static int
read_model(const char *fn)
{
	char *line = NULL;
	size_t llen = 0U;
	ssize_t nrd;
	int rc = 0;
	FILE *fp;

	if (UNLIKELY((fp = fopen(fn, "r")) == NULL)) {
		/* just be abrupt */
		return -1;
	}

	/* first line is special */
	if ((nrd = getline(&line, &llen, fp)) <= 0) {
		/* yeah right */
		rc = -1;
		goto out;
	}
	/* count number of tabs, which should be our arity */
	arity = 0U;
	for (char *on = line;
	     (on = memchr(on, '\t', nrd - (on - line)));
	     arity++, on++);

	/* instantiate chain tensor */
	make_chain();

	do {
		char *on = line;
		long unsigned int cnt;

		for (size_t i = 0U; i < arity; i++, on++) {
			px_t s = strtopx(on, &on);
			int q;

			if (UNLIKELY(*on != '\t')) {
				on = memchr(on, '\t', nrd - (on - line));
				s = 0.df;
			} else if (UNLIKELY(qunt != (q = quantexpd32(s)))) {
				if (LIKELY(qunt)) {
					errno = 0, serror("\
Warning: quantum has changed  %i -> %i", qunt, q);
				}
				qunt = q;
			}
			hist[i] = pxtohx(s, qunt);
		}
		/* now read the count */
		cnt = strtoul(on, &on, 10U);

		/* bang to chain */
		add_chain(cnt);
	} while ((nrd = getline(&line, &llen, fp)) > 0);

	/* we don't trust the totals passed to us, recalc */
	for (size_t i = 0U, n = 1ULL << (5U * arity); i < n; i += 32U) {
		chain[i] = 0U;
		for (size_t j = 1U; j < 32U; j++) {
			chain[i] += chain[i + j];
		}
		for (size_t j = 1U; j < 32U; j++) {
			chain[i + j] = chain[i + j - 1U] - chain[i + j];
		}
	}
out:
	fclose(fp);
	return rc;
}

static int
boot_state(void)
{
	/* fit a model */
	char *line = NULL;
	size_t llen = 0UL;
	ssize_t nrd;

	while ((nrd = getline(&line, &llen, stdin)) > 0 && nhist < arity) {
		(void)push_beef(line, nrd);
	}

	free(line);
	fclose(stdin);
	return 0;
}

static int
init_rng(void)
{
	/* not quite random */
	pcg32_srandom(42U, 54U);
	return 0;
}

static int
infer1(void)
{
	size_t k = 0U;
	uint32_t r;

	for (size_t i = nhist; i < arity; i++) {
		k += hist[i];
		k *= 32U;
	}
	for (size_t i = 0U; i < nhist; i++) {
		k += hist[i];
		k *= 32U;
	}
	if (UNLIKELY(!chain[k])) {
		/* we're in a state we don't know about */
		return -1;
	}
	r = pcg32_boundedrand(chain[k]);
	for (hx_t j = 1U; j < 32U; j++) {
		if (r >= chain[k + j]) {
			char buf[32U];
			size_t len;

			len = pxtostr(buf, sizeof(buf), hxtopx(j, qunt));
			buf[len++] = '\n';
			fwrite(buf, 1, len, stdout);

			/* keep state */
			add_hist(j);
			break;
		}
	}
	return 0;
}


#include "mcmc.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	} else if (!argi->nargs) {
		errno = 0, serror("\
Error: model file mandatory");
		rc = 1;
		goto out;
	} else if (UNLIKELY(read_model(*argi->args) < 0)) {
		errno = 0, serror("\
Error: cannot read model file `%s'", *argi->args);
		rc = 1;
		goto out;
	}

	/* kick off bootstrapping */
	memset(hist, 0, sizeof(hist));
	nhist = 0U;
	/* arity becomes the input dimension */
	arity--;
	if (!isatty(STDIN_FILENO)) {
		boot_state();
	}

	/* seed the RNG */
	init_rng();
	/* run the simulation */
	while (!(infer1() < 0));

	/* no need for that chain no more */
	free_chain();

out:
	yuck_free(argi);
	return rc;
}

/* mcmc.c ends here */
