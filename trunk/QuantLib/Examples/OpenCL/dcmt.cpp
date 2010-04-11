/* Copyright (C) 2001-2009 Makoto Matsumoto and Takuji Nishimura.  */
/* Copyright (C) 2009 Mutsuo Saito                                 */
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
/* 02111-1307  USA                                                 */

#include "dcmt.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <cstring> //required for memcmp()
#include <boost/cstdint.hpp>

static const int irredpolylist[NIRREDPOLY][MAX_IRRED_DEG+1] = {
    {0,1,0,0,0,0,0,0,0,0,},{1,1,0,0,0,0,0,0,0,0,},{1,1,1,0,0,0,0,0,0,0,},
    {1,1,0,1,0,0,0,0,0,0,},{1,0,1,1,0,0,0,0,0,0,},{1,1,0,0,1,0,0,0,0,0,},
    {1,0,0,1,1,0,0,0,0,0,},{1,1,1,1,1,0,0,0,0,0,},{1,0,1,0,0,1,0,0,0,0,},
    {1,0,0,1,0,1,0,0,0,0,},{1,1,1,1,0,1,0,0,0,0,},{1,1,1,0,1,1,0,0,0,0,},
    {1,1,0,1,1,1,0,0,0,0,},{1,0,1,1,1,1,0,0,0,0,},{1,1,0,0,0,0,1,0,0,0,},
    {1,0,0,1,0,0,1,0,0,0,},{1,1,1,0,1,0,1,0,0,0,},{1,1,0,1,1,0,1,0,0,0,},
    {1,0,0,0,0,1,1,0,0,0,},{1,1,1,0,0,1,1,0,0,0,},{1,0,1,1,0,1,1,0,0,0,},
    {1,1,0,0,1,1,1,0,0,0,},{1,0,1,0,1,1,1,0,0,0,},{1,1,0,0,0,0,0,1,0,0,},
    {1,0,0,1,0,0,0,1,0,0,},{1,1,1,1,0,0,0,1,0,0,},{1,0,0,0,1,0,0,1,0,0,},
    {1,0,1,1,1,0,0,1,0,0,},{1,1,1,0,0,1,0,1,0,0,},{1,1,0,1,0,1,0,1,0,0,},
    {1,0,0,1,1,1,0,1,0,0,},{1,1,1,1,1,1,0,1,0,0,},{1,0,0,0,0,0,1,1,0,0,},
    {1,1,0,1,0,0,1,1,0,0,},{1,1,0,0,1,0,1,1,0,0,},{1,0,1,0,1,0,1,1,0,0,},
    {1,0,1,0,0,1,1,1,0,0,},{1,1,1,1,0,1,1,1,0,0,},{1,0,0,0,1,1,1,1,0,0,},
    {1,1,1,0,1,1,1,1,0,0,},{1,0,1,1,1,1,1,1,0,0,},{1,1,0,1,1,0,0,0,1,0,},
    {1,0,1,1,1,0,0,0,1,0,},{1,1,0,1,0,1,0,0,1,0,},{1,0,1,1,0,1,0,0,1,0,},
    {1,0,0,1,1,1,0,0,1,0,},{1,1,1,1,1,1,0,0,1,0,},{1,0,1,1,0,0,1,0,1,0,},
    {1,1,1,1,1,0,1,0,1,0,},{1,1,0,0,0,1,1,0,1,0,},{1,0,1,0,0,1,1,0,1,0,},
    {1,0,0,1,0,1,1,0,1,0,},{1,0,0,0,1,1,1,0,1,0,},{1,1,1,0,1,1,1,0,1,0,},
    {1,1,0,1,1,1,1,0,1,0,},{1,1,1,0,0,0,0,1,1,0,},{1,1,0,1,0,0,0,1,1,0,},
    {1,0,1,1,0,0,0,1,1,0,},{1,1,1,1,1,0,0,1,1,0,},{1,1,0,0,0,1,0,1,1,0,},
    {1,0,0,1,0,1,0,1,1,0,},{1,0,0,0,1,1,0,1,1,0,},{1,0,1,1,1,1,0,1,1,0,},
    {1,1,0,0,0,0,1,1,1,0,},{1,1,1,1,0,0,1,1,1,0,},{1,1,1,0,1,0,1,1,1,0,},
    {1,0,1,1,1,0,1,1,1,0,},{1,1,1,0,0,1,1,1,1,0,},{1,1,0,0,1,1,1,1,1,0,},
    {1,0,1,0,1,1,1,1,1,0,},{1,0,0,1,1,1,1,1,1,0,},{1,1,0,0,0,0,0,0,0,1,},
    {1,0,0,0,1,0,0,0,0,1,},{1,1,1,0,1,0,0,0,0,1,},{1,1,0,1,1,0,0,0,0,1,},
    {1,0,0,0,0,1,0,0,0,1,},{1,0,1,1,0,1,0,0,0,1,},{1,1,0,0,1,1,0,0,0,1,},
    {1,1,0,1,0,0,1,0,0,1,},{1,0,0,1,1,0,1,0,0,1,},{1,1,1,1,1,0,1,0,0,1,},
    {1,0,1,0,0,1,1,0,0,1,},{1,0,0,1,0,1,1,0,0,1,},{1,1,1,1,0,1,1,0,0,1,},
    {1,1,1,0,1,1,1,0,0,1,},{1,0,1,1,1,1,1,0,0,1,},{1,1,1,0,0,0,0,1,0,1,},
    {1,0,1,0,1,0,0,1,0,1,},{1,0,0,1,1,0,0,1,0,1,},{1,1,0,0,0,1,0,1,0,1,},
    {1,0,1,0,0,1,0,1,0,1,},{1,1,1,1,0,1,0,1,0,1,},{1,1,1,0,1,1,0,1,0,1,},
    {1,0,1,1,1,1,0,1,0,1,},{1,1,1,1,0,0,1,1,0,1,},{1,0,0,0,1,0,1,1,0,1,},
    {1,1,0,1,1,0,1,1,0,1,},{1,0,1,0,1,1,1,1,0,1,},{1,0,0,1,1,1,1,1,0,1,},
    {1,0,0,0,0,0,0,0,1,1,},{1,1,0,0,1,0,0,0,1,1,},{1,0,1,0,1,0,0,0,1,1,},
    {1,1,1,1,1,0,0,0,1,1,},{1,1,0,0,0,1,0,0,1,1,},{1,0,0,0,1,1,0,0,1,1,},
    {1,1,0,1,1,1,0,0,1,1,},{1,0,0,1,0,0,1,0,1,1,},{1,1,1,1,0,0,1,0,1,1,},
    {1,1,0,1,1,0,1,0,1,1,},{1,0,0,0,0,1,1,0,1,1,},{1,1,0,1,0,1,1,0,1,1,},
    {1,0,1,1,0,1,1,0,1,1,},{1,1,0,0,1,1,1,0,1,1,},{1,1,1,1,1,1,1,0,1,1,},
    {1,0,1,0,0,0,0,1,1,1,},{1,1,1,1,0,0,0,1,1,1,},{1,0,0,0,0,1,0,1,1,1,},
    {1,0,1,0,1,1,0,1,1,1,},{1,0,0,1,1,1,0,1,1,1,},{1,1,1,0,0,0,1,1,1,1,},
    {1,1,0,1,0,0,1,1,1,1,},{1,0,1,1,0,0,1,1,1,1,},{1,0,1,0,1,0,1,1,1,1,},
    {1,0,0,1,1,0,1,1,1,1,},{1,1,0,0,0,1,1,1,1,1,},{1,0,0,1,0,1,1,1,1,1,},
    {1,1,0,1,1,1,1,1,1,1,},
};

static const boost::uint8_t pivot_calc_tbl[256] = {
    0, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    2, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    1, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    2, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    3, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
    4, 8, 7, 8, 6, 8, 7, 8, 5, 8, 7, 8, 6, 8, 7, 8,
};

static int proper_mersenne_exponent(int p)
{
    switch(p) {
    case 521:
    case 607:
    case 1279:
    case 2203:
    case 2281:
    case 3217:
    case 4253:
    case 4423:
    case 9689:
    case 9941:
    case 11213:
    case 19937:
    case 21701:
    case 23209:
    case 44497:
	return 1;
    default:
	return 0;
    }
}

static mt_params* alloc_mt_params(int n)
{
    mt_params *mtp;

    mtp = (mt_params*)malloc(sizeof(mt_params));
    if (NULL == mtp) return NULL;
    mtp->state = (uint32_t*)malloc(n*sizeof(uint32_t));
    if (NULL == mtp->state) {
	free(mtp);
	return NULL;
    }

    return mtp;
}

static void make_masks(int r, int w, mt_params *mtp)
{
    int i;
    uint32_t ut, wm, um, lm;

    wm = 0xFFFFFFFF;
    wm >>= (WORDLENGTH - w);

    ut = 0;
    for (i=0; i<r; i++) {
	ut <<= 1;
	ut |= LSB;
    }

    lm = ut;
    um = (~ut) & wm;

    mtp->wmask = wm;
    mtp->umask = um;
    mtp->lmask = lm;
}

static dcmt_Polynomial *NewPoly(int degree)
{
    dcmt_Polynomial *p;

    p = (dcmt_Polynomial *)calloc( 1, sizeof(dcmt_Polynomial));
    if( p==NULL ){
	printf("calloc error in \"NewPoly()\"\n");
	exit(1);
    }
    p->deg = degree;

    if (degree < 0) {
	p->x = NULL;
	return p;
    }

    p->x = (int *)calloc( degree + 1, sizeof(int));
    if( p->x == NULL ){
	printf("calloc error\n");
	exit(1);
    }

    return p;
}

static void FreePoly( dcmt_Polynomial *p)
{
    if (p->x != NULL)
	free( p->x );
    free( p );
}

static dcmt_Polynomial *make_tntm( int n, int m)
{
    dcmt_Polynomial *p;

    p = NewPoly(n);
    p->x[n] = p->x[m] = 1;

    return p;
}

static dcmt_Polynomial *PolynomialDup(dcmt_Polynomial *pl)
{
    dcmt_Polynomial *pt;
    int i;

    pt = NewPoly(pl->deg);
    for (i=pl->deg; i>=0; i--)
	pt->x[i] = pl->x[i];

    return pt;
}

static dcmt_Polynomial *PolynomialMult(dcmt_Polynomial *p0,dcmt_Polynomial *p1)
{
    int i, j;
    dcmt_Polynomial *p;

    /* if either p0 or p1 is 0, return 0 */
    if ( (p0->deg < 0) || (p1->deg < 0) ) {
	p = NewPoly(-1);
	return p;
    }

    p = NewPoly(p0->deg + p1->deg);
    for( i=0; i<=p1->deg; i++){
	if( p1->x[i] ){
	    for( j=0; j<=p0->deg; j++){
		p->x[i+j] ^= p0->x[j];
	    }
	}
    }

    return p;
}

static void MakepreModPolys(prescr_t *pre, int mm, int nn, int rr, int ww)
{
    dcmt_Polynomial *t, *t0, *t1, *s, *s0, *s1;
    int i,j;

    j = 0;
    t = NewPoly(0);
    t->deg = 0;
    t->x[0] = 1;
    pre->preModPolys[j++] = t;

    t = make_tntm (nn, mm);
    t0 = make_tntm (nn, mm);
    s = make_tntm (nn-1, mm-1);

    for( i=1; i<(ww - rr); i++){
	pre->preModPolys[j++] = PolynomialDup(t0);
	t1 = t0;
	t0 = PolynomialMult(t0, t);
	FreePoly(t1);
    }

    pre->preModPolys[j++] = PolynomialDup(t0);

    s0 =PolynomialMult( t0, s);
    FreePoly(t0);	FreePoly(t);
    for( i=(rr-2); i>=0; i--){
	pre->preModPolys[j++] = PolynomialDup(s0);
	s1 = s0;
	s0 = PolynomialMult( s0, s);
	FreePoly(s1);
    }

    pre->preModPolys[j++] = PolynomialDup(s0);

    FreePoly(s0); FreePoly(s);
}

void NextIrredPoly(dcmt_Polynomial *pl, int nth)
{
    int i, max_deg;

    for (max_deg=0,i=0; i<=MAX_IRRED_DEG; i++) {
	if ( irredpolylist[nth][i] )
	    max_deg = i;
	pl->x[i] = irredpolylist[nth][i];
    }

    pl->deg = max_deg;
}

static void PolynomialMod( dcmt_Polynomial *wara, const dcmt_Polynomial *waru)
{
    int i;
    int deg_diff;

    while( wara->deg >= waru->deg  ){
	deg_diff = wara->deg - waru->deg;
	for( i=0; i<=waru->deg; i++){
	    wara->x[ i+deg_diff ] ^= waru->x[i];
	}

	for( i=wara->deg; i>=0; i--){
	    if( wara->x[i] ) break;
	}
	wara->deg=i;

    }
}

static uint32_t word2bit(dcmt_Polynomial *pl)
{
    int i;
    uint32_t bx;

    bx = 0;
    for (i=pl->deg; i>0; i--) {
	if (pl->x[i]) bx |= 0x1;
	bx <<= 1;
    }
    if (pl->x[0]) bx |= 0x1;

    return bx;
}

static void makemodlist(prescr_t *pre, dcmt_Polynomial *pl, int nPoly)
{
    dcmt_Polynomial *tmpPl;
    int i;

    for (i=0; i<=pre->sizeofA; i++) {
	tmpPl = PolynomialDup(pre->preModPolys[i]);
	PolynomialMod(tmpPl,pl);
	pre->modlist[nPoly][i] = word2bit(tmpPl);
	FreePoly(tmpPl);
    }
}

void _InitPrescreening_dc(prescr_t *pre, int m, int n, int r, int w)
{
    int i;
    dcmt_Polynomial *pl;

    pre->sizeofA = w;

    pre->preModPolys = (dcmt_Polynomial **)malloc(
	(pre->sizeofA+1)*(sizeof(dcmt_Polynomial*)));
    if (NULL == pre->preModPolys) {
	printf ("malloc error in \"InitPrescreening\"\n");
	exit(1);
    }
    MakepreModPolys(pre, m,n,r,w);

    pre->modlist = (uint32_t**)malloc(NIRREDPOLY * sizeof(uint32_t*));
    if (NULL == pre->modlist) {
	printf ("malloc error in \"InitPrescreening()\"\n");
	exit(1);
    }
    for (i=0; i<NIRREDPOLY; i++) {
	pre->modlist[i]
	    = (uint32_t*)malloc( (pre->sizeofA + 1) * (sizeof(uint32_t)) );
	if (NULL == pre->modlist[i]) {
	    printf ("malloc error in \"InitPrescreening()\"\n");
	    exit(1);
	}
    }

    for (i=0; i<NIRREDPOLY; i++) {
	pl = NewPoly(MAX_IRRED_DEG);
	NextIrredPoly(pl,i);
	makemodlist(pre, pl, i);
	FreePoly(pl);
    }

    for (i=pre->sizeofA; i>=0; i--)
	FreePoly(pre->preModPolys[i]);
    free(pre->preModPolys);
}

void _InitCheck32_dc(check32_t *ck, int r, int w)
{
    int i;

    /* word_mask (least significant w bits) */
    ck->word_mask = 0xFFFFFFFF;
    ck->word_mask <<= WORDLENGTH - w;
    ck->word_mask >>= WORDLENGTH - w;
    /* lower_mask (least significant r bits) */
    for (ck->lower_mask=0,i=0; i<r; ++i) {
	ck->lower_mask <<= 1;
	ck->lower_mask |= LSB;
    }
    /* upper_mask (most significant (w-r) bits */
    ck->upper_mask = (~ck->lower_mask) & ck->word_mask;
}

static mt_params *init_mt_search(check32_t *ck, prescr_t *pre, int w, int p)
{
    int n, m, r;
    mt_params *mtp;

    if ( (w>32) || (w<31) ) {
	printf ("Sorry, currently only w = 32 or 31 is allowded.\n");
	return NULL;
    }

    if ( !proper_mersenne_exponent(p) ) {
	if (p<521) {
	    printf ("\"p\" is too small.\n");
	    return NULL;
	}
	else if (p>44497){
	    printf ("\"p\" is too large.\n");
	    return NULL;
	}
	else {
	    printf ("\"p\" is not a Mersenne exponent.\n");
	    return NULL;
	}
    }

    n = p/w + 1; /* since p is Mersenne Exponent, w never divids p */
    mtp = alloc_mt_params(n);
    if (NULL == mtp) return NULL;

    m = n/2;
    if (m < 2) m = n-1;
    r = n * w - p;

    make_masks(r, w, mtp);
    _InitPrescreening_dc(pre, m, n, r, w);
    _InitCheck32_dc(ck, r, w);

    mtp->mm = m;
    mtp->nn = n;
    mtp->rr = r;
    mtp->ww = w;

    return mtp;
}

/* Initializing the array with a seed */
void _sgenrand_dc(_org_state *st, uint32_t seed) {
    int i;

    for (i=0;i<MT_NN;i++) {
		st->mt[i] = seed;
        seed = (1812433253UL * (seed  ^ (seed >> 30))) + i + 1;
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
    }
    st->mti = MT_NN;
}

static void copy_mt_params(mt_params *src, mt_params *dst)
{
    dst->nn = src->nn;
    dst->mm = src->mm;
    dst->rr = src->rr;
    dst->ww = src->ww;
    dst->wmask = src->wmask;
    dst->umask = src->umask;
    dst->lmask = src->lmask;
}

uint32_t _genrand_dc(_org_state *st)
{
    uint32_t x;
    static const uint32_t mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (st->mti >= MT_NN) { /* generate N words at one time */
        int kk;

        for (kk=0;kk<MT_NN-MT_MM;kk++) {
            x = (st->mt[kk]&MT_UMASK)|(st->mt[kk+1]&MT_LMASK);
            st->mt[kk] = st->mt[kk+MT_MM] ^ (x >> 1) ^ mag01[x & 0x1];
        }
        for (;kk<MT_NN-1;kk++) {
            x = (st->mt[kk]&MT_UMASK)|(st->mt[kk+1]&MT_LMASK);
            st->mt[kk] = st->mt[kk+(MT_MM-MT_NN)] ^ (x >> 1) ^ mag01[x & 0x1];
        }
        x = (st->mt[MT_NN-1]&MT_UMASK)|(st->mt[0]&MT_LMASK);
        st->mt[MT_NN-1] = st->mt[MT_MM-1] ^ (x >> 1) ^ mag01[x & 0x1];

        st->mti = 0;
    }

    x = st->mt[st->mti++];
    x ^= TEMPERING_SHIFT_U(x);
    x ^= TEMPERING_SHIFT_S(x) & TEMPERING_MASK_B;
    x ^= TEMPERING_SHIFT_T(x) & TEMPERING_MASK_C;
    x ^= TEMPERING_SHIFT_L(x);

    return x;
}

static uint32_t nextA(_org_state *org, int w)
{
    uint32_t x, word_mask;

    word_mask = 0xFFFFFFFF;
    word_mask <<= WORDLENGTH - w;
    word_mask >>= WORDLENGTH - w;

    x = _genrand_dc(org);
    x &= word_mask;
    x |= (LSB << (w-1));

    return x;
}

static uint32_t nextA_id(_org_state *org, int w, int id, int idw)
{
    uint32_t x, word_mask;

    word_mask = 0xFFFFFFFF;
    word_mask <<= WORDLENGTH - w;
    word_mask >>= WORDLENGTH - w;
    word_mask >>= idw;
    word_mask <<= idw;

    x = _genrand_dc(org);
    x &= word_mask;
    x |= (LSB << (w-1));
    x |= (uint32_t)id; /* embedding id */

    return x;
}

static int IsReducible(prescr_t *pre, uint32_t aaa, uint32_t *polylist)
{
    int i;
    uint32_t x;

    x = polylist[pre->sizeofA];
    for (i=pre->sizeofA-1; i>=0; i--) {
	if (aaa&0x1)
	    x ^= polylist[i];
	aaa >>= 1;
    }

    if ( x == 0 ) return REDU;
    else return NONREDU;
}

int _prescreening_dc(prescr_t *pre, uint32_t aaa)
{

    int i;

    for (i=0; i<NIRREDPOLY; i++) {
	if (IsReducible(pre, aaa,pre->modlist[i])==REDU)
	    return REJECTED;
    }
    return NOT_REJECTED;
}

int _CheckPeriod_dc(check32_t *ck, _org_state *st, uint32_t a, int m, int n, int r, int w)
{
    int i, j, p, pp;
    uint32_t y, *x, *init, mat[2];


    p = n*w-r;
    x = (uint32_t*) malloc (2*p*sizeof(uint32_t));
    if (NULL==x) {
	printf("malloc error in \"_CheckPeriod_dc()\"\n");
	exit(1);
    }

    init = (uint32_t*) malloc (n*sizeof(uint32_t));
    if (NULL==init) {
	printf("malloc error \"_CheckPeriod_dc()\"\n");
	free(x);
	exit(1);
    }

    /* set initial values */
    for (i=0; i<n; ++i)
	x[i] = init[i] = (ck->word_mask & _genrand_dc(st));
    /* it is better that LSBs of x[2] and x[3] are different */
    if ( (x[2]&LSB) == (x[3]&LSB) ) {
	x[3] ^= 1;
	init[3] ^= 1;
    }

    pp = 2*p-n;
    mat[0] = 0; mat[1] = a;
    for (j=0; j<p; ++j) {

	/* generate */
	for (i=0; i<pp; ++i){
	    y = (x[i]&ck->upper_mask) | (x[i+1]&ck->lower_mask);
	    x[i+n] = x[i+m] ^ ( (y>>1) ^ mat[y&LSB] );
	}

	/* pick up odd subscritpt elements */
	for (i=2; i<=p; ++i)
	    x[i] = x[(i<<1)-1];

	/* reverse generate */
	for (i=p-n; i>=0; --i) {
	    y = x[i+n] ^ x[i+m] ^ mat[ x[i+1]&LSB ];
	    y <<=1; y |= x[i+1]&LSB;

	    x[i+1] = (x[i+1]&ck->upper_mask) | (y&ck->lower_mask);
	    x[i] = (y&ck->upper_mask) | (x[i]&ck->lower_mask);
	}

    }

    if ((x[0]&ck->upper_mask)==(init[0]&ck->upper_mask)) {
	for (i=1; i<n; ++i) {
	    if (x[i] != init[i])
		break;
	}
	if (i==n) {
	    free(x); free(init);
	    return IRRED;
	}
    }


    free(x); free(init);
    return REDU;
}

static int get_irred_param(check32_t *ck, prescr_t *pre, _org_state *org, mt_params *mtp, int id, int idw)
{
    int i;
    uint32_t a;

    for (i=0; i<MAX_SEARCH; i++) {
	if (idw == 0)
	    a = nextA(org, mtp->ww);
	else
	    a = nextA_id(org, mtp->ww, id, idw);
	if (NOT_REJECTED == _prescreening_dc(pre, a) ) {
	    if (IRRED
		== _CheckPeriod_dc(ck, org, a,mtp->mm,mtp->nn,mtp->rr,mtp->ww)) {
		mtp->aaa = a;
		break;
	    }
	}
    }

    if (MAX_SEARCH == i) return NOT_FOUND;
    return FOUND;
}

void free_mt_params(mt_params *mtp)
{
    free(mtp->state);
    free(mtp);
}

static void init_tempering(eqdeg_t *eq, mt_params *mtp)
{
    int i;

    eq->mmm = mtp->mm;
    eq->nnn = mtp->nn;
    eq->rrr = mtp->rr;
    eq->www = mtp->ww;
    eq->shift_0 = S00;
    eq->shift_1 = S01;
    eq->shift_s = SSS;
    eq->shift_t = TTT;
    eq->ggap = WORDLENGTH - eq->www;
    /* bits are filled in mts->aaa from MSB */
    eq->aaa[0] = 0; eq->aaa[1] = (mtp->aaa) << eq->ggap;


    for( i=0; i<WORDLENGTH; i++)
        eq->bitmask[i] = 0x80000000UL >> i;

    for( i=0, eq->glower_mask=0; i<eq->rrr; i++)
	eq->glower_mask = (eq->glower_mask<<1)| 0x1;

    eq->gupper_mask = ~eq->glower_mask;
    eq->gupper_mask <<= eq->ggap;
    eq->glower_mask <<= eq->ggap;

    eq->greal_mask = (eq->gupper_mask | eq->glower_mask);
}

static int push_mask(eqdeg_t *eq, int l, int v, uint32_t b, uint32_t c, uint32_t *bbb, uint32_t *ccc)
{
    int i, j, k, nbv, nbvt;
    uint32_t bmask, bv_buf[2], bvt_buf[2];

    k = l;
    if( (eq->shift_s+v) >= eq->www ){
        nbv = 1; bv_buf[0] = 0;
    }
    else if( (v>=eq->shift_t) && (c&eq->bitmask[v-eq->shift_t] ) ){
        nbv = 1; bv_buf[0] = b&eq->bitmask[v];
    }
    else {
        nbv = 2; bv_buf[0] = eq->bitmask[v]; bv_buf[1] = 0;
    }

    if( ((v+eq->shift_t+eq->shift_s) < eq->www) && (c&eq->bitmask[v]) ){
        nbvt = 2; bvt_buf[0] = eq->bitmask[v+eq->shift_t]; bvt_buf[1] = 0;
    }
    else {
        nbvt = 1; bvt_buf[0] = 0;
    }

    bmask = eq->bitmask[v];
    if( (v+eq->shift_t) < eq->www )
        bmask |= eq->bitmask[v+eq->shift_t];
    bmask = ~bmask;
    for( i=0; i<nbvt; ++i){
        for( j=0; j<nbv; ++j){
            bbb[k] = (b&bmask) | bv_buf[j] | bvt_buf[i];
            ccc[k] = c;
            ++k;
        }
    }

    return k-l;
}

static int push_stack(eqdeg_t *eq, uint32_t b, uint32_t c, int v, uint32_t *bbb, uint32_t *ccc)
{
    int i, ll, ncv;
    uint32_t cv_buf[2];

    ll = 0;

    if( (v+eq->shift_t) < eq->www ){
        ncv = 2; cv_buf[0] = c | eq->bitmask[v]; cv_buf[1] = c;
    }
    else {
        ncv = 1; cv_buf[0] = c;
    }

    for( i=0; i<ncv; ++i)
        ll += push_mask(eq, ll, v, b, cv_buf[i], bbb, ccc);

    return ll;
}

static int pivot_reduction(eqdeg_t *eq, int v)
{
    Vector **lattice, *ltmp;
    int i;
    int pivot;
    int count;
    int min;

    eq->upper_v_bits = 0;
    for( i=0; i<v; i++) {
        eq->upper_v_bits |= eq->bitmask[i];
    }

    lattice = make_lattice(eq, v );

    for (;;) {
	pivot = calc_pivot(lattice[v]->next);
	if (lattice[pivot]->count < lattice[v]->count) {
	    ltmp = lattice[pivot];
	    lattice[pivot] = lattice[v];
	    lattice[v] = ltmp;
	}
	add(eq->nnn, lattice[v], lattice[pivot]);
	if (lattice[v]->next == 0) {
	    count = 0;
	    next_state(eq, lattice[v], &count);
	    if (lattice[v]->next == 0) {
		if (is_zero(eq->nnn, lattice[v])) {
		    break;
		}
		while (lattice[v]->next == 0) {
		    count++;
		    next_state(eq, lattice[v], &count);
		    if (count > eq->nnn * (eq->www-1) - eq->rrr) {
			break;
		    }
		}
		if (lattice[v]->next == 0) {
		    break;
		}
	    }
	}
    }

    min = lattice[0]->count;
    for (i = 1; i < v; i++) {
	if (min > lattice[i]->count) {
	    min = lattice[i]->count;
	}
    }
    free_lattice( lattice, v );
    return min;
}

static MaskNode *optimize_v_hard(eqdeg_t *eq, int v, MaskNode *prev_masks)
{
    int i, ll, t;
    uint32_t bbb[8], ccc[8];
    MaskNode *cur_masks;

    cur_masks = NULL;

    while (prev_masks != NULL) {

	ll = push_stack(eq, prev_masks->b,prev_masks->c,v,bbb,ccc);

	for (i=0; i<ll; ++i) {
	    eq->mask_b = bbb[i];
	    eq->mask_c = ccc[i];
	    t = pivot_reduction(eq, v+1);
	    if (t >= eq->gcur_maxlengs[v]) {
		eq->gcur_maxlengs[v] = t;
		eq->gmax_b = eq->mask_b;
		eq->gmax_c = eq->mask_c;
		cur_masks = cons_MaskNode(cur_masks, eq->mask_b, eq->mask_c, t);
	    }
	}
	prev_masks = prev_masks->next;
    }

    cur_masks = delete_lower_MaskNodes(cur_masks, eq->gcur_maxlengs[v]);

    return cur_masks;
}

void _get_tempering_parameter_dc(mt_params *mtp)
{
    eqdeg_t eq;
    init_tempering(&eq, mtp);
    optimize_v(&eq, 0, 0, 0);
    mtp->shift0 = eq.shift_0;
    mtp->shift1 = eq.shift_1;
    mtp->shiftB = eq.shift_s;
    mtp->shiftC = eq.shift_t;
    mtp->maskB = eq.mask_b >> eq.ggap;
    mtp->maskC = eq.mask_c >> eq.ggap;
}

mt_params** get_mt_parameters_st(int w, int p, int start_id, int max_id, uint32_t seed, int *count) {
    mt_params **mtpp, *template_mtp;
    int i;
    prescr_t pre;
    _org_state org;
    check32_t ck;

    if ((start_id > max_id) || (max_id > 0xffff) || (start_id < 0)) {
	printf("\"id\" error\n");
	return NULL;
    }

    _sgenrand_dc(&org, seed);
    mtpp = (mt_params**)malloc(sizeof(mt_params*)*(max_id-start_id+1));
    if (NULL == mtpp) return NULL;

    template_mtp = init_mt_search(&ck, &pre, w, p);
    if (template_mtp == NULL) {
	free(mtpp);
	return NULL;
    }
    *count = 0;
    for (i=0; i<=max_id-start_id; i++) {
	mtpp[i] = alloc_mt_params(template_mtp->nn);
	if (NULL == mtpp[i]) {
	    break;
	}

	copy_mt_params(template_mtp, mtpp[i]);

	if ( NOT_FOUND == get_irred_param(&ck, &pre, &org, mtpp[i],
					  i+start_id,DEFAULT_ID_SIZE) ) {
	    free_mt_params(mtpp[i]);
	    break;
	}
	_get_tempering_parameter_hard_dc(mtpp[i]);
	++(*count);
    }

    free_mt_params(template_mtp);
    end_mt_search(&pre);
    if (*count > 0) {
	return mtpp;
    } else {
	free(mtpp);
	return NULL;
    }
}

static void free_lattice( Vector **lattice, int v)
{
    int i;

    for( i=0; i<=v; i++)
        free_Vector( lattice[i] );
    free( lattice );
}

/* adds v to u (then u will change) */
static void add(int nnn, Vector *u, Vector *v)
{
    int i;
    int diff = (v->start - u->start + nnn) % nnn;
    for (i = 0; i < nnn - diff; i++) {
	u->cf[i] ^= v->cf[i + diff];
    }
    diff = diff - nnn;
    for (; i < nnn; i++) {
	u->cf[i] ^= v->cf[i + diff];
    }
    u->next ^=  v->next;
}

static Vector **make_lattice(eqdeg_t *eq, int v)
{
    int i;
    int count;
    Vector **lattice, *bottom;

    lattice = (Vector **)malloc( (v+1) * sizeof( Vector *) );
    if( NULL == lattice ){
        printf("malloc error in \"make_lattice\"\n");
        exit(1);
    }

    for( i=0; i<v; i++){ /* from 0th row to v-1-th row */
        lattice[i] = new_Vector(eq->nnn);
        lattice[i]->next = eq->bitmask[i];
        lattice[i]->start = 0;
        lattice[i]->count = 0;
    }

    bottom = new_Vector(eq->nnn); /* last row */
    for(i=0; i< eq->nnn; i++) {
	bottom->cf[i] = 0;
    }
    bottom->cf[eq->nnn -1] = 0xc0000000 & eq->greal_mask;
    bottom->start = 0;
    bottom->count = 0;
    count = 0;
    do {
	next_state(eq, bottom, &count);
    } while (bottom->next == 0);
//    degree_of_vector(eq, top );
    lattice[v] = bottom;

    return lattice;
}

static void next_state(eqdeg_t *eq, Vector *v, int *count) {
    uint32_t tmp;

    do {
	tmp = ( v->cf[v->start] & eq->gupper_mask )
	    | ( v->cf[(v->start + 1) % eq->nnn] & eq->glower_mask );
	v->cf[v->start] = v->cf[(v->start + eq->mmm) % eq->nnn]
	    ^ ( (tmp>>1) ^ eq->aaa[lsb(eq, tmp)] );
	v->cf[v->start] &= eq->greal_mask;
	tmp = v->cf[v->start];
	v->start = (v->start + 1) % eq->nnn;
	v->count++;
	tmp = trnstmp(eq, tmp);
	tmp = masktmp(eq, tmp);
	v->next = tmp & eq->upper_v_bits;
	(*count)++;
	if (*count > eq->nnn * (eq->www-1) - eq->rrr) {
	    break;
	}
    } while (v->next == 0);
}

/***********/
static MaskNode *cons_MaskNode(MaskNode *head, uint32_t b, uint32_t c, int leng)
{
    MaskNode *t;

    t = (MaskNode*)malloc(sizeof(MaskNode));
    if (t == NULL) {
	printf("malloc error in \"cons_MaskNode\"\n");
        exit(1);
    }

    t->b = b;
    t->c = c;
    t->leng = leng;
    t->next = head;

    return t;
}

static int calc_pivot(uint32_t v) {
    int p1, p2, p3, p4;

    p1 = pivot_calc_tbl[v & 0xff];
    if (p1) {
	return p1 + 24 - 1;
    }
    p2 = pivot_calc_tbl[(v >> 8) & 0xff];
    if (p2) {
	return p2 + 16 - 1;
    }
    p3 = pivot_calc_tbl[(v >> 16) & 0xff];
    if (p3) {
	return p3 + 8 - 1;
    }
    p4 = pivot_calc_tbl[(v >> 24) & 0xff];
    if (p4) {
	return p4 - 1;
    }
    return -1;
}

static int is_zero(int size, Vector *v) {
    if (v->cf[0] != 0) {
	return 0;
    } else {
	return (memcmp(v->cf, v->cf + 1, sizeof(uint32_t) * (size - 1)) == 0);
    }
}

static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l)
{
    MaskNode *s, *t, *tail;

    s = head;
    while(1) { /* heading */
	if (s == NULL)
	    return NULL;
	if (s->leng >= l)
	    break;
	t = s->next;
	free(s);
	s = t;
    }

    head = tail = s;

    while (head != NULL) {
	t = head->next;
	if (head->leng < l) {
	    free(head);
	}
	else {
	    tail->next = head;
	    tail = head;
	}
	head = t;
    }

    tail->next = NULL;
    return s;
}

static void optimize_v(eqdeg_t *eq, uint32_t b, uint32_t c, int v)
{
    int i, max_len, max_i, ll, t;
    uint32_t bbb[8], ccc[8];

    ll = push_stack(eq, b,c,v,bbb,ccc);

    max_len = max_i = 0;
    if (ll > 1) {
	for (i=0; i<ll; ++i) {
	    eq->mask_b = bbb[i];
	    eq->mask_c = ccc[i];
	    t = pivot_reduction(eq, v+1);
	    if (t > max_len) {
		max_len = t;
		max_i = i;
	    }
	}
    }

    if ( v >= eq->www-1 ) {
	eq->mask_b = bbb[max_i];
	eq->mask_c = ccc[max_i];
	return;
    }

    optimize_v(eq, bbb[max_i], ccc[max_i], v+1);
}

void _get_tempering_parameter_hard_dc(mt_params *mtp)
{
    int i;
    MaskNode mn0, *cur, *next;
    eqdeg_t eq;

    init_tempering(&eq, mtp);

    for (i=0; i<eq.www; i++)
	eq.gcur_maxlengs[i] = -1;

    mn0.b = mn0.c = mn0.leng = 0;
    mn0.next = NULL;

    cur = &mn0;
    for (i=0; i<LIMIT_V_BEST_OPT; i++) {
	next = optimize_v_hard(&eq, i, cur);
	if (i > 0)
	    delete_MaskNodes(cur);
	cur = next;
    }
    delete_MaskNodes(cur);

    optimize_v(&eq, eq.gmax_b, eq.gmax_c,i);
    mtp->shift0 = eq.shift_0;
    mtp->shift1 = eq.shift_1;
    mtp->shiftB = eq.shift_s;
    mtp->shiftC = eq.shift_t;
    mtp->maskB = eq.mask_b >> eq.ggap;
    mtp->maskC = eq.mask_c >> eq.ggap;
}

static Vector *new_Vector(int nnn)
{
    Vector *v;

    v = (Vector *)malloc( sizeof( Vector ) );
    if( v == NULL ){
        printf("malloc error in \"new_Vector()\"\n");
        exit(1);
    }

    v->cf = (uint32_t *)calloc( nnn, sizeof( uint32_t ) );
    if( v->cf == NULL ){
        printf("calloc error in \"new_Vector()\"\n");
        exit(1);
    }

    v->start = 0;

    return v;
}

static void free_Vector( Vector *v )
{
    if( NULL != v->cf ) free( v->cf );
    if( NULL != v ) free( v );
}

static inline uint32_t trnstmp(eqdeg_t *eq, uint32_t tmp) {
    tmp ^= (tmp >> eq->shift_0) & eq->greal_mask;
    return tmp;
}

static inline uint32_t masktmp(eqdeg_t *eq, uint32_t tmp) {
    tmp ^= (tmp << eq->shift_s) & eq->mask_b;
    tmp ^= (tmp << eq->shift_t) & eq->mask_c;
    return tmp;
}

static inline uint32_t lsb(eqdeg_t *eq, uint32_t x) {
    return (x >> eq->ggap) & 1;
}

static void delete_MaskNodes(MaskNode *head)
{
    MaskNode *t;

    while(head != NULL) {
	t = head->next;
	free(head);
	head = t;
    }
}

static void end_mt_search(prescr_t *pre)
{
    _EndPrescreening_dc(pre);
}

void _EndPrescreening_dc(prescr_t *pre)
{
    int i;

    for (i=0; i<NIRREDPOLY; i++)
      free(pre->modlist[i]);
    free(pre->modlist);
}
