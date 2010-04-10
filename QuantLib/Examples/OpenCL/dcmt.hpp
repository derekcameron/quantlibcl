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

#ifndef dcmt_h
#define dcmt_h

#include "typedefs.h"
#include "mt_params.h"

#define WORDLENGTH 32
#define LSB 0x1
#define NIRREDPOLY 127
#define MAX_IRRED_DEG 9
#define LIMIT_IRRED_DEG 31
#define NOT_REJECTED 1
#define REJECTED 0
#define REDU 0
#define IRRED 1
#define NONREDU 1
#define FOUND 1
#define NOT_FOUND 0
#define DEFAULT_ID_SIZE 16
#define MAX_SEARCH 10000
#define LIMIT_V_BEST_OPT 15

typedef struct {
	int *x;
	int deg;
} dcmt_Polynomial;

typedef struct {
    uint32_t *cf;  /* fraction part */              // status
    int start;     /* beginning of fraction part */ // idx
    int count;	   /* maximum (degree) */
    uint32_t next; /* (bp) rm (shifted&bitmasked) at the maximum degree */
} Vector;

typedef struct mask_node {
    uint32_t b,c;
    int v,leng;
    struct mask_node *next;
} MaskNode;

typedef struct {
    int sizeofA; /* parameter size */
    uint32_t **modlist;
    dcmt_Polynomial **preModPolys;
} prescr_t;

typedef struct {
    uint32_t mt[MT_NN];
    int mti;
} _org_state;

typedef struct {
    uint32_t upper_mask;
    uint32_t lower_mask;
    uint32_t word_mask;
} check32_t;

typedef struct {
    uint32_t bitmask[32];
    uint32_t mask_b;
    uint32_t mask_c;
    uint32_t upper_v_bits;
    int shift_0;
    int shift_1;
    int shift_s;
    int shift_t;
    int mmm;
    int nnn;
    int rrr;
    int www;
    uint32_t aaa[2];
    uint32_t gupper_mask;   /** most significant  (WWW - RRR) bits **/
    uint32_t glower_mask;	/** least significant RRR bits **/
    uint32_t greal_mask;	/** upper WWW bitmask **/
    int ggap; /** difference between machine wordsize and dest wordsize **/
    int gcur_maxlengs[32];	/** for optimize_v_hard **/
    uint32_t gmax_b, gmax_c;
} eqdeg_t;

static mt_params* alloc_mt_params(int n);
static void make_masks(int r, int w, mt_params *mtp);
static dcmt_Polynomial *NewPoly(int degree);
static void FreePoly( dcmt_Polynomial *p);
static dcmt_Polynomial *make_tntm( int n, int m);
static dcmt_Polynomial *PolynomialDup(dcmt_Polynomial *pl);
static dcmt_Polynomial *PolynomialMult(dcmt_Polynomial *p0,dcmt_Polynomial *p1);
static void MakepreModPolys(prescr_t *pre, int mm, int nn, int rr, int ww);
void NextIrredPoly(dcmt_Polynomial *pl, int nth);
static void PolynomialMod( dcmt_Polynomial *wara, const dcmt_Polynomial *waru);
static uint32_t word2bit(dcmt_Polynomial *pl);
static void makemodlist(prescr_t *pre, dcmt_Polynomial *pl, int nPoly);
void _InitPrescreening_dc(prescr_t *pre, int m, int n, int r, int w);
void _InitCheck32_dc(check32_t *ck, int r, int w);
static mt_params *init_mt_search(check32_t *ck, prescr_t *pre, int w, int p);
void _sgenrand_dc(_org_state *st, uint32_t seed);
static void copy_mt_params(mt_params *src, mt_params *dst);
static int get_irred_param(check32_t *ck, prescr_t *pre, _org_state *org, mt_params *mtp, int id, int idw);
mt_params** get_mt_parameters_st(int w, int p, int start_id, int max_id, uint32_t seed, int *count);
static uint32_t nextA(_org_state *org, int w);
uint32_t _genrand_dc(_org_state *st);
static uint32_t nextA_id(_org_state *org, int w, int id, int idw);
int _prescreening_dc(prescr_t *pre, uint32_t aaa);
static int IsReducible(prescr_t *pre, uint32_t aaa, uint32_t *polylist);
int _CheckPeriod_dc(check32_t *ck, _org_state *st, uint32_t a, int m, int n, int r, int w);
void free_mt_params(mt_params *mtp);
void _get_tempering_parameter_dc(mt_params *mtp);
static MaskNode *optimize_v_hard(eqdeg_t *eq, int v, MaskNode *prev_masks);
static int push_stack(eqdeg_t *eq, uint32_t b, uint32_t c, int v, uint32_t *bbb, uint32_t *ccc);
static int push_mask(eqdeg_t *eq, int l, int v, uint32_t b, uint32_t c, uint32_t *bbb, uint32_t *ccc);
static int pivot_reduction(eqdeg_t *eq, int v);
static Vector **make_lattice(eqdeg_t *eq, int v);
static void add(int nnn, Vector *u, Vector *v);
static void free_lattice( Vector **lattice, int v);
static MaskNode *cons_MaskNode(MaskNode *head, uint32_t b, uint32_t c, int leng);
static void next_state(eqdeg_t *eq, Vector *v, int *count);
static int is_zero(int size, Vector *v);
static int calc_pivot(uint32_t v);
static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l);
static void optimize_v(eqdeg_t *eq, uint32_t b, uint32_t c, int v);
void _get_tempering_parameter_hard_dc(mt_params *mtp);
static void free_Vector( Vector *v );
static Vector *new_Vector(int nnn);
static inline uint32_t lsb(eqdeg_t *eq, uint32_t x);
static inline uint32_t masktmp(eqdeg_t *eq, uint32_t tmp);
static inline uint32_t trnstmp(eqdeg_t *eq, uint32_t tmp);
static void delete_MaskNodes(MaskNode *head);
static void end_mt_search(prescr_t *pre);
void _EndPrescreening_dc(prescr_t *pre);

#endif