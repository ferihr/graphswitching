/* Copyright (C) 2020--2023, Ferdinand Ihringer
 * 
 * This is program applies WQH switching to
 * a regular graph G of order v at most MAX_VTS
 * with a partition of type p,p,v-2p with
 * p <= MAX_SET.
 * 
 * Arguments: 
 * First argument is v.
 * Second argument is p.
 * 
 * Input: VTS^2 0's and 1's, the adjacency matrix of G.
 * All other characters are ignored.
 * 
 * Output: "n=VTS", followed by the adjacency matrices of all
 * graphs which can be obtained from G via
 * applying WQH switching once.        
 * 
 * This version has no intrinsic hard limit on v.                                                                        */


#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>

#define MAX_VTS 448

#define INTT int32_t
#define BLKS 14 /* we need MAX_VTS <= sizeof(INTT)*8*BLKS */
#define BLKSZ 5 /* ld(sizeof(INTT)*8) */
#define LOWBLK 31 /* has to be 2^BLKS-1 */

#define MAX_SET 8 /* maximal size of C1 */

/* Here
 * 
 * el: vertex of partition
 * co1: el*BLKS
 * co2a: el >> BLKSZ
 * co2b: el & LOWBLK
 */

typedef struct {
        int el;
        int co1;
        int co2a;
        int co2b;
} partition_el;


/* apply_wqh applies WQH switching with a given partition c[2m] to the
 * adjacency matrix am[VTS*BLKS].
 * 
 * good_sums[i]: number if neighbours of i in C1
 * good_sums[i+VTS]: number if neighbours of i in C2
 */

void apply_wqh(int vts, 
               INTT am[MAX_VTS*BLKS], 
               partition_el c[2*MAX_SET], 
               int set,
               int good_sums[2*MAX_VTS]);

/* The following chooses a partition for WQH switching with |C1| = set.
 * It chooses element number cur, then calls itself recursively.
 */

void choose_partition(int vts, 
                      INTT am[MAX_VTS*BLKS], 
                      partition_el c[2*MAX_SET], 
                      int set, int cur);

/* Test global conditions and apply switching.
 */

void apply_partition_global(int vts,
                            INTT am[MAX_VTS*BLKS], 
                            partition_el c[2*MAX_SET],
                            int set);

/* Test if a partition on C1 can still be regular.
 */
int tst_partition_reg11(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur);

int tst_partition_reg12(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur);

int tst_partition_reg21(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur);

int tst_partition_reg22(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur);

/* Functions which calculate the degrees on the induced subgraph.
 * Naming is hopefully self-explanatory. Last argument is for degrees.
 */
void get_degs11(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[]);

void get_degs12(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[]);

void get_degs21(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[]);

void get_degs22(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[]);

/* Little helpers
 */

int is_in_set(partition_el c[2*MAX_SET], int max, int j);


/* The following prints the adjacency matrix.
 */

void print_am(int vts, INTT am[MAX_VTS*BLKS]);



int main(int argc, char *argv[]) 
{
        INTT am[MAX_VTS*BLKS]; 
        int vts, set, co1, co2a, co2b;
        char newc, cont;
        partition_el c[2*MAX_SET];
        
        /* Some error handling.
         */
        if (argc != 3) {
                printf("Two arguments needed.\n");
                return 1;
        }
        
        vts = atoi(argv[1]);
        set = atoi(argv[2]);
        
        if (!vts || !set) {
                printf("Arguments no valid.\n");
                return 1;
        }
        if (vts > MAX_VTS) {
                printf("More vertices than MAX_VTS.\n");
                return 1;
        }
        if (set > MAX_SET) {
                printf("Larger partition than MAX_SET.\n");
                return 1;
        }
        
        /* initialize */
        for (int i = 0; i < BLKS*vts; i++) {
                am[i] = (INTT) 0;
        }
        
        /* Read input */
        for (int i = 0; i < vts; i++) {
                co1 = i*BLKS;
                for (int j = 0; j < vts; j++) {
                        cont = 0;
                        while (!cont) {
                                if (!scanf("%c", &newc)) {
                                        return 1;
                                }
                                /* we already have 0 everywhere */
                                if (newc == '0') {
                                        cont = 1;
                                } else if (newc == '1') {
                                        co2a = (j >> BLKSZ);
                                        co2b = (j & LOWBLK);
                                        am[co1 + co2a] |= ((INTT) 1 << co2b);
                                        cont = 1;
                                }
                        }
                }
        }
        
        /* Print n for nauty */
        printf("n=%d\n", vts);
        
        choose_partition(vts, am, c, set, 0);
        
        return 0;
}

/* Apply WQH switching.
 */

void apply_wqh(int vts, 
               INTT am[MAX_VTS*BLKS], 
               partition_el c[2*MAX_SET], 
               int set,
               int good_sums[2*MAX_VTS])
{
        int co1, co2a, co2b;
        
        for(int i = 0; i < vts; i++) {
                if (is_in_set(c, 2*set-1, i))
                        continue;
                
                if ((good_sums[2*i] == set && good_sums[2*i+1] == 0) 
                        || (good_sums[2*i] == 0 && good_sums[2*i+1] == set)) {
                        
                        co1 = i*BLKS;
                        co2a = (i >> BLKSZ);
                        co2b = (i & LOWBLK);
      
                        for(int j = 0; j < set; j++) {
                                am[co1 + c[j].co2a] ^= ((INTT) 1 << c[j].co2b);
                                am[co1 + c[j+set].co2a] ^= ((INTT) 1 << c[j+set].co2b);
                                am[c[j].co1 + co2a] ^= ((INTT) 1 << co2b);
                                am[c[j+set].co1 + co2a] ^= ((INTT) 1 << co2b);
                        }
                }
        }
        
        return;
}

void choose_partition(int vts, 
                      INTT am[MAX_VTS*BLKS], 
                      partition_el c[2*MAX_SET], 
                      int set, int cur) 
{
        int start = 0;
        
        if (cur == 2*set) {
                apply_partition_global(vts, am, c, set);
                return;
        }
        
        if (cur > 0) {
                start = c[cur-1].el+1;
                if(cur == set)
                        start = c[0].el+1;
        }
        
        /* Try to find a good new tuple. */
        for (c[cur].el = start; c[cur].el < vts; c[cur].el++) {
                if (is_in_set(c, cur-1, c[cur].el) )
                        continue;
                
                c[cur].co1 = c[cur].el*BLKS;
                c[cur].co2a = c[cur].el >> BLKSZ;
                c[cur].co2b = c[cur].el & LOWBLK;
                
                if (cur <= set/2 || cur == set) {
                        /* do nothing */
                } else if (cur < set) {
                        /* If set/2 < cur < set, then we can only
                        * test if regularity is still feasible.
                        */
                        if(!tst_partition_reg11(am, c, set, cur))
                                continue;
                } else {
                        /* Last case, we can test all local conditions. */
                        if (!tst_partition_reg21(am, c, set, cur))
                                continue;
                        if (!tst_partition_reg12(am, c, set, cur))
                                continue;
                        if (!tst_partition_reg22(am, c, set, cur))
                                continue;
                }
                
                choose_partition(vts, am, c, set, cur+1);
        }
}

int tst_partition_reg11(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur)
{
        int r_sum[MAX_SET], min, max;
        
        /* degree of vertices is r_sum */
        get_degs11(am, c, set, cur, r_sum);
        min = r_sum[0];
        max = r_sum[0];
        
        for (int j = 1; j <= cur; j++) {
                if (r_sum[j] < min)
                        min = r_sum[j];
                else if (r_sum[j] > max)
                        max = r_sum[j];
                
                /* To be regular after set, we need max-min <= set-cur-1. */
                if (max-min > set-cur-1)
                        return 0;
        }
        
        return 1;
}

int tst_partition_reg22(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur)
{
        int reg[MAX_SET], r_sum[MAX_SET];
        
        /* Degree of vertices is r_sum.
         * Degree on C1, so also C2, is reg[0].
         */
        get_degs11(am, c, set, set-1, reg);
        get_degs22(am, c, set, cur, r_sum);
        
        for (int j = 0; j <= cur-set; j++) {
                if (reg[0] < r_sum[j] || reg[0] > r_sum[j]+2*set-cur-1)
                        return 0;
        }
        
        return 1;
}

int tst_partition_reg12(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur)
{
        int r_sum[MAX_SET], min, max;
        
        get_degs12(am, c, set, cur, r_sum);
        min = r_sum[0];
        max = r_sum[0];
        
        for (int j = 1; j <= set-1; j++) {
                if (r_sum[j] < min)
                        min = r_sum[j];
                else if (r_sum[j] > max)
                        max = r_sum[j];
                
                if (max-min > 2*set-cur-1)
                        return 0;
        }
        
        return 1;
}

int tst_partition_reg21(INTT am[MAX_VTS*BLKS], 
                        partition_el c[2*MAX_SET], 
                        int set, int cur)
{
        int r_sum[MAX_SET];
        
        get_degs21(am, c, set, cur, r_sum);
        
        for (int j = 1; j <= cur-set; j++) {
                /* As we know C1 completely,
                 * all degrees have to be the same.
                 */
                if (r_sum[j] != r_sum[j-1])
                        return 0;
        }
        
        return 1;
}

void get_degs11(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[])
{
        for (int j = 0; j <= cur && j <= set-1; j++) {
                degs[j] = 0;
                
                for (int i = 0; i <= cur && j <= set-1; i++) {
                        degs[j] += (am[c[j].co1 + c[i].co2a] >> c[i].co2b) & 1;
                }
        }
}


void get_degs22(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[])
{
        for (int j = set; j <= cur; j++) {
                degs[j-set] = 0;
                
                for (int i = set; i <= cur; i++) {
                        degs[j-set] += (am[c[j].co1 + c[i].co2a] >> c[i].co2b) & 1;
                }
        }
}


void get_degs12(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[])
{
        for (int j = 0; j <= set-1; j++) {
                degs[j] = 0;
                
                for (int i = set; i <= cur; i++) {
                        degs[j] += (am[c[j].co1 + c[i].co2a] >> c[i].co2b) & 1;
                }
        }
}


void get_degs21(INTT am[MAX_VTS*BLKS],
                partition_el c[2*MAX_SET],
                int set, int cur, int degs[])
{
        for (int j = set; j <= cur; j++) {
                degs[j-set] = 0;
                
                for (int i = 0; i <= set-1; i++) {
                        degs[j-set] += (am[c[j].co1 + c[i].co2a] >> c[i].co2b) & 1;
                }
        }
}

int is_in_set(partition_el c[2*MAX_SET], int max, int i)
{
       for(int j = 0; j <= max; j++)
               if (i == c[j].el)
                       return 1;
        return 0;
}


void apply_partition_global(int vts,
                            INTT am[MAX_VTS*BLKS], 
                            partition_el c[2*MAX_SET],
                            int set)
{
        int good_sums[2*MAX_VTS];
        int any, co1;
        
        any = 0;
        for (int i = 0; i < vts; i++) {
                if (is_in_set(c, 2*set-1, i))
                        continue;
                
                co1 = i*BLKS;
                good_sums[2*i] = 0;
                good_sums[2*i+1] = 0;
                for (int j = 0; j < set; j++) {
                        good_sums[2*i] += (am[co1 + c[j].co2a] >> c[j].co2b) & 1;
                        good_sums[2*i+1] += (am[co1 + c[j+set].co2a] >> c[j+set].co2b) & 1;
                }
                
                if (good_sums[2*i] == good_sums[2*i+1])
                        continue;
                
                if ((good_sums[2*i] == set && good_sums[2*i+1] == 0)
                        || (good_sums[2*i] == 0 && good_sums[2*i+1] == set)) {
                        
                        if (!any)
                                any = 1;
                        continue;
                }
                
                return;
        }
        
        if (!any)
                return;
        
        apply_wqh(vts, am, c, set, good_sums);
        print_am(vts, am);
        apply_wqh(vts, am, c, set, good_sums);
}


void print_am(int vts, INTT am[MAX_VTS*BLKS])
{
        int co1, co2a, co2b;
        
        for(int i = 0; i < vts; i++) {
                co1 = i*BLKS;
          
                for(int j = 0; j < vts; j++) {
                        co2a = (j >> BLKSZ);
                        co2b = (j & LOWBLK);
            
                        if((am[co1 + co2a] >> co2b) & 1) {
                                printf("1");
                        } else {
                        printf("0");
                        }
                }
                printf("\n");
        }
        printf("\n");
}
