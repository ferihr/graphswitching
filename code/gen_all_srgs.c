/* Copyright (C) 2020, Ferdinand Ihringer
 * 
 * This is program applies Godsil-McKay switching to
 * a regular graph G on VTS vertices with a partition
 * of type 4,VTS-4.
 * 
 * Input: VTS^2 0's and 1's, the adjacency matrix of G.
 *        All other characters are ignored.
 * 
 * Output: "n=VTS", followed by the adjacency matrices of all
 *         graphs which can be obtained from G via
 *         applying GM switching once.
 * 
 * Note that this program is primarily for my own use.
 * 
 * This version has no hard limit on VTS.                    */


#include <stdio.h>

#define VTS 63

void apply_gm(char am[VTS][VTS], int c[4], int good_sums[VTS]);

int main() {
  char am[VTS][VTS];
  char newc, cont, all, any;
  int r_sum, c[4], good_sums[VTS];
  
  // Read input
  for(int i = 0; i < VTS; i++) {
    for(int j = 0; j < VTS; j++) {
      cont = 0;
      while(!cont) {
        if(!scanf("%c", &newc)) {
          return 1;
        }
        if(newc == '0') {
          am[i][j] = 0;
          cont = 1;
        }
        else if(newc == '1') {
          am[i][j] = 1;
          cont = 1;
        }
      }
    }
  }
  
  // Print n for nauty
  printf("n=%d\n", VTS);
  
  // All combinations
  for(c[0] = 0; c[0] < VTS; c[0]++) {
    for(c[1] = c[0]+1; c[1] < VTS; c[1]++) {
      for(c[2] = c[1]+1; c[2] < VTS; c[2]++) {
        for(c[3] = c[2]+1; c[3] < VTS; c[3]++) {
          // First condition: induced subgraph has to be regular
          r_sum = am[c[0]][c[1]] + am[c[0]][c[2]] + am[c[0]][c[3]];
          all = 1;
          for(int i = 1; i < 4 && all; i++) {
            if(r_sum != am[c[i]][c[0]] + am[c[i]][c[1]] + am[c[i]][c[2]] + am[c[i]][c[3]]) {
              all = 0;
            }
          }
          if(!all) {
            continue;
          }
          
          // Second condition: all vertices have 0, 2, or 4 neighbours in GM set 
          // Calculate all neighbourhoods
          for(int i = 0; i < VTS; i++) {
            good_sums[i] = am[i][c[0]] + am[i][c[1]] + am[i][c[2]] + am[i][c[3]];
          }
          // Test condition
          all = 1;
          for(int i = 0; i < VTS && all; i++) {
            if(i == c[0] || i == c[1] || i == c[2] || i == c[3]) {
              continue;
            }
            if(good_sums[i] == 0 || good_sums[i] == 2 || good_sums[i] == 4) {
              continue;
            }
            all = 0;
          }
          if(!all) {
            continue;
          }
          
          // Last condition: one vertex should have 2 neighbours
          // Otherwise, new graph is old graph.
          any = 0;
          for(int i = 0; i < VTS && (!any); i++) {
            if(i == c[0] || i == c[1] || i == c[2] || i == c[3]) {
              continue;
            }
            if(good_sums[i] == 2) {
              any = 1;
            }
          }
          if(!any) {
            continue;
          }
          
          // Now we are good
          // Apply switching, print, apply switching back
          apply_gm(am, c, good_sums);
          // Print New Graph
          for(int i = 0; i < VTS; i++) {
            for(int j = 0; j < VTS; j++) {
              printf("%d", am[i][j]);
            }
            printf("\n");
          }
          printf("\n");
          // Revert back
          apply_gm(am, c, good_sums);
        }
      }
    }
  }
  
  return 0;
}

// Apply GM switching, in place
void apply_gm(char am[VTS][VTS], int c[4], int good_sums[VTS]) {
  for(int i = 0; i < VTS; i++) {
    if(i == c[0] || i == c[1] || i == c[2] || i == c[3]) {
      continue;
    }
    
    if(good_sums[i] == 2) {
      for(int j = 0; j < 4; j++) {
        am[i][c[j]] ^= (char) 1;
        am[c[j]][i] ^= (char) 1;
      }
    }
  }
  
  return;
}



