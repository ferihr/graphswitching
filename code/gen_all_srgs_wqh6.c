/* Copyright (C) 2020, Ferdinand Ihringer
 * 
 * This is program applies WQH switching to
 * a regular graph G on VTS vertices with a partition
 * of type 3,3,VTS-6.
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

#define VTS 81

void apply_wqh(char am[VTS][VTS], int c[6], int good_sums[VTS], int good_sums2[VTS]);

int main() {
  char am[VTS][VTS];
  char newc, cont, all, any;
  int r_sum, r_sum2, c[6], good_sums[VTS], good_sums2[VTS];
  
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
  // Partition is {c[0],c[1],c[2]},{c[3],c[4],c[5]}
  for(c[0] = 0; c[0] < VTS; c[0]++) {
    for(c[1] = c[0]+1; c[1] < VTS; c[1]++) {
      for(c[2] = c[1]+1; c[2] < VTS; c[2]++) {
        // First condition (b): induced subgraphs on parts have to be regular
        // First half
        r_sum = 0;
        for(int i = 1; i < 3; i++) {
          r_sum += am[c[0]][c[i]];
        }
        all = 1;
        for(int i = 1; i < 3 && all; i++) {
          r_sum2 = 0;
          for(int j = 0; j < 3; j++) {
            r_sum2 += am[c[j]][c[i]];
          }
          if(r_sum != r_sum2) {
            all = 0;
          }
        }
        if(!all) {
          continue;
        }
        
        for(c[3] = c[0]+1; c[3] < VTS; c[3]++) {
          if(c[3] == c[1] || c[3] == c[2]) {
            continue;
          }
          for(c[4] = c[3]+1; c[4] < VTS; c[4]++) {
            if(c[4] == c[1] || c[4] == c[2]) {
              continue;
            }
            
            // Test if we can still reach condition for second half of condition (b)
            r_sum = 0;
            for(int i = 1; i < 3; i++) {
              r_sum += am[c[0]][c[i]];
            }
            if((am[c[3]][c[4]] > r_sum) || (r_sum > am[c[3]][c[4]]+1)) {
              continue;
            }
            
            // First condition (a): induced subgraph on all 6 has to be regular
            // Test if this is still achievable.
            all = 1;
            r_sum = 0;
            r_sum2 = 0;
            for(int i = 0; i < 3; i++) {
              r_sum += am[c[3]][c[i]];
              r_sum2 += am[c[4]][c[i]];
            }
            if(r_sum != r_sum2) {
              all = 0;
            }
            if(!all) {
              continue;
            }
            
            all = 1;
            r_sum = 0;
            for(int i = 3; i < 5; i++) {
              r_sum += am[c[0]][c[i]];
            }
            for(int i = 1; i < 3 && all; i++) {
              r_sum2 = 0;
              for(int j = 3; j < 5; j++) {
                r_sum2 += am[c[j]][c[i]];
              }
              if(r_sum != r_sum2 && r_sum != r_sum2+1 && r_sum != r_sum2-1) {
                all = 0;
              }
            }
            if(!all) {
              continue;
            }
            
            for(c[5] = c[4]+1; c[5] < VTS; c[5]++) {
              if(c[5] == c[1] || c[5] == c[2]) {
                continue;
              }
              
              // Test first condition (a).
              r_sum = 0;
              for(int i = 0; i < 6; i++) {
                r_sum += am[c[0]][c[i]];
              }
              all = 1;
              for(int i = 1; i < 6 && all; i++) {
                r_sum2 = 0;
                for(int j = 0; j < 6; j++) {
                  r_sum2 += am[c[j]][c[i]];
                }
                if(r_sum != r_sum2) {
                  all = 0;
                }
              }
              if(!all) {
                continue;
              }
              
              // Second half of (b)
              r_sum = 0;
              for(int i = 0; i < 3; i++) {
                r_sum += am[c[0]][c[i]];
              }
              all = 1;
              for(int i = 3; i < 6 && all; i++) {
                r_sum2 = 0;
                for(int j = 3; j < 6; j++) {
                  r_sum2 += am[c[j]][c[i]];
                }
                if(r_sum != r_sum2) {
                  all = 0;
                }
              }
              if(!all) {
                continue;
              }
              
              // Second condition: all are adjacent to both parts GM set equally
              // or the whole of one half
              // Calculate all neighbourhoods
              all = 1;
              for(int i = 0; i < VTS && all; i++) {
                good_sums[i] = am[i][c[0]] + am[i][c[1]] + am[i][c[2]];
                good_sums2[i] = am[i][c[3]] + am[i][c[4]] + am[i][c[5]];
                
                // Test condition
                if(i == c[0] || i == c[1] || i == c[2] || i == c[3] || i == c[4] || i == c[5]) {
                  continue;
                }
                if(good_sums[i] == good_sums2[i] || (good_sums[i] == 3 && good_sums2[i] == 0) 
                    || (good_sums[i] == 0 && good_sums2[i] == 3)) {
                  continue;
                }
                all = 0;
              }
              
              if(!all) {
                continue;
              }
              
              // Last condition: one vertex should be adjacenct to one whole half and none of the other
              // Otherwise, new graph is old graph.
              any = 0;
              for(int i = 0; i < VTS && (!any); i++) {
                if(i == c[0] || i == c[1] || i == c[2] || i == c[3] || i == c[4] || i == c[5]) {
                  continue;
                }
                if(good_sums[i] == 3 && good_sums2[i] == 0) {
                  any = 1;
                }
                if(good_sums2[i] == 3 && good_sums[i] == 0) {
                  any = 1;
                }
              }
              if(!any) {
                continue;
              }
              
              // Now we are good
              // Apply switching, print, apply switching back
              apply_wqh(am, c, good_sums, good_sums2);
              // Print New Graph
              for(int i = 0; i < VTS; i++) {
                for(int j = 0; j < VTS; j++) {
                  printf("%d", am[i][j]);
                }
                printf("\n");
              }
              printf("\n");
              // Revert back
              apply_wqh(am, c, good_sums, good_sums2);
            }
          }
        }
      }
    }
  }
  
  return 0;
}

// Apply GM switching, in place
void apply_wqh(char am[VTS][VTS], int c[6], int good_sums[VTS], int good_sums2[VTS]) {
  for(int i = 0; i < VTS; i++) {
    if(i == c[0] || i == c[1] || i == c[2] || i == c[3] || i == c[4] || i == c[5]) {
      continue;
    }
    
    if((good_sums[i] == 3 && good_sums2[i] == 0) || (good_sums[i] == 0 && good_sums2[i] == 3)) {
      for(int j = 0; j < 3; j++) {
          am[i][c[j]] ^= (char) 1;
          am[i][c[j+3]] ^= (char) 1;
          am[c[j]][i] ^= (char) 1;
          am[c[j+3]][i] ^= (char) 1;
      }
    }
  }
  
  return;
}



