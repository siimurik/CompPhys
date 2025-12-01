/* pointer.c */
/*
	Simple illustration of the action of pointers
*/

#include <stdio.h>

int main() {

	int u = 5;
	int v;
	unsigned int *pu; 	// Declare pointer to an integer variable 
	unsigned int *pv; 	// Declare pointer to an integer variable 
	
	pu = &u; 	// Assign address of u to pu
	v = *pu; 	// Assign value of u to v
	pv = &v; 	// Assign address of v to pv

	printf("\nu = %d  &u = %p  pu = %p  *pu = %d", u, &u, pu, *pu);
  	printf("\nv = %d  &v = %p  pv = %p  *pv = %d\n", v, &v, pv, *pv);
	
	return 0;
}
