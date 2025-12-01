#include <stdio.h>
#include <ctype.h>

int main() {

	char unit;
	float temp;

	printf("\nIs the temperature in (F) or (C)?: ");
	scanf("%c", &unit);
	
	unit = toupper(unit); 	// Capitalizes a lower case
							// letter. Ex: c -> C
	if (unit == 'C') {
		printf("The temp is currently in C.\n");
		scanf("%f", &temp);
		temp = (temp * 9./5.) + 32.0;
		printf("The temp in Farenheit is: %.1f\n", temp);
	}
	else if (unit == 'F') {
		printf("The temp is currently in F.\n");
		scanf("%f", &temp);
		temp = (temp - 32.0) * 5./9.;
		printf("The temp in Celsius is: %.1f\n", temp);
	}
	else{
		printf("%c is nor a valid unit.\n", unit);
	}

	return 0;
}
