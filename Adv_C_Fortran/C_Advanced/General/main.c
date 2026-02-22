#include <stdio.h>   // Standard library - in and out
#include <stdint.h>  // Standard library - integer pack
#include <stdbool.h> // Standard library - bool pack
#include <stdlib.h>  // For malloc

#include "foo.h"

typedef int8_t  i8;
typedef int16_t i16;
typedef int32_t i32; // -2e9 < i32 < 2e9
typedef int64_t i64;

// Only positive integer values
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32; // for things that are less than 4 billion // 0 < i32 < 4e9
typedef uint64_t u64; // typically meant for dealing with memory

typedef i8 b8;
typedef i32 b32;

typedef float  f32;
typedef double f64;

#define MIN(a, b) (a < b ? a : b)

typedef struct // typedef is needed so we do not need to 
{              // add struct in front of vec2f in main().
    f32 x;     // by default it would be struct vec2f{}.
    f32 y;
} vec2f;

// Function has to come above where you are using it
//void foo();

int main(void) {

    i32 x = 1;
    i32 y = 123;

    printf("%d\n", MIN(x, y));

    /*
    u32 num_vectors = 1000;
    vec2f* vectors = (vec2f*)malloc(sizeof(vec2f) * num_vectors);

    foo();

    for (u32 i = 0; i < num_vectors; i++){
        vectors[i].x = i * 1;
        vectors[i].x = i * 2;
    }

    free(vectors);
    */

    /*
    // 'nums' is actually just a pointer. Hence you cannot
    // just get its length with len() or sth equiv. 
    i32 nums[] = { 1, 2, 3, 4, 5 };

    //nums[123123] = 123; // will fail but we can find it easlitly if we compile with
    // gcc main.c -fsanitize=address // on Unix systems

    for (u32 i = 0; i < 5; i++){
        printf("%d ", nums[i]);
    }
    printf("\n");

    for (u32 i = 0; i < 5; i++){
        //*(nums + i) = 5 - i;
        nums[i] = 5 - i; // simpler equivalent
        // funnily enough, this means that we can also use
        // i[nums] = 5 - i to get the same thing. 
        // obv don't write code like this.
    }

    for (u32 i = 0; i < 5; i++){
        printf("%d ", nums[i]);
    }
    printf("\n");

    *nums = 6; // modifies the first number
    *(nums + 1) = 8; // '+1' offsets it by one (modifies the next value in nums)
                     // aka, increments the pointer 4 bytes in memory

    for (u32 i = 0; i < 5; i++){
        printf("%d ", nums[i]);
    }
    printf("\n");
    */

    /*
    vec2f v = { 1, 2 };
    vec2f* pv = &v;

    printf("%f %f\n", v.x, v.y);

    // How it looks under the hood
    // (*pv).x = 3;
    // (*pv).y = 2;
    pv->x = 3;  // Easier on the eyes
    pv->y = 2;

    printf("%f %f\n", v.x, v.y);
    */

    /*
    i32 x = 123;
    // The asterix shows that 'px' is a pointer
    i32* px = &x; // px corresponds to the address of x

    printf("%d %p\n", x, px);

    // Pointer dereferencing
    // Asterix says:
    // "Go to wherever this place is in memory 
    // and set its value to ..."
    *px = 321;
    printf("%d %p\n", x, px);
    */

    /*
    vec2f v = {
        .x = 2, // anything not specified here gets set to zero 
        .y = 1, 
    };

    printf("Vector = < %f, %f >\n", v.x, v.y);
    
    v = (vec2f){
        .x = 3, .y = 4,
    };
    
    printf("Vector = < %f, %f >\n", v.x, v.y);
    */
    return 0;
}

// void foo() {
//     printf("hello\n");
// }