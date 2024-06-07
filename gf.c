#include <stdio.h>
#include <stdlib.h>

//#define GF_DEGREE	2
//#define PRIMITIVE_POLY	0x7	// x^2 + x + 1

//#define GF_DEGREE	3
//#define PRIMITIVE_POLY	0xb	// x^3 + x + 1

//#define GF_DEGREE	4
//#define PRIMITIVE_POLY	0x13	// x^4 + x + 1

#define GF_DEGREE	6
#define PRIMITIVE_POLY	0x43	// x^6 + x + 1
//#define PRIMITIVE_POLY	0x67	// x^6 + x^5 + x^2 + x + 1
//#define PRIMITIVE_POLY	0x6D	// x^6 + x^5 + x^3 + x^2 + 1

//#define GF_DEGREE       8
//#define PRIMITIVE_POLY 	0x11D 	// x^8 + x^4 + x^3 + x^2 + 1
//#define PRIMITIVE_POLY	0x12B	// x^8 + x^5 + x^3 + x^1 + 1

#define GF_SIZE         (1<<GF_DEGREE)
typedef unsigned long long      gf_binary;
typedef unsigned long long      gf_index;

gf_binary log_table[GF_SIZE] = {0};
gf_index exp_table[GF_SIZE] = {0};

gf_binary gf_add(gf_binary a, gf_binary b)
{
    return a ^ b;
}

gf_binary gf_mul(gf_binary a, gf_binary b)
{
    if (a == 0 || b == 0) return 0;
    return exp_table[(log_table[a] + log_table[b])%(GF_SIZE-1)];
}

gf_binary gf_pow(gf_binary a, int p)
{
	gf_binary result = 1;
	for (int i=0; i<p; i++)
	{
		result = gf_mul(result, a);
	}
	return result;
}

unsigned char gf_inv(unsigned char a)
{
    if (a == 0) return 0;
    return exp_table[(GF_SIZE-1) - log_table[a]];
}

void generate_gf_tables()
{
    gf_binary x = 1;
    for (unsigned long long i = 0; i < (GF_SIZE-1); i++) {   // index
        exp_table[i] = x;
        log_table[x] = i;
	printf("%-8lld: %llX\n", i, x);
        x = ((x << 1) ^ (x & (1<<(GF_DEGREE-1)) ? PRIMITIVE_POLY : 0)) & (GF_SIZE-1);
    }
    exp_table[GF_SIZE-1] = exp_table[0]; // Circular list for exp_table
}

void print_gf_table(unsigned long long *table, const char *name) {
    printf("%s:\n", name);
    for (int i = 0; i < GF_SIZE; i++) {
        printf("%08llX ", table[i]);
	printf("\n");
        //if ((i + 1) % 16 == 0) printf("\n");
    }
    printf("\n");
}

int main() {
    generate_gf_tables(log_table, exp_table);

    print_gf_table(log_table, "Log Table");
    //print_gf_table(exp_table, "Exp Table");

    //gf_binary a = 0x53;
    //gf_binary b = 0xCA;

    //gf_binary sum = gf_add(a, b);
    //gf_binary product = gf_mul(a, b);
    //gf_binary inverse = gf_inv(a);

    //printf("Addition: %08llX + %08llX = %08llX\n", a, b, sum);
    //printf("Multiplication: %08llX * %08llX = %08llX\n", a, b, product);
    //printf("Inverse: %08llX^-1 = %08llX\n", a, inverse);
    //printf("Power: = %08llX\n", gf_pow(1, 1));
    //printf("Power: = %08llX\n", gf_pow(1, 2));
    //printf("Power: = %08llX\n", gf_pow(2, 1));
    //printf("Power: = %08llX\n", gf_pow(2, 2));

    // 1st symbol
    printf("First data (index / binary) \n");
    for (gf_binary i=1; i<16; i++) {
	    printf("%3lld / %4llx\n", log_table[i], i);
    }
    // 2nd symbol
    printf("Second data (index / binary) \n");
    for (gf_binary i=1; i<16; i++) {
	    unsigned long long temp = (i&0x7) + 8 + ((i&0x8)<<2);
	    temp = gf_mul(temp, 2);
	    printf("%3lld / %4llx\n", log_table[temp], temp);
    }
    // 3rd symbol
    printf("Third data (index / binary) \n");
    for (gf_binary i=1; i<16; i++) {
	    unsigned long long temp = (i&0x3) + 8 + ((i&0xC)<<2);
	    temp = gf_mul(temp, 4);
	    printf("%3lld / %4llx\n", log_table[temp], temp);
    }
    // 4th symbol
    printf("Fourth data (index / binary) \n");
    for (gf_binary i=1; i<16; i++) {
	    unsigned long long temp = (i&0x1) + 6 + ((i&0xe)<<2);
	    temp = gf_mul(temp, 8);
	    printf("%3lld / %4llx\n", log_table[temp], temp);
    }
    /*
    for (gf_binary i=0; i < (GF_SIZE-1); i++) {
	    //gf_binary result = gf_pow(i, 4) ^ i ^ 1;
	    gf_binary result = gf_pow(i, 2) ^ i ^ 1;

	    //printf("Binary %lld / Index %lld : %llx ^ %llx ^ 1 = %llx\n", i, log_table[i], gf_pow(i, 4), i, result);
	    if (result==0) {
		    printf("Found: %llX: %llX\n", i, log_table[i]);
	    }
    }
    */

    return 0;
}
