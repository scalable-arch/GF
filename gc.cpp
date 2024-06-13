#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <cstdio>
#include <iostream>

typedef unsigned long long      gf_binary;
typedef unsigned long long      gf_index;

class GF
{
public:
    GF(unsigned long long primitive_poly)
    {
        this->primitive_poly = primitive_poly;

        degree = -1;
        while (primitive_poly)
        {
            degree++;
            primitive_poly >>= 1;
        }
        this->size           = 1ull<<degree;
    }
    ~GF() {}

    unsigned long long getSize() { return size; }
    int getDegree() { return degree; }

    gf_binary add(gf_binary a, gf_binary b)
    {
        return a^b;
    }

    virtual gf_binary mul(gf_binary a, gf_binary b) = 0;
    virtual gf_binary div(gf_binary a, gf_binary b) = 0;
    virtual gf_binary get_binary(gf_index i) = 0;
    virtual void print() = 0;

    gf_binary pow(gf_binary a, int p)
    {
        gf_binary result = 1;
        for (int i=0; i<p; i++){
            result = mul(result, a);
        }
        return result;
    }

    void hasRoot(GF *gf)
    {
        for (unsigned long long a=1; a<gf->getSize(); a++) {
            gf_binary result = 1;
            for (unsigned long long i=63; i>0; i--)
            {
                if (primitive_poly&(1ull<<i)) {
                    result = gf->add(result, gf->pow(a, i));
                }
            }
            if (result==0) {
                //printf("Found!: %#6llx\n", a);
                std::cout << "Found!: %#6llx" << a << std::endl;
            }
            if ((a%1000000)==999999)
            {
                //printf("Processing %lld\n", a);
                std::cout << "Processing %lld" << a << std::endl;
            }
        }
    }
public:
    int                 degree;
    unsigned long long  size;
    unsigned long long  primitive_poly;
};

class GFSmall : public GF
{
public:
    GFSmall(unsigned long long primitive_poly): GF(primitive_poly)
    {
        log_table       = (gf_index *) malloc(sizeof(gf_index)*size);
        exp_table       = (gf_binary *) malloc(sizeof(gf_binary)*size);

        generate_table();

        if (!check_table())
        {
            printf("FATAL\n");
            exit(1);
        }
    }
    ~GFSmall()
    {
        free(log_table);
        free(exp_table);
    }

    gf_binary mul(gf_binary a, gf_binary b)
    {
        if (a == 0 || b == 0)
        {
            return 0;
        }
        return exp_table[(log_table[a]+log_table[b])%(size-1)];
    }

    gf_binary div(gf_binary a, gf_binary b)
    {
        if (a == 0 || b == 0)
        {
            return 0;
        }
        return exp_table[(log_table[a]-log_table[b])%(size-1)];
    }

    gf_binary get_binary(gf_index i)
    {
        return exp_table[i];
    }

    void print()
    {
        printf("GF(%llu) with ", size);
        for (unsigned long long i=63; i>0; i--)
        {
            if (primitive_poly&(1ull<<i)) {
                printf("x^%lld+", i);
            }
        }
        printf("1\n");

        printf("Index : Binary\n");
        for (unsigned long long i=0; i<size; i++) {
            printf("%4llu:%4llx    ", i, exp_table[i]);
            if (i%4==3) { printf("\n"); }
        }

        printf("Binary : Index\n");
        for (unsigned long long i=0; i<size; i++) {
            printf("%4llx:%4llu    ", i, log_table[i]);
            if (i%4==3) { printf("\n"); }
        }
    }

private:
    void generate_table()
    {
        gf_binary x = 1;
        for (gf_index i = 0; i < (size-1); i++) {   // index
            exp_table[i] = x;
            log_table[x] = i;
            //printf("%-8lld: %llX\n", i, x);
            x = ((x << 1) ^ (x & (1<<(degree-1)) ? primitive_poly: 0)) & (size-1);
        }
        exp_table[size-1] = exp_table[0]; // Circular list for exp_table
    }

    bool check_table()
    {
        char *visited_array;
        bool visited_all = true;

        visited_array = (char *)malloc(size * sizeof(char));
        memset(visited_array, 0, size * sizeof(char));

        for (gf_binary i=0; i<size; i++)
        {
            visited_array[log_table[i]] = 1;
        }

        for (gf_index i=0; i<(size-1); i++)
        {
            if (visited_array[i]!=1) {
                visited_all = false;
                printf("%lld not visited\n", i);
            }
        }

        free(visited_array);
        if (visited_all)
            return true;
        else
            return false;
    }

private:
    gf_binary          *log_table;
    gf_binary          *exp_table;
};

class GFBig : public GF
{
public:
    GFBig(unsigned long long primitive_poly) : GF(primitive_poly) {}
    ~GFBig() {}

    gf_binary mul(gf_binary a, gf_binary b) {
        gf_binary result = 0;
        gf_binary temp = b;

        while (a) {
            if (a & 1ull) {
                result ^= temp;
            }
            temp <<= 1;
            if (temp & (1ull<<degree)) {  // Check if overflow occurs
                temp ^= primitive_poly;
            }
            a >>= 1;
        }

        return result;
    }
    void print()
    {
        printf("GF(%llu) with ", size);
        for (unsigned long long i=63; i>0; i--)
        {
            if (primitive_poly&(1ull<<i)) {
                printf("x^%lld+", i);
            }
        }
        printf("1\n");
    }
};

void generate_syndromes(GF *gf)
{
    gf_binary e;

    printf("1st chip\n");
    for (e=1; e<16; e++)
    {
        printf("%4llx ", gf->mul(e, 1));
    }
    printf("\n");
    printf("2nd chip\n");
    for (e=1; e<16; e++)
    {
        printf("%4llx ", gf->mul(e, 2));
    }
    printf("\n");
    printf("3rd chip\n");
    for (e=1; e<16; e++)
    {
        printf("%4llx ", gf->mul(e, 4));
    }
    printf("\n");
}

bool check_degree4_polynomial(gf_binary a, GF *gf)
{
    gf_binary result = gf->pow(a, 4);
    result = gf->add(result, a);
    result = gf->add(result, 1);
    return (result==0);
}

bool check_degree8_polynomial(gf_binary a, GF *gf)
{
    //GF *gf = new GF(0x11D);   // x^8+x^4+x^3+x^2+1   --> Found D4
    gf_binary result = gf->pow(a, 8);
    result = gf->add(result, gf->pow(a, 4));
    result = gf->add(result, gf->pow(a, 3));
    result = gf->add(result, gf->pow(a, 2));
    result = gf->add(result, 1);
    return (result==0);
}

int main() {
    // Irreducible polynomial list
    // https://www.partow.net/programming/polynomials/index.html

    // Degree 2
    //GF *gf = new GFSmall(0x7);     // x^2+x^1+1

    // Degree 3
    //GF *gf = new GFSmall(0xb);      // x^3+x^1+1

    // Degree 4
    //GF *gf = new GFSmall(0x13);     // x^4+x^1+1

    // Degree 6
    //GF *gf = new GFSmall(0x43);     // x^6+x^1+1
    //GF *gf = new GFSmall(0x67);     // x^6+x^5+x^2+x^1+1
    //GF *gf = new GFSmall(0x6D);     // x^6+x^5+x^3+x^2+1

    // Degree 8
    //GF *gf = new GFSmall(0x11D);    // x^8+x^4+x^3+x^2+1   --> Found D4
    //GF *gf = new GFSmall(0x12B);    // x^8+x^5+x^3+x^1+1   --> Found D4
    //GF *gf = new GFSmall(0x15F);    // x^8+x^6+x^4+x^3+x^2+x^1+1   --> Found D4
    //GF *gf = new GFSmall(0x163);    // x^8+x^6+x^5+x^1+1   --> Found D4
    //GF *gf = new GFSmall(0x165);    // x^8+x^6+x^5+x^2+1   --> Found D4
    //GF *gf = new GFSmall(0x169);    // x^8+x^6+x^5+x^3+1   --> Found D4
    //GF *gf = new GFSmall(0x1E7);    // x^8+x^7+x^6+x^5+x^2+x^1+1   --> Found D4

    // Degree 11
    //GF *gf = new GFSmall(0x805);    // x^11+x^2+1
    //GF *gf = new GFSmall(0x82b);    // x^11+x^5+x^3+x^1+1
    //GF *gf = new GFSmall(0x82d);    // x^11+x^5+x^3+x^2+1
    //GF *gf = new GFSmall(0x863);    // x^11+x^6+x^5+x^1+1
    //GF *gf = new GFSmall(0x88d);    //x^11+x^7+x^3+x^2+1
    //GF *gf = new GFSmall(0x925);    //x^11+x^8+x^5+x^2+1
    //GF *gf = new GFSmall(0x973);    //x^11+x^8+x^6+x^5+x^4+x^1+1
    //GF *gf = new GFSmall(0x97F);    //x^11+x^8+x^6+x^5+x^4+x^3+x^2+x^1+1
    //GF *gf = new GFSmall(0xA13);    //x^11+x^9+x^4+x^1+1
    //GF *gf = new GFSmall(0xB93);    //x^11+x^9+x^8+x^7+x^4+x^1+1
    //GF *gf = new GFSmall(0xC0D);    //x^11+x^10+x^3+x^2+1
    //GF *gf = new GFSmall(0xC9B);    //x^11+x^10+x^7+x^4+x^3+x^1+1
    //GF *gf = new GFSmall(0xDBB);    //x^11+x^10+x^8+x^7+x^5+x^4+x^3+x^1+1
    //GF *gf = new GFSmall(0xF0B);    //x^11+x^10+x^9+x^8+x^3+x^1+1

    // Degree 12
    //GF *gf = new GFSmall(0x1053);   // x^12+x^6+x^4+x^1+1
    //GF *gf = new GFSmall(0x120D);   // x^12+x^9+x^3+x^2+1
    //GF *gf = new GFSmall(0x130F);   // x^12+x^9+x^8+x^3+x^2+x^1+1
    //GF *gf = new GFSmall(0x1745);   // x^12+x^10+x^9+x^8+x^6+x^2+1
    //GF *gf = new GFSmall(0x1775);   // x^12+x^10+x^9+x^8+x^6+x^5+x^4+x^2+1
    //GF *gf = new GFSmall(0x1857);   // x^12+x^11+x^6+x^4+x^2+x^1+1
    //GF *gf = new GFSmall(0x1A2B);   // x^12+x^11+x^9+x^5+x^3+x^1+1
    //GF *gf = new GFSmall(0x1AD1);   // x^12+x^11+x^9+x^7+x^6+x^4+1
    //GF *gf = new GFSmall(0x1AE1);   // x^12+x^11+x^9+x^7+x^6+x^5+1
    //GF *gf = new GFSmall(0x1B91);   // x^12+x^11+x^9+x^8+x^7+x^4+1
    //GF *gf = new GFSmall(0x1BA7);   // x^12+x^11+x^9+x^8+x^7+x^5+x^2+x^1+1
    //GF *gf = new GFSmall(0x1C27);   // x^12+x^11+x^10+x^5+x^2+x^1+1
    //GF *gf = new GFSmall(0x1D5B);   // x^12+x^11+x^10+x^8+x^6+x^4+x^3+x^1+1
    //GF *gf = new GFSmall(0x1FBB);   // x^12+x^11+x^10+x^9+x^8+x^7+x^5+x^4+x^3+x^1+1

    // Degree 13
    //GF *gf = new GFSmall(0x201B);   // x^13+x^4+x^3+x^1+1
    //GF *gf = new GFSmall(0x22BF);   // x^13+x^9+x^7+x^5+x^4+x^3+x^2+x^1+1
    //GF *gf = new GFSmall(0x23A3);   // x^13+x^9+x^8+x^7+x^5+x^1+1
    //GF *gf = new GFSmall(0x26B1);   // x^13+x^10+x^9+x^7+x^5+x^4+1
    //GF *gf = new GFSmall(0x274F);   // x^13+x^10+x^9+x^8+x^6+x^3+x^2+x^1+1
    //GF *gf = new GFSmall(0x2993);   //x^13+x^11+x^8+x^7+x^4+x^1+1
    //GF *gf = new GFSmall(0x2FFF);   //x^13+x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x^1+1
    //GF *gf = new GFSmall(0x3079);   //x^13+x^12+x^6+x^5+x^4+x^3+1
    //GF *gf = new GFSmall(0x31E1);   //x^13+x^12+x^8+x^7+x^6+x^5+1
    //GF *gf = new GFSmall(0x3315);   //x^13+x^12+x^9+x^8+x^4+x^2+1
    //GF *gf = new GFSmall(0x355D);   //x^13+x^12+x^10+x^8+x^6+x^4+x^3+x^2+1
    //GF *gf = new GFSmall(0x3827);   //x^13+x^12+x^11+x^5+x^2+x^1+1
    //GF *gf = new GFSmall(0x39D3);   //x^13+x^12+x^11+x^8+x^7+x^6+x^4+x^1+1
    //GF *gf = new GFSmall(0x3A29);   //x^13+x^12+x^11+x^9+x^5+x^3+1

    // Degree 14
    //GF *gf = new GFSmall(0x4143);   // x^14+x^8+x^6+x^1+1
    //GF *gf = new GFSmall(0x4443);   // x^14+x^10+x^6+x^1+1
    //GF *gf = new GFSmall(0x46DB);   // x^14+x^10+x^9+x^7+x^6+x^4+x^3+x^1+1
    //GF *gf = new GFSmall(0x4843);   // x^14+x^11+x^6+x^1+1
    //GF *gf = new GFSmall(0x4A65);   // x^14+x^11+x^9+x^6+x^5+x^2+1
// x^14+x^12+x^9+x^8+x^7+x^6+x^5+x^4+1
// x^14+x^12+x^11+x^9+x^8+x^7+x^6+x^5+x^3+x^1+1
// x^14+x^12+x^11+x^10+x^9+x^7+x^4+x^3+1
// x^14+x^13+x^6+x^5+x^3+x^1+1
// x^14+x^13+x^10+x^8+x^7+x^5+x^4+x^3+x^2+x^1+1
// x^14+x^13+x^11+x^6+x^5+x^4+x^2+x^1+1
// x^14+x^13+x^11+x^8+x^5+x^3+x^2+x^1+1
// x^14+x^13+x^12+x^11+x^10+x^7+x^6+x^1+1
// x^14+x^13+x^12+x^11+x^10+x^9+x^6+x^5+1

    // Degree 15
    //GF *gf = new GFSmall(0x8003);   // x^15+x^1+1
    //GF *gf = new GFSmall(0x8011);   // x^15+x^4+1
    //GF *gf = new GFSmall(0x8081);   // x^15+x^7+1
// x^15+x^7+x^6+x^3+x^2+x^1+1
// x^15+x^10+x^5+x^1+1
// x^15+x^10+x^5+x^4+1
// x^15+x^10+x^5+x^4+x^2+x^1+1
// x^15+x^10+x^9+x^7+x^5+x^3+1
// x^15+x^10+x^9+x^8+x^5+x^3+1
// x^15+x^11+x^7+x^6+x^2+x^1+1
// x^15+x^12+x^3+x^1+1
// x^15+x^12+x^5+x^4+x^3+x^2+1
// x^15+x^12+x^11+x^8+x^7+x^6+x^4+x^2+1
// x^15+x^14+x^13+x^12+x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+1

    // Degree 16
    //GF *gf = new GFSmall(0x103DD);  // x^16+x^9+x^8+x^7+x^6+x^4+x^3+x^2+1
    //GF *gf = new GFSmall(0x1100B);  // x^16+x^12+x^3+x^1+1
    //GF *gf = new GFSmall(0x11085);  // x^16+x^12+x^7+x^2+1
// x^16+x^13+x^12+x^10+x^9+x^7+x^6+x^1+1
// x^16+x^13+x^12+x^11+x^7+x^6+x^3+x^1+1
// x^16+x^13+x^12+x^11+x^10+x^6+x^2+x^1+1
// x^16+x^14+x^10+x^8+x^3+x^1+1
// x^16+x^14+x^13+x^12+x^6+x^5+x^3+x^2+1
// x^16+x^14+x^13+x^12+x^10+x^7+1
// x^16+x^15+x^10+x^6+x^5+x^3+x^2+x^1+1
// x^16+x^15+x^11+x^9+x^8+x^7+x^5+x^4+x^2+x^1+1
// x^16+x^15+x^11+x^10+x^7+x^6+x^5+x^3+x^2+x^1+1
// x^16+x^15+x^11+x^10+x^9+x^6+x^2+x^1+1
// x^16+x^15+x^11+x^10+x^9+x^8+x^6+x^4+x^2+x^1+1


    // Degree 32
    //GF *gf = new GFBig(0x100400007ull);    //x^32+x^22+x^2+x^1+1
    //GF *gf = new GF(0x10076B553ull);    // x^32+x^22+x^21+x^20+x^18+x^17+x^15+x^13+x^12+x^10+x^8+x^6+x^4+x^1+1
    //GF *gf = new GF(0x1008345E9ull);    // x^32+x^23+x^17+x^16+x^14+x^10+x^8+x^7+x^6+x^5+x^3+1
    //GF *gf = new GF(0x104C11BB7ull);    // x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x^1+1
    //GF *gf = new GF(0x10FC22F86ull);    // x^32+x^27+x^26+x^25+x^24+x^23+x^22+x^17+x^13+x^11+x^10+x^9+x^8+x^7+x^2+x^1+1
    //GF *gf = new GF(0x1100D4E63ull);    // x^32+x^28+x^19+x^18+x^16+x^14+x^11+x^10+x^9+x^6+x^5+x^1+1

    /*
    gf2->hasRoot(gf);
    gf->print();
    gf2->print();
    */

    gf->print();


    int cnt_per_bit[64] = {0};

    for (gf_index i=0; i<gf->getSize(); i+=8)
    {
        gf_binary binary = gf->get_binary(i);;
        printf("a^%-3lld ", i);

        for (int i=gf->getDegree()-1; i>=0; i--)
        {
            if ((binary>>i)&1) {
                printf("1");
                cnt_per_bit[i]++;
            }
            else
            {
                printf("0");
            }
            if (((i%4)==0) && (i!=0)) {
                printf("_");
            }
        }
        printf("\n");
    }

    for (int i=gf->getDegree()-1; i>=0; i--)
    {
        printf("%10d ", cnt_per_bit[i]);
    }
    printf("\n");

    for (int i=gf->getDegree()-1; i>=0; i--)
    {
        if ((cnt_per_bit[i]%2)==1) {
            printf("1");
        }
        else
        {
            printf("0");
        }
        if (((i%4)==0) && (i!=0)) {
            printf("_");
        }
    }
    printf("\n");
/*
    gf_binary received[8];
    gf_binary syndrome[2];

    for (int pos=0; pos<8; pos++)
    {
        for (int i=0; i<8; i++)
        {
            received[i] = 0;
        }

        for (int error=1; error<16; error++)
        {
            received[pos] = error;

            syndrome[0] = 0;
            syndrome[1] = 0;
            for (int j=0; j<8; j++)
            {
                syndrome[0] ^= received[j];
                //syndrome[1] = gf->add(syndrome[1], gf->mul(received[j], 1ull<j));
                syndrome[1] = gf->mul(received[pos], pos+1);
            }

            //printf("%x %llx %llx\n", error, syndrome[0], syndrome[1]);

            assert(syndrome[0]==error);
            gf_binary cpos = gf->div(syndrome[1], syndrome[0]);
            printf("- %x @%d -> %llx %llx -> %lld\n", error, pos, syndrome[0], syndrome[1], cpos);
        }
    }
    */
    /*
    for (int x = 0; x<gf->getSize(); x++)
    {
        gf_binary tmp;

        //printf("%4x %4llx %4x\n", x, gf->pow(x, 4), x);
        tmp = gf->add(gf->pow(x, 4), gf->pow(x, 1));
        tmp = gf->add(tmp, 1);
        if (tmp==0) {
            printf("Found: %x\n", x);
        }
    }
    */

    /*
    for (gf_binary a=1; a<gf->getSize(); a++)
    {
        if (check_degree8_polynomial(a, gf)) {
            printf("Found: %llX\n", a);
        }
    }*/

    //generate_syndromes(gf);

    return 0;
}

