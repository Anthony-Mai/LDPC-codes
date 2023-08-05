// This program is used to generate the new LDPC coding matrix using in the GNSS Project.
// To build:
//    g++ -c -O     genmax.cpp
//    g++           genmax.o -lm -o genmax

#include <stdio.h>
#include <stdint.h>
#include <string.h>

void Usage(int argc, char* argv[]) {
    printf("Usage: %s -[S##/G##] myfile.pchk\n", argv[0]);
    printf("Examples:\n");
    printf("    %s -S272 myfile.pchk\n        Generate the Star 272 encoding matrix\n", argv[0]);
    printf("Please refer to this link for reference of *.pchk files:\n");
    printf("    https://glizen.com/radfordneal/ftp/LDPC-2012-02-11/index.html\n");
    printf("Once a valida PCHK file is generated you can do the following:\n");
    printf("Generate a *.gen generation matrix for encoding a frame:\n");
    printf("    ./make-gen star272b.pchk star272b.gen sparse\n");
    printf("Encode raw input data frame into encoded frame (origin272.txt is a text file of 0 and 1, one char for one bit.):\n");
    printf("    ./encode -f star272b.pchk star272b.gen origin272.txt star272b.bin\n");
    printf("Run multiple iterations of decoding for benchmarking performance:\n");
    printf("    ./mbench star272b.bin -N 1000 -S 1234579 -PCHK star272b.pchk -AWGN 0.60 -M prprp -it 100\n");
    printf("Explanation of the command line above:\n");
    printf("    star272b.bin (encoded block to be decoded)\n");
    printf("    -N 1000 (Run 1000 rounds of decoding)\n");
    printf("    -S 1234579 (Use this seed to generate pseudo random numbers)\n");
    printf("    -PCHK star272b.pchk (The parity check file needed for decoding)\n");
    printf("    -AWGN 0.60 (AWGN noise level 0.60. Es/N0 = 20*log10(sqrt(0.5)/0.60) = 1.426675 dB\n");
    printf("    -M prprp (Decode method is probability propagation. We always use this)\n");
    printf("    -it 100 (Decode iteration limited to 100 steps\n");
}

union u320 {
#define popcnt64 __builtin_popcount
    uint64_t d[5];
    void reset() { d[0]=d[1]=d[2]=d[3]=d[4]=0llu; }
    void set(uint32_t idx) { d[idx>>6] |= (1llu<<(idx&63)); }
    void eor(uint32_t idx) { d[idx>>6] ^= (1llu<<(idx&63)); }
    uint32_t popcnt() const { return (popcnt64(d[0])+popcnt64(d[1])+popcnt64(d[2])+popcnt64(d[3])+popcnt64(d[4])); }
    u320 operator & (const u320& s) const {
        return u320{d[0]&s.d[0], d[1]&s.d[1], d[2]&s.d[2], d[3]&s.d[3], d[4]&s.d[4]};
    }
    u320& operator ^= (const u320& s) {
        d[0]^=s.d[0]; d[1]^=s.d[1]; d[2]^=s.d[2]; d[3]^=s.d[3]; d[4]^=s.d[4];
        return (*this);
    }
};

void GenStar272(const char* fileName);
void SavePchk(const char* fname, u320 (&D)[548]);

int main(int argc, char* argv[]) {
    printf("Genesis Utility Program to Generate New LDPC coding matrixes\n");
    if (argc <= 1) {
        Usage(argc, argv); return 0;
    }
    const char* outFileName = nullptr;
    if (argv[1][0] == '-') switch (argv[1][1]) {
    case 's':
    case 'S':
        if (argc >= 3) outFileName = argv[2];
        GenStar272(outFileName);
        break;
    default:
        break;
    }
    return 0;
}

void SavePchk(const char* fname, u320 (&D)[548]) {
    int32_t i,j,N=272, M=544, x, y;
    uint8_t c1='P', c0=0x80, cc=0x00;
    FILE* fout = fopen(fname, "wb");
    fwrite(&c0, 1, 1, fout); fwrite(&c1, 1, 1, fout); fwrite(&cc, 1, 1, fout); fwrite(&cc, 1, 1, fout);
    fwrite(&N, sizeof(N), 1, fout); fwrite(&M, sizeof(M), 1, fout);
    for (j=0; j<N; j++) {
        y = 0 - j - 1;
        fwrite(&y, sizeof(y), 1, fout);
        for (i=0; i<M; i++) {
            if (i>=N) {
                if (i == j+N) {
                    x = i+1;
                    fwrite(&x, sizeof(x), 1, fout);
                }
                continue;
            }
            if (D[j+N].d[i>>6]&(1llu<<(i&63))) {
                x = i+1;
                fwrite(&x, sizeof(x), 1, fout);
            }
        }
    }
    i=0;  fwrite(&i, sizeof(i), 1, fout);
    fclose(fout);
}


// Validate extended Star-Codes for a 17x16 MDS coding matrix. Previous one may be wrong.
// LDPC-codes$ ./make-gen star272b.pchk star272b.gen sparse
// Number of 1s per check in L is 1.0, U is 1.0, B is 31.1, total is 33.1
// LDPC-codes$ ./encode -f star272b.pchk star272b.gen origin272.txt star272b.bin
// Encoded 1 blocks, source block size 272, encoded block size 544
// LDPC-codes$ ./mbench star272b.bin -N 1000 -S 1234579 -PCHK star272b.pchk -AWGN 0.60 -M prprp -it 100
void GenStar272(const char* fileName) {
    uint32_t i,j,k,x,y,s;
    u320 D[548];
    u320 S[18];

    printf("Will generate encoding matric file %s.\n", fileName);

    memset(S, 0, sizeof(S));
    memset(D, 0, sizeof(D));
    S[0].reset();
    for (i=0; i<17; i++) {
        for (j=0; j<16; j++) {
            k = i*16+j;
            D[k].eor(k);
            //S[i].eor(k);
        }
    }
    for (i=0; i<17; i++) {
        s = i; // Slope.
        S[0].reset();
        for (j=0; j<17; j++) {
            u320 t; t.reset();
            for (x=0,y=j; x<17; x++,y+=s) {
                if (y>=17) { y -= 17;}
                if (y<16) t.eor(x*16+y);
                else t ^= S[x];
            }
            if (j<16) {
                D[(i+17)*16+j] = t;
            } else {
                for (k=0; k<16; k++) {
                    D[(i+17)*16+k] ^= t;
                }
            }
        }
    }

    printf("S=\n");
    for (k=0; k<17; k++) {
        printf("[%03d] %016llX %016llX %016llX %016llX %016llX\n", k, S[k].d[4], S[k].d[3], S[k].d[2], S[k].d[1], S[k].d[0]);
    }
    printf("\n");
    printf("Result:\n");
    for (i=0; i<34; i++) {
        for (j=0; j<16; j++) {
            k = i*16+j;
            printf("[%03d] %016llX %016llX %016llX %016llX %016llX\n", k, D[k].d[4], D[k].d[3], D[k].d[2], D[k].d[1], D[k].d[0]);
        }
        printf("\n");
    }

    if (1) {
        SavePchk(fileName, D);
        printf("PCHK file saved to %s\n", fileName);
    }

    printf("All Done\n");
}

