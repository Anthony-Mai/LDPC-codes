/* TRANSMIT.C - Simulate transmission of bits through a channel. */

/* Copyright (c) 1995-2012 by Radford M. Neal.
 *
 * Permission is granted for anyone to copy, use, modify, and distribute
 * these programs and accompanying documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, and note
 * is made of any changes made to these programs.  These programs and
 * documents are distributed without any warranty, express or implied.
 * As the programs were written for research purposes only, they have not
 * been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own
 * risk.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "channel.h"
#include "open.h"
#include "rand.h"

#include "alloc.h"
#include "blockio.h"
#include "check.h"
#include "dec.h"
#include "mod2convert.h"
#include "mod2dense.h"
#include "mod2sparse.h"
#include "open.h"

#include "rcode.h"

static void usage(char* argv[]);

/* MAIN PROGRAM. */

int main(int argc, char **argv) {
    const char* pchk_file = NULL;
    const char* tfile = NULL; //*rfile = NULL;
    FILE *tf = NULL, *rf = NULL;
    int i, j, it, block_size=0, n_bits=0, nIt=100;
    char junk;
    int seed;
    int cnt;
    int n, b;
    char* dblk, * pchk;
    double* lratio;
    double* bitpr;

    double* awn_data = NULL; /* Places to store channel data */
    int* bsc_data = NULL;

    unsigned iters;            /* Unsigned because can be huge for enum */
    double tot_iter;           /* Double because can be huge for enum */
    double chngd, tot_changed; /* Double because can be fraction if lratio==1*/

    int tot_valid, valid;

    uint64_t inBlock[16];
    float tranBlock[1024];

    if (argc == 1) {
        usage(argv); return 0;
    }

    memset(inBlock, 0, sizeof(inBlock));
    memset(tranBlock, 0, sizeof(tranBlock));

    dec_method = Prprp; // Prprp "prprp" Enum_block "enum-block" Enum_bit "enum-bit"
    max_iter = 20;

    // Look at arguments.  The arguments specifying the channel are looked at by channel_parse in channel.c
    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            if (tfile == nullptr) { tfile = argv[i]; }
            continue;
        }
        if (i + 2 > argc) continue;
        switch (argv[i][1]) {
        case 'i':
        case 'I':
            sscanf(argv[++i], "%d", &max_iter); break;
        case 'n':
        case 'N':
            sscanf(argv[++i], "%d", &nIt); break;
        case 'a':
        case 'A':
            if ((argv[i][3] & 0xDF) == 'L') {
                channel = AWLN;
                sscanf(argv[++i], "%lf", &lwidth);
            } else {
                channel = AWGN;
                sscanf(argv[++i], "%lf", &std_dev);
            }
            break;
        case 'b': // bsc
        case 'B':
            channel = BSC;
            sscanf(argv[++i], "%lf", &error_prob);
            break;
        case 'p':
        case 'P':
            pchk_file = argv[++i]; break;
        case 'm':
        case 'M':
            if (0 == strcasecmp(argv[i + 1], "prprp")) dec_method = Prprp;
            else if (0 == strcasecmp(argv[i + 1], "enum_block")) dec_method = Enum_block;
            else if (0 == strcasecmp(argv[i + 1], "enum_bit")) dec_method = Enum_bit;
            else dec_method = Prprp;
            i++; break;
        case 's':
        case 'S':
            sscanf(argv[++i], "%d", &seed); break;
        }
    }

    //if (!(tfile = argv[1]) || !argv[2] ||
    //    sscanf(argv[2], "%d%c", &seed, &junk) != 1) {
    //    usage();
    //}

    //n = channel_parse(argv + 4, argc - 4);
    //if (n <= 0 /* || argc - 4 - n != 0*/) {
    //    usage();
    //}

    /* See if the source is all zeros or a file. */

    if (sscanf(tfile, "%d%c", &n_bits, &junk) == 1 && n_bits > 0) {
        block_size = 1;
        tf = NULL;
    } else if (sscanf(tfile, "%dx%d%c", &block_size, &n_bits, &junk) == 2 &&
               block_size > 0 && n_bits > 0) {
        n_bits *= block_size;
        tf = NULL;
    } else {
        tf = open_file_std(tfile, "r");
        if (tf == NULL) {
            fprintf(stderr, "Can't open encoded file to transmit: %s\n", tfile);
            exit(1);
        }
    }

    // Read original input file, containing one block.
    for (n_bits = 0; n_bits < sizeof(inBlock) * 8; ) {
        b = getc(tf);
        if (b == '1') {
            inBlock[n_bits >> 6] ^= (1llu << (n_bits & 63));
            n_bits++; continue;
        }

        if (b == '0') {
            n_bits++; continue;
        }
        if (b == EOF) break;
    }
    fclose(tf); tf = NULL;

    // Read parity check file.
    pchk_file = "frame3.pchk"; // TODO: get from command line.

    read_pchk(pchk_file);

    if (N <= M) {
        fprintf(stderr,
            "Number of bits (%d) should be greater than number of checks "
            "(%d)\n",
            N, M);
        exit(1);
    }

    /* Allocate space for data from channel. */

    switch (channel) {
    case BSC: {
        bsc_data = (int*)chk_alloc(N, sizeof * bsc_data);
        break;
    }
    case AWGN:
    case AWLN: {
        awn_data = (double*)chk_alloc(N, sizeof * awn_data);
        break;
    }
    default: {
        abort();
    }
    }

    /* Allocate other space. */

    dblk = (char*)chk_alloc(N, sizeof * dblk);
    lratio = (double*)chk_alloc(N, sizeof * lratio);
    pchk = (char*)chk_alloc(M, sizeof * pchk);
    bitpr = (double*)chk_alloc(N, sizeof * bitpr);

    /* Do the setup for the decoding method. */

    switch (dec_method) {
    case Prprp: {
        prprp_decode_setup();
        break;
    }
    case Enum_block:
    case Enum_bit: {
        enum_decode_setup();
        break;
    }
    default:
        abort();
    }

    // Set random seed to avoid duplications with other programs.
    rand_seed(10 * seed + 3);

    tot_iter = 0;
    tot_valid = 0;
    tot_changed = 0;

    // Create transmission and apply noise
    for (it = 0; it < nIt; it++) {
        switch (channel) {
        case BSC: {
            int32_t* pOut = (int32_t*)tranBlock;
            int bsc_noise;
            for (i = 0; i < n_bits; i++) {
                b = (inBlock[i >> 6] >> (i & 63)) & 1;
                bsc_noise = (rand_uniform() < error_prob) ? 1 : 0;
                pOut[i] = (b ^ bsc_noise);
            }
            break;
        }
        case AWGN: {
            double awgn_noise;
            for (i = 0; i < n_bits; i++) {
                b = (inBlock[i >> 6] >> (i & 63)) & 1;
                awgn_noise = std_dev * rand_gaussian();
                tranBlock[i] = b ? 1.0f + awgn_noise : -1.0f + awgn_noise;
            }
            break;
        }
        case AWLN: {
            double awln_noise;
            for (i = 0; i < n_bits; i++) {
                b = (inBlock[i >> 6] >> (i & 63)) & 1;
                awln_noise = lwidth * rand_logistic();
                tranBlock[i] = b ? 1.0f + awln_noise : -1.0f + awln_noise;
            }
            break;
        }
        default: {
            abort();
        }
        }

        // The transmit data is filled. Now we try to do error correction decoding.

        // Read received blocks, decode, and write decoded blocks.

        /// Read block from received signal data block

        for (i = 0; i < N; i++) {
            int c;
            switch (channel) {
            case BSC: {
                bsc_data[i] = ((int32_t*)tranBlock)[i];
                break;
            }
            case AWGN:
            case AWLN: {
                awn_data[i] = tranBlock[i];
                break;
            }
            default:
                abort();
            }
        }

        /* Find likelihood ratio for each bit. */

        switch (channel) {
        case BSC: {
            for (i = 0; i < N; i++) {
                lratio[i] = bsc_data[i] == 1 ? (1 - error_prob) / error_prob
                    : error_prob / (1 - error_prob);
            }
            break;
        }
        case AWGN: {
            for (i = 0; i < N; i++) {
                lratio[i] = exp(2 * awn_data[i] / (std_dev * std_dev));
            }
            break;
        }
        case AWLN: {
            for (i = 0; i < N; i++) {
                double e, d1, d0;
                e = exp(-(awn_data[i] - 1) / lwidth);
                d1 = 1 / ((1 + e) * (1 + 1 / e));
                e = exp(-(awn_data[i] + 1) / lwidth);
                d0 = 1 / ((1 + e) * (1 + 1 / e));
                lratio[i] = d1 / d0;
            }
            break;
        }
        default:
            abort();
        }

        // Try to decode using the specified method.
        switch (dec_method) {
        case Prprp: {
            iters = prprp_decode(H, lratio, dblk, pchk, bitpr);
            break;
        }
        case Enum_block:
        case Enum_bit: {
            iters = enum_decode(lratio, dblk, bitpr, dec_method == Enum_block);
            break;
        }
        default:
            abort();
        }

        // See if it worked, and how many bits were changed.
        valid = check(H, dblk, pchk) == 0;

        chngd = changed(lratio, dblk, N);

        printf("Try No. %d: block decoded %s, bit changed %d\n", it, valid?"valid":"invalid", (int)chngd);

        tot_iter += iters;
        tot_valid += valid;
        tot_changed += chngd;

        /* Print summary table entry. */

        if (table == 1) {
            printf("%7d %10f    %d  %8.1f\n", block_no, (double)iters, valid,
                (double)chngd);
            /* iters is printed as a double to avoid problems if it's >= 2^31 */
            fflush(stdout);
        }

        /* Write decoded block. */

        //blockio_write(df, dblk, N);


        // Finish a round of decoding
    }

    printf("\nStd_dev=%1.3f. Total valid blocks %d out of %d (%2.3f%%). Total bits changed %d.\n", std_dev, tot_valid, nIt, float(tot_valid)*100.0f/nIt, (int)tot_changed);

    return 0;
}

/* PRINT USAGE MESSAGE AND EXIT. */

void usage(char* argv[]) {
    printf(
        "Usage: %s encode_file [-N rounds] [-Seed ###] [-PCHK pchk_file] [-AWGN/AWLN/BSC #.#] [-prprp/enum_block/enum_bit] [-it ##]\n", argv[0]);
    printf("    -N 1000    Run 1000 rounds of random decoding.\n");
    printf("    -S 1234    Give a pseudo random number seed.\n");
    printf("    -PCHK file.pchk    Provide a *.pchk file.\n");
    printf("    -[AWGN/AWLN/BSC] 0.10    Provide a channel mode with parameter.\n");
    printf("    -[PRPRP/ENUM_BLOCK/ENUM_BIT]    Provide the decode method to use.\n");
    printf("    -IT 20    Specify number of maxium decode iterations.\n");
    channel_usage();
    exit(1);
}
